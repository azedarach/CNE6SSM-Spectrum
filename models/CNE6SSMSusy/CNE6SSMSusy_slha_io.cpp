// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Sun 19 Apr 2015 20:31:38

#include "CNE6SSMSusy_slha_io.hpp"
#include "CNE6SSMSusy_two_scale_input_parameters.hpp"
#include "logger.hpp"
#include "wrappers.hpp"
#include "numerics2.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "config.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/bind.hpp>

using namespace softsusy;

namespace flexiblesusy {

char const * const CNE6SSMSusy_slha_io::drbar_blocks[NUMBER_OF_DRBAR_BLOCKS] =
   { "gauge", "Yu", "Yd", "Ye", "HMIX", "ESIXHEYUK", "ESIXRUN", "ESIXGDYUK",
     "ESIXFUYUK", "ESIXFDYUK", "ESIXKAPPA", "ESIXLAMBDA"
   }
;

CNE6SSMSusy_slha_io::CNE6SSMSusy_slha_io()
   : slha_io()
{
}

void CNE6SSMSusy_slha_io::clear()
{
   slha_io.clear();
}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CNE6SSMSusy_slha_io::set_extpar(const CNE6SSMSusy_input_parameters<Two_scale>& input)
{
   std::ostringstream extpar;

   extpar << "Block EXTPAR\n";
   extpar << FORMAT_ELEMENT(65, input.sInput, "sInput");
   extpar << FORMAT_ELEMENT(66, input.TanTheta, "TanTheta");
   extpar << FORMAT_ELEMENT(67, input.vphiInput, "vphiInput");
   extpar << FORMAT_ELEMENT(72, input.QSInput, "QSInput");
   slha_io.set_block(extpar);

}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CNE6SSMSusy_slha_io::set_minpar(const CNE6SSMSusy_input_parameters<Two_scale>& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(3, input.TanBeta, "TanBeta");
   slha_io.set_block(minpar);

}

/**
 * Stores the SMINPUTS input parameters in the SLHA object.
 *
 * @param qedqcd class of Standard Model parameters
 */
void CNE6SSMSusy_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void CNE6SSMSusy_slha_io::set_spinfo(const Problems<CNE6SSMSusy_info::NUMBER_OF_PARTICLES>& problems)
{
   std::ostringstream spinfo;
   spinfo << "Block SPINFO\n"
          << FORMAT_SPINFO(1, PKGNAME)
          << FORMAT_SPINFO(2, FLEXIBLESUSY_VERSION);

   if (problems.have_warning()) {
      std::ostringstream warnings;
      problems.print_warnings(warnings);
      spinfo << FORMAT_SPINFO(3, warnings.str());
   }

   if (problems.have_problem()) {
      std::ostringstream problems_str;
      problems.print_problems(problems_str);
      spinfo << FORMAT_SPINFO(4, problems_str.str());
   }

   slha_io.set_block(spinfo, SLHA_io::front);
}

/**
 * Write SLHA object to file.
 *
 * @param file_name file name
 */
void CNE6SSMSusy_slha_io::write_to_file(const std::string& file_name)
{
   slha_io.write_to_file(file_name);
}

/**
 * Read (DR-bar) model parameter output scale from MODSEL entry 12
 */
double CNE6SSMSusy_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

/**
 * Read SLHA object from file
 *
 * @param file_name file name
 */
void CNE6SSMSusy_slha_io::read_from_file(const std::string& file_name)
{
   slha_io.read_from_file(file_name);
   slha_io.read_modsel();
}

/**
 * Fill struct of model input parameters from SLHA object (MINPAR and
 * EXTPAR blocks)
 *
 * @param input struct of model input parameters
 */
void CNE6SSMSusy_slha_io::fill(CNE6SSMSusy_input_parameters<Two_scale>& input) const
{
   // part of the reason the workaround is not ideal - need
   // to disambiguate between which overload is used...
   void (*fill_two_scale_minpar_tuple) (CNE6SSMSusy_input_parameters<Two_scale>&,
                                        int, double) = &CNE6SSMSusy_slha_io::fill_minpar_tuple;
   void (*fill_two_scale_extpar_tuple) (CNE6SSMSusy_input_parameters<Two_scale>&,
                                        int, double) = &CNE6SSMSusy_slha_io::fill_extpar_tuple;
   SLHA_io::Tuple_processor minpar_processor
      = boost::bind(fill_two_scale_minpar_tuple, boost::ref(input), _1, _2);
   SLHA_io::Tuple_processor extpar_processor
      = boost::bind(fill_two_scale_extpar_tuple, boost::ref(input), _1, _2);

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);

   input.MuPhiInput = slha_io.read_entry("HMIXIN", 31);
   input.KappaPrInput = slha_io.read_entry("HMIXIN", 32);
   input.SigmaxInput = slha_io.read_entry("HMIXIN", 33);
   slha_io.read_block("ESIXHEYUKIN", input.hEInput);
   input.SigmaLInput = slha_io.read_entry("ESIXRUNIN", 42);
   slha_io.read_block("ESIXGDYUKIN", input.gDInput);
   slha_io.read_block("ESIXFUYUKIN", input.fuInput);
   slha_io.read_block("ESIXFDYUKIN", input.fdInput);
   slha_io.read_block("ESIXKAPPAIN", input.KappaInput);
   slha_io.read_block("ESIXLAMBDAIN", input.Lambda12Input);
   input.MuPrInput = slha_io.read_entry("ESIXRUNIN", 0);
   input.LambdaxInput = slha_io.read_entry("ESIXRUNIN", 1);
   input.XiFInput = slha_io.read_entry("HMIXIN", 30);
}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings
 */
void CNE6SSMSusy_slha_io::fill(Spectrum_generator_settings& settings) const
{
   SLHA_io::Tuple_processor flexiblesusy_processor
      = boost::bind(&CNE6SSMSusy_slha_io::fill_flexiblesusy_tuple, boost::ref(settings), _1, _2);

   slha_io.read_block("FlexibleSUSY", flexiblesusy_processor);
}

void CNE6SSMSusy_slha_io::fill_minpar_tuple(CNE6SSMSusy_input_parameters<Two_scale>& input,
                                            int key, double value)
{
   switch (key) {
   case 3: input.TanBeta = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CNE6SSMSusy_slha_io::fill_extpar_tuple(CNE6SSMSusy_input_parameters<Two_scale>& input,
                                            int key, double value)
{
   switch (key) {
   case 65: input.sInput = value; break;
   case 66: input.TanTheta = value; break;
   case 67: input.vphiInput = value; break;
   case 72: input.QSInput = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CNE6SSMSusy_slha_io::fill_flexiblesusy_tuple(Spectrum_generator_settings& settings,
                                                  int key, double value)
{
   if (0 <= key && key < static_cast<int>(Spectrum_generator_settings::NUMBER_OF_OPTIONS)) {
      settings.set((Spectrum_generator_settings::Settings)key, value);
   } else {
      WARNING("Unrecognized key in block FlexibleSUSY: " << key);
   }
}

/**
 * Reads the renormalization scales from all DR-bar parameter blocks.
 * If blocks with different scales are found the last scale is
 * returned and a warning is printed.
 *
 * @return common renormalization scale
 */
double CNE6SSMSusy_slha_io::read_scale() const
{
   double scale = 0.;

   for (unsigned i = 0; i < NUMBER_OF_DRBAR_BLOCKS; i++) {
      const double block_scale = slha_io.read_scale(drbar_blocks[i]);
      if (!is_zero(block_scale)) {
         if (!is_zero(scale) && !is_equal(scale, block_scale))
            WARNING("DR-bar parameters defined at different scales");
         scale = block_scale;
      }
   }

   return scale;
}

} // namespace flexiblesusy
