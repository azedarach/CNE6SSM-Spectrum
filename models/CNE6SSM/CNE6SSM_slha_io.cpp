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

// File generated at Mon 15 Sep 2014 17:29:09

#include "CNE6SSM_slha_io.hpp"
#include "CNE6SSM_input_parameters.hpp"
#include "logger.hpp"
#include "wrappers.hpp"
#include "numerics.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "config.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/bind.hpp>

using namespace softsusy;

namespace flexiblesusy {

char const * const CNE6SSM_slha_io::drbar_blocks[NUMBER_OF_DRBAR_BLOCKS] =
   { "gauge", "Yu", "Yd", "Ye", "Te", "Td", "Tu", "HMIX", "ESIXHEYUK",
   "ESIXRUN", "ESIXGDYUK", "ESIXFUYUK", "ESIXFDYUK", "ESIXTHETRI", "ESIXTGDTRI"
   , "ESIXTFUTRI", "ESIXTFDTRI", "MSQ2", "MSE2", "MSL2", "MSU2", "MSD2",
   "MSOFT", "mX2", "mXBar2", "ESIXKAPPA", "ESIXTKAPPA", "ESIXLAMBDA",
   "ESIXTLAMBDA" }
;

CNE6SSM_slha_io::CNE6SSM_slha_io()
   : slha_io()
{
}

void CNE6SSM_slha_io::clear()
{
   slha_io.clear();
}

void CNE6SSM_slha_io::set_extpar(const CNE6SSM_input_parameters& input)
{
   std::ostringstream extpar;

   extpar << "Block EXTPAR\n";
   extpar << FORMAT_ELEMENT(65, input.ssumInput, "ssumInput");
   extpar << FORMAT_ELEMENT(72, input.QS, "QS");
   slha_io.set_block(extpar);

}

void CNE6SSM_slha_io::set_minpar(const CNE6SSM_input_parameters& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(1, input.m0, "m0");
   minpar << FORMAT_ELEMENT(2, input.m12, "m12");
   minpar << FORMAT_ELEMENT(3, input.TanBeta, "TanBeta");
   minpar << FORMAT_ELEMENT(5, input.Azero, "Azero");
   slha_io.set_block(minpar);

}

void CNE6SSM_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

void CNE6SSM_slha_io::set_spinfo(const Problems<CNE6SSM_info::NUMBER_OF_PARTICLES>& problems)
{
   std::ostringstream spinfo;
   spinfo << "# FlexibleSUSY " FLEXIBLESUSY_VERSION " SLHA compliant output\n"
             "# P. Athron, Jae-hyeon Park, D. Stöckinger, A. Voigt\n"
             "Block SPINFO\n"
          << FORMAT_SPINFO(1, PKGNAME)
          << FORMAT_SPINFO(2, FLEXIBLESUSY_VERSION);

   if (problems.have_serious_problem()) {
      std::ostringstream serious_problems;
      problems.print(serious_problems);
      spinfo << FORMAT_SPINFO(4, serious_problems.str());
   }

   slha_io.set_block(spinfo, SLHA_io::front);
}

void CNE6SSM_slha_io::write_to_file(const std::string& file_name)
{
   slha_io.write_to_file(file_name);
}

double CNE6SSM_slha_io::get_input_scale() const
{
   return slha_io.get_extpar().input_scale;
}

double CNE6SSM_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

void CNE6SSM_slha_io::read_from_file(const std::string& file_name)
{
   slha_io.read_from_file(file_name);
   slha_io.read_modsel();
   slha_io.read_extpar();
}

void CNE6SSM_slha_io::fill(CNE6SSM_input_parameters& input) const
{
   SLHA_io::Tuple_processor minpar_processor
      = boost::bind(&CNE6SSM_slha_io::fill_minpar_tuple, boost::ref(input), _1, _2);
   SLHA_io::Tuple_processor extpar_processor
      = boost::bind(&CNE6SSM_slha_io::fill_extpar_tuple, boost::ref(input), _1, _2);

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
   input.BMuPhiInput = slha_io.read_entry("ESIXRUNIN", 30);
   slha_io.read_block("ESIXKAPPAIN", input.KappaInput);
   input.LambdaxInput = slha_io.read_entry("ESIXRUNIN", 1);
   slha_io.read_block("ESIXLAMBDAIN", input.Lambda12Input);
   input.MuPrInput = slha_io.read_entry("ESIXRUNIN", 0);
   input.BMuPrInput = slha_io.read_entry("ESIXRUNIN", 101);

}

void CNE6SSM_slha_io::fill(Spectrum_generator_settings& settings) const
{
   SLHA_io::Tuple_processor flexiblesusy_processor
      = boost::bind(&CNE6SSM_slha_io::fill_flexiblesusy_tuple, boost::ref(settings), _1, _2);

   slha_io.read_block("FlexibleSUSY", flexiblesusy_processor);
}

void CNE6SSM_slha_io::fill_minpar_tuple(CNE6SSM_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   case 1: input.m0 = value; break;
   case 2: input.m12 = value; break;
   case 3: input.TanBeta = value; break;
   case 5: input.Azero = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CNE6SSM_slha_io::fill_extpar_tuple(CNE6SSM_input_parameters& input,
                                                int key, double value)
{
   // key 0 is the model parameter input scale, which is read in
   // slha_io.{hpp,cpp}
   if (key == 0)
      return;

   switch (key) {
   case 65: input.ssumInput = value; break;
   case 72: input.QS = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CNE6SSM_slha_io::fill_flexiblesusy_tuple(Spectrum_generator_settings& settings,
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
double CNE6SSM_slha_io::read_scale() const
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
