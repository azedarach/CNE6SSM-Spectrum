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

#include "CNE6SSM_slha_io.hpp"
#include "CNE6SSM_semi_two_scale_input_parameters.hpp"
#include "CNE6SSM_two_scale_input_parameters.hpp"
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

#define LOCALPHYSICAL(p) physical.p

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

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CNE6SSM_slha_io::set_extpar(const CNE6SSM_input_parameters<Two_scale>& input)
{
   std::ostringstream extpar;

   extpar << "Block EXTPAR\n";
   extpar << FORMAT_ELEMENT(65, input.ssumInput, "ssumInput");
   extpar << FORMAT_ELEMENT(72, input.QSInput, "QSInput");
   slha_io.set_block(extpar);

}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CNE6SSM_slha_io::set_extpar(const CNE6SSM_semianalytic_input_parameters<Two_scale>& input)
{
   std::ostringstream extpar;

   extpar << "Block EXTPAR\n";
   extpar << FORMAT_ELEMENT(65, input.sInput, "sInput");
   extpar << FORMAT_ELEMENT(72, input.QSInput, "QSInput");
   slha_io.set_block(extpar);

}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CNE6SSM_slha_io::set_minpar(const CNE6SSM_input_parameters<Two_scale>& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(1, input.m0, "m0");
   minpar << FORMAT_ELEMENT(2, input.m12, "m12");
   minpar << FORMAT_ELEMENT(3, input.TanBeta, "TanBeta");
   minpar << FORMAT_ELEMENT(4, input.SignLambdax, "SignLambdax");
   minpar << FORMAT_ELEMENT(5, input.Azero, "Azero");
   slha_io.set_block(minpar);

}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void CNE6SSM_slha_io::set_minpar(const CNE6SSM_semianalytic_input_parameters<Two_scale>& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(2, input.m12, "m12");
   minpar << FORMAT_ELEMENT(3, input.TanBeta, "TanBeta");
   minpar << FORMAT_ELEMENT(5, input.Azero, "Azero");
   slha_io.set_block(minpar);

}

/**
 * Stores the SMINPUTS input parameters in the SLHA object.
 *
 * @param qedqcd class of Standard Model parameters
 */
void CNE6SSM_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void CNE6SSM_slha_io::set_spinfo(const Problems<CNE6SSM_info::NUMBER_OF_PARTICLES>& problems)
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
 * Stores the particle masses in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_masses flag to indicate if Standard Model
 *    particle masses should be written as well
 */
void CNE6SSM_slha_io::set_mass(const CNE6SSM_physical& physical,
                                   bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block MASS\n"
      << FORMAT_MASS(1000021, LOCALPHYSICAL(MGlu), "Glu")
      << FORMAT_MASS(24, LOCALPHYSICAL(MVWm), "VWm")
      << FORMAT_MASS(1000091, LOCALPHYSICAL(MChaP), "ChaP")
      << FORMAT_MASS(31, LOCALPHYSICAL(MVZp), "VZp")
      << FORMAT_MASS(1000092, LOCALPHYSICAL(MChiP(0)), "ChiP(1)")
      << FORMAT_MASS(1000094, LOCALPHYSICAL(MChiP(1)), "ChiP(2)")
      << FORMAT_MASS(1000024, LOCALPHYSICAL(MCha(0)), "Cha(1)")
      << FORMAT_MASS(1000037, LOCALPHYSICAL(MCha(1)), "Cha(2)")
      << FORMAT_MASS(37, LOCALPHYSICAL(MHpm(1)), "Hpm(2)")
      << FORMAT_MASS(92, LOCALPHYSICAL(MSHp0(0)), "SHp0(1)")
      << FORMAT_MASS(94, LOCALPHYSICAL(MSHp0(1)), "SHp0(2)")
      << FORMAT_MASS(91, LOCALPHYSICAL(MSHpp(0)), "SHpp(1)")
      << FORMAT_MASS(93, LOCALPHYSICAL(MSHpp(1)), "SHpp(2)")
      << FORMAT_MASS(1000088, LOCALPHYSICAL(MChaI(0)), "ChaI(1)")
      << FORMAT_MASS(1000089, LOCALPHYSICAL(MChaI(1)), "ChaI(2)")
      << FORMAT_MASS(1000012, LOCALPHYSICAL(MSv(0)), "Sv(1)")
      << FORMAT_MASS(1000014, LOCALPHYSICAL(MSv(1)), "Sv(2)")
      << FORMAT_MASS(1000016, LOCALPHYSICAL(MSv(2)), "Sv(3)")
      << FORMAT_MASS(51, LOCALPHYSICAL(MFDX(0)), "FDX(1)")
      << FORMAT_MASS(52, LOCALPHYSICAL(MFDX(1)), "FDX(2)")
      << FORMAT_MASS(53, LOCALPHYSICAL(MFDX(2)), "FDX(3)")
      << FORMAT_MASS(81, LOCALPHYSICAL(MSHIPM(0)), "SHIPM(1)")
      << FORMAT_MASS(85, LOCALPHYSICAL(MSHIPM(1)), "SHIPM(2)")
      << FORMAT_MASS(83, LOCALPHYSICAL(MSHIPM(2)), "SHIPM(3)")
      << FORMAT_MASS(87, LOCALPHYSICAL(MSHIPM(3)), "SHIPM(4)")
      << FORMAT_MASS(25, LOCALPHYSICAL(Mhh(0)), "hh(1)")
      << FORMAT_MASS(35, LOCALPHYSICAL(Mhh(1)), "hh(2)")
      << FORMAT_MASS(45, LOCALPHYSICAL(Mhh(2)), "hh(3)")
      << FORMAT_MASS(55, LOCALPHYSICAL(Mhh(3)), "hh(4)")
      << FORMAT_MASS(65, LOCALPHYSICAL(Mhh(4)), "hh(5)")
      << FORMAT_MASS(91191138, LOCALPHYSICAL(MAh(2)), "Ah(3)")
      << FORMAT_MASS(36, LOCALPHYSICAL(MAh(3)), "Ah(4)")
      << FORMAT_MASS(91191137, LOCALPHYSICAL(MAh(4)), "Ah(5)")
      << FORMAT_MASS(1000001, LOCALPHYSICAL(MSd(0)), "Sd(1)")
      << FORMAT_MASS(1000003, LOCALPHYSICAL(MSd(1)), "Sd(2)")
      << FORMAT_MASS(1000005, LOCALPHYSICAL(MSd(2)), "Sd(3)")
      << FORMAT_MASS(2000001, LOCALPHYSICAL(MSd(3)), "Sd(4)")
      << FORMAT_MASS(2000003, LOCALPHYSICAL(MSd(4)), "Sd(5)")
      << FORMAT_MASS(2000005, LOCALPHYSICAL(MSd(5)), "Sd(6)")
      << FORMAT_MASS(1000011, LOCALPHYSICAL(MSe(0)), "Se(1)")
      << FORMAT_MASS(1000013, LOCALPHYSICAL(MSe(1)), "Se(2)")
      << FORMAT_MASS(1000015, LOCALPHYSICAL(MSe(2)), "Se(3)")
      << FORMAT_MASS(2000011, LOCALPHYSICAL(MSe(3)), "Se(4)")
      << FORMAT_MASS(2000013, LOCALPHYSICAL(MSe(4)), "Se(5)")
      << FORMAT_MASS(2000015, LOCALPHYSICAL(MSe(5)), "Se(6)")
      << FORMAT_MASS(1000002, LOCALPHYSICAL(MSu(0)), "Su(1)")
      << FORMAT_MASS(1000004, LOCALPHYSICAL(MSu(1)), "Su(2)")
      << FORMAT_MASS(1000006, LOCALPHYSICAL(MSu(2)), "Su(3)")
      << FORMAT_MASS(2000002, LOCALPHYSICAL(MSu(3)), "Su(4)")
      << FORMAT_MASS(2000004, LOCALPHYSICAL(MSu(4)), "Su(5)")
      << FORMAT_MASS(2000006, LOCALPHYSICAL(MSu(5)), "Su(6)")
      << FORMAT_MASS(1000051, LOCALPHYSICAL(MSDX(0)), "SDX(1)")
      << FORMAT_MASS(2000051, LOCALPHYSICAL(MSDX(1)), "SDX(2)")
      << FORMAT_MASS(1000052, LOCALPHYSICAL(MSDX(2)), "SDX(3)")
      << FORMAT_MASS(2000052, LOCALPHYSICAL(MSDX(3)), "SDX(4)")
      << FORMAT_MASS(1000053, LOCALPHYSICAL(MSDX(4)), "SDX(5)")
      << FORMAT_MASS(2000053, LOCALPHYSICAL(MSDX(5)), "SDX(6)")
      << FORMAT_MASS(1000081, LOCALPHYSICAL(MChiI(0)), "ChiI(1)")
      << FORMAT_MASS(1000082, LOCALPHYSICAL(MChiI(1)), "ChiI(2)")
      << FORMAT_MASS(1000083, LOCALPHYSICAL(MChiI(2)), "ChiI(3)")
      << FORMAT_MASS(1000084, LOCALPHYSICAL(MChiI(3)), "ChiI(4)")
      << FORMAT_MASS(1000085, LOCALPHYSICAL(MChiI(4)), "ChiI(5)")
      << FORMAT_MASS(1000086, LOCALPHYSICAL(MChiI(5)), "ChiI(6)")
      << FORMAT_MASS(1000087, LOCALPHYSICAL(MChiI(6)), "ChiI(7)")
      << FORMAT_MASS(82, LOCALPHYSICAL(MSHI0(0)), "SHI0(1)")
      << FORMAT_MASS(86, LOCALPHYSICAL(MSHI0(1)), "SHI0(2)")
      << FORMAT_MASS(84, LOCALPHYSICAL(MSHI0(2)), "SHI0(3)")
      << FORMAT_MASS(88, LOCALPHYSICAL(MSHI0(3)), "SHI0(4)")
      << FORMAT_MASS(9994453, LOCALPHYSICAL(MSHI0(4)), "SHI0(5)")
      << FORMAT_MASS(9994454, LOCALPHYSICAL(MSHI0(5)), "SHI0(6)")
      << FORMAT_MASS(9994455, LOCALPHYSICAL(MSHI0(6)), "SHI0(7)")
      << FORMAT_MASS(1000022, LOCALPHYSICAL(MChi(0)), "Chi(1)")
      << FORMAT_MASS(1000023, LOCALPHYSICAL(MChi(1)), "Chi(2)")
      << FORMAT_MASS(1000025, LOCALPHYSICAL(MChi(2)), "Chi(3)")
      << FORMAT_MASS(1000035, LOCALPHYSICAL(MChi(3)), "Chi(4)")
      << FORMAT_MASS(1000045, LOCALPHYSICAL(MChi(4)), "Chi(5)")
      << FORMAT_MASS(1000055, LOCALPHYSICAL(MChi(5)), "Chi(6)")
      << FORMAT_MASS(1000065, LOCALPHYSICAL(MChi(6)), "Chi(7)")
      << FORMAT_MASS(1000075, LOCALPHYSICAL(MChi(7)), "Chi(8)")
   ;

   if (write_sm_masses) {
      mass
         << FORMAT_MASS(21, LOCALPHYSICAL(MVG), "VG")
         << FORMAT_MASS(12, LOCALPHYSICAL(MFv(0)), "Fv(1)")
         << FORMAT_MASS(14, LOCALPHYSICAL(MFv(1)), "Fv(2)")
         << FORMAT_MASS(16, LOCALPHYSICAL(MFv(2)), "Fv(3)")
         << FORMAT_MASS(22, LOCALPHYSICAL(MVP), "VP")
         << FORMAT_MASS(23, LOCALPHYSICAL(MVZ), "VZ")
         << FORMAT_MASS(11, LOCALPHYSICAL(MFe(0)), "Fe(1)")
         << FORMAT_MASS(13, LOCALPHYSICAL(MFe(1)), "Fe(2)")
         << FORMAT_MASS(15, LOCALPHYSICAL(MFe(2)), "Fe(3)")
         << FORMAT_MASS(1, LOCALPHYSICAL(MFd(0)), "Fd(1)")
         << FORMAT_MASS(3, LOCALPHYSICAL(MFd(1)), "Fd(2)")
         << FORMAT_MASS(5, LOCALPHYSICAL(MFd(2)), "Fd(3)")
         << FORMAT_MASS(2, LOCALPHYSICAL(MFu(0)), "Fu(1)")
         << FORMAT_MASS(4, LOCALPHYSICAL(MFu(1)), "Fu(2)")
         << FORMAT_MASS(6, LOCALPHYSICAL(MFu(2)), "Fu(3)")
      ;
   }

   slha_io.set_block(mass);

}

/**
 * Stores the mixing matrices in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_mixing_matrics flag to indicate if Standard Model
 *    particle mixing matrices should be written as well
 */
void CNE6SSM_slha_io::set_mixing_matrices(const CNE6SSM_physical& physical,
                                              bool write_sm_mixing_matrics)
{
   slha_io.set_block("UHNIMIX", LOCALPHYSICAL(UHI0), "UHI0");
   slha_io.set_block("UHPPMIX", LOCALPHYSICAL(UHIPM), "UHIPM");
   slha_io.set_block("UHNPMIX", LOCALPHYSICAL(UHp0), "UHp0");
   slha_io.set_block("UHPPMIX", LOCALPHYSICAL(UHpp), "UHpp");
   slha_io.set_block("UMIX", LOCALPHYSICAL(UM), "UM");
   slha_io.set_block("VMIX", LOCALPHYSICAL(UP), "UP");
   slha_io.set_block("NMAMIX", LOCALPHYSICAL(ZA), "ZA");
   slha_io.set_block("DSQMIX", LOCALPHYSICAL(ZD), "ZD");
   slha_io.set_block("ESIXZDX", LOCALPHYSICAL(ZDX), "ZDX");
   slha_io.set_block("ESIXZXL", LOCALPHYSICAL(ZDXL), "ZDXL");
   slha_io.set_block("ESIXZXR", LOCALPHYSICAL(ZDXR), "ZDXR");
   slha_io.set_block("SELMIX", LOCALPHYSICAL(ZE), "ZE");
   slha_io.set_block("NMHMIX", LOCALPHYSICAL(ZH), "ZH");
   slha_io.set_block("ESIXZMI", LOCALPHYSICAL(ZMI), "ZMI");
   slha_io.set_block("NMNMIX", LOCALPHYSICAL(ZN), "ZN");
   slha_io.set_block("ZNIMIX", LOCALPHYSICAL(ZNI), "ZNI");
   slha_io.set_block("ZNPMIX", LOCALPHYSICAL(ZNp), "ZNp");
   slha_io.set_block("CHARGEMIX", LOCALPHYSICAL(ZP), "ZP");
   slha_io.set_block("ESIXZPI", LOCALPHYSICAL(ZPI), "ZPI");
   slha_io.set_block("USQMIX", LOCALPHYSICAL(ZU), "ZU");
   slha_io.set_block("SNUMIX", LOCALPHYSICAL(ZV), "ZV");

   if (write_sm_mixing_matrics) {
      slha_io.set_block("UELMIX", LOCALPHYSICAL(ZEL), "ZEL");
      slha_io.set_block("UERMIX", LOCALPHYSICAL(ZER), "ZER");
      slha_io.set_block("UDLMIX", LOCALPHYSICAL(ZDL), "ZDL");
      slha_io.set_block("UDRMIX", LOCALPHYSICAL(ZDR), "ZDR");
      slha_io.set_block("UULMIX", LOCALPHYSICAL(ZUL), "ZUL");
      slha_io.set_block("UURMIX", LOCALPHYSICAL(ZUR), "ZUR");
   }

}

void CNE6SSM_slha_io::set_ckm(
   const Eigen::Matrix<std::complex<double>,3,3>& ckm_matrix,
   double scale)
{
   slha_io.set_block("VCKM"  , ckm_matrix.real(), "Re(CKM)", scale);
   slha_io.set_block("IMVCKM", ckm_matrix.imag(), "Im(CKM)", scale);
}

void CNE6SSM_slha_io::set_pmns(
   const Eigen::Matrix<std::complex<double>,3,3>& pmns_matrix,
   double scale)
{
   slha_io.set_block("VPMNS"  , pmns_matrix.real(), "Re(PMNS)", scale);
   slha_io.set_block("IMVPMNS", pmns_matrix.imag(), "Im(PMNS)", scale);
}

/**
 * Write SLHA object to file.
 *
 * @param file_name file name
 */
void CNE6SSM_slha_io::write_to_file(const std::string& file_name)
{
   slha_io.write_to_file(file_name);
}

/**
 * Read (DR-bar) model parameter output scale from MODSEL entry 12
 */
double CNE6SSM_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

/**
 * Read SLHA object from file
 *
 * @param file_name file name
 */
void CNE6SSM_slha_io::read_from_file(const std::string& file_name)
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
void CNE6SSM_slha_io::fill(CNE6SSM_input_parameters<Two_scale>& input) const
{
   void (*fill_two_scale_minpar_tuple) (CNE6SSM_input_parameters<Two_scale>&,
                                        int, double) = &CNE6SSM_slha_io::fill_minpar_tuple;
   void (*fill_two_scale_extpar_tuple) (CNE6SSM_input_parameters<Two_scale>&,
                                        int, double) = &CNE6SSM_slha_io::fill_extpar_tuple;
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
   input.BMuPhiInput = slha_io.read_entry("ESIXRUNIN", 30);
   slha_io.read_block("ESIXKAPPAIN", input.KappaInput);
   slha_io.read_block("ESIXLAMBDAIN", input.Lambda12Input);
   input.MuPrInput = slha_io.read_entry("ESIXRUNIN", 0);
   input.BMuPrInput = slha_io.read_entry("ESIXRUNIN", 101);

}

/**
 * Fill struct of model input parameters from SLHA object (MINPAR and
 * EXTPAR blocks)
 *
 * @param input struct of model input parameters
 */
void CNE6SSM_slha_io::fill(CNE6SSM_semianalytic_input_parameters<Two_scale>& input) const
{
   void (*fill_semianalytic_minpar_tuple) (CNE6SSM_semianalytic_input_parameters<Two_scale>&,
                                           int, double) = &CNE6SSM_slha_io::fill_minpar_tuple;
   void (*fill_semianalytic_extpar_tuple) (CNE6SSM_semianalytic_input_parameters<Two_scale>&,
                                           int, double) = &CNE6SSM_slha_io::fill_extpar_tuple;
   SLHA_io::Tuple_processor minpar_processor
      = boost::bind(fill_semianalytic_minpar_tuple, boost::ref(input), _1, _2);
   SLHA_io::Tuple_processor extpar_processor
      = boost::bind(fill_semianalytic_extpar_tuple, boost::ref(input), _1, _2);

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);

   input.LambdaxInput = slha_io.read_entry("ESIXRUNIN", 1); 
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
   slha_io.read_block("ESIXLAMBDAIN", input.Lambda12Input);
   input.MuPrInput = slha_io.read_entry("ESIXRUNIN", 0);
   input.BMuPrInput = slha_io.read_entry("ESIXRUNIN", 101);

}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings
 */
void CNE6SSM_slha_io::fill(Spectrum_generator_settings& settings) const
{
   SLHA_io::Tuple_processor flexiblesusy_processor
      = boost::bind(&CNE6SSM_slha_io::fill_flexiblesusy_tuple, boost::ref(settings), _1, _2);

   slha_io.read_block("FlexibleSUSY", flexiblesusy_processor);
}

void CNE6SSM_slha_io::fill_minpar_tuple(CNE6SSM_input_parameters<Two_scale>& input,
                                                int key, double value)
{
   switch (key) {
   case 1: input.m0 = value; break;
   case 2: input.m12 = value; break;
   case 3: input.TanBeta = value; break;
   case 4: input.SignLambdax = value; break;
   case 5: input.Azero = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CNE6SSM_slha_io::fill_minpar_tuple(CNE6SSM_semianalytic_input_parameters<Two_scale>& input,
                                                int key, double value)
{
   switch (key) {
   case 2: input.m12 = value; break;
   case 3: input.TanBeta = value; break;
   case 5: input.Azero = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CNE6SSM_slha_io::fill_extpar_tuple(CNE6SSM_input_parameters<Two_scale>& input,
                                                int key, double value)
{
   switch (key) {
   case 65: input.ssumInput = value; break;
   case 72: input.QSInput = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void CNE6SSM_slha_io::fill_extpar_tuple(CNE6SSM_semianalytic_input_parameters<Two_scale>& input,
                                                int key, double value)
{
   switch (key) {
   case 65: input.sInput = value; break;
   case 72: input.QSInput = value; break;
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
