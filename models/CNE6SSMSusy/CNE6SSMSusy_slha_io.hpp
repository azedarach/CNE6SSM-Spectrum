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

#ifndef CNE6SSMSusy_SLHA_IO_H
#define CNE6SSMSusy_SLHA_IO_H

#include "CNE6SSMSusy_two_scale_model_slha.hpp"
#include "CNE6SSMSusy_info.hpp"
#include "slha_io.hpp"
#include "ew_input.hpp"
#include "lowe.h"

#include <Eigen/Core>
#include <string>
#include <utility>

#define Pole(p) physical.p
#define MODELPARAMETER(p) model.get_##p()
#define DEFINE_PARAMETER(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(MODELPARAMETER(p))>::type>::type p;
#define LowEnergyConstant(p) Electroweak_constants::p
#define SCALES(p) scales.p

namespace flexiblesusy {

template <class T>
struct CNE6SSMSusy_input_parameters;

class Two_scale;
class Spectrum_generator_settings;

struct CNE6SSMSusy_scales {
   CNE6SSMSusy_scales() : HighScale(0.), LowScale(0.) {}
   double HighScale, LowScale;
};

class CNE6SSMSusy_slha_io {
public:
   CNE6SSMSusy_slha_io();
   ~CNE6SSMSusy_slha_io() {}

   void clear();

   void fill(QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(CNE6SSMSusy_input_parameters<Two_scale>&) const;
   template <class T> void fill(CNE6SSMSusy_slha<T>&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void set_extpar(const CNE6SSMSusy_input_parameters<Two_scale>&);
   template <class T> void set_extra(const CNE6SSMSusy_slha<T>&, const CNE6SSMSusy_scales&);
   void set_minpar(const CNE6SSMSusy_input_parameters<Two_scale>&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class T> void set_spectrum(const CNE6SSMSusy_slha<T>&);
   template <class T> void set_spectrum(const CNE6SSMSusy<T>&);
   void set_spinfo(const Problems<CNE6SSMSusy_info::NUMBER_OF_PARTICLES>&);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& ostr = std::cout) { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(CNE6SSMSusy_input_parameters<Two_scale>&, int, double);
   static void fill_extpar_tuple(CNE6SSMSusy_input_parameters<Two_scale>&, int, double);
   static void fill_flexiblesusy_tuple(Spectrum_generator_settings&, int, double);

   template <class T>
   static void fill_slhaea(SLHAea::Coll&, const CNE6SSMSusy_slha<T>&, const QedQcd&, const CNE6SSMSusy_scales&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const CNE6SSMSusy_slha<T>&, const QedQcd&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const CNE6SSMSusy_slha<T>&, const QedQcd&, const CNE6SSMSusy_scales&);

private:
   SLHA_io slha_io; ///< SLHA io class
   static unsigned const NUMBER_OF_DRBAR_BLOCKS = 12;
   static char const * const drbar_blocks[NUMBER_OF_DRBAR_BLOCKS];

   template <class T> void set_model_parameters(const CNE6SSMSusy_slha<T>&);
   double read_scale() const;
   template <class T> void fill_drbar_parameters(CNE6SSMSusy_slha<T>&) const;
};

/**
 * Reads SUSY parameters from a SLHA output file.
 */
template <class T>
void CNE6SSMSusy_slha_io::fill(CNE6SSMSusy_slha<T>& model) const
{
   fill_drbar_parameters(model);
}

/**
 * Reads DR-bar parameters from a SLHA output file.
 */
template <class T>
void CNE6SSMSusy_slha_io::fill_drbar_parameters(CNE6SSMSusy_slha<T>& model) const
{
   model.set_g1(slha_io.read_entry("gauge", 1) * 1.2909944487358056);
   model.set_g2(slha_io.read_entry("gauge", 2));
   model.set_g3(slha_io.read_entry("gauge", 3));
   model.set_g1p(slha_io.read_entry("gauge", 4));
   {
      DEFINE_PARAMETER(Yu);
      slha_io.read_block("Yu", Yu);
      model.set_Yu(Yu);
   }
   {
      DEFINE_PARAMETER(Yd);
      slha_io.read_block("Yd", Yd);
      model.set_Yd(Yd);
   }
   {
      DEFINE_PARAMETER(Ye);
      slha_io.read_block("Ye", Ye);
      model.set_Ye(Ye);
   }
   model.set_XiF(slha_io.read_entry("HMIX", 30));
   model.set_MuPhi(slha_io.read_entry("HMIX", 31));
   model.set_KappaPr(slha_io.read_entry("HMIX", 32));
   model.set_Sigmax(slha_io.read_entry("HMIX", 33));
   {
      DEFINE_PARAMETER(hE);
      slha_io.read_block("ESIXHEYUK", hE);
      model.set_hE(hE);
   }
   model.set_SigmaL(slha_io.read_entry("ESIXRUN", 42));
   {
      DEFINE_PARAMETER(gD);
      slha_io.read_block("ESIXGDYUK", gD);
      model.set_gD(gD);
   }
   {
      DEFINE_PARAMETER(fu);
      slha_io.read_block("ESIXFUYUK", fu);
      model.set_fu(fu);
   }
   {
      DEFINE_PARAMETER(fd);
      slha_io.read_block("ESIXFDYUK", fd);
      model.set_fd(fd);
   }
   model.set_vd(slha_io.read_entry("HMIX", 102));
   model.set_vu(slha_io.read_entry("HMIX", 103));
   model.set_vs(slha_io.read_entry("ESIXRUN", 11));
   model.set_vsb(slha_io.read_entry("ESIXRUN", 12));
   model.set_vphi(slha_io.read_entry("ESIXRUN", 13));
   {
      DEFINE_PARAMETER(Kappa);
      slha_io.read_block("ESIXKAPPA", Kappa);
      model.set_Kappa(Kappa);
   }
   model.set_Lambdax(slha_io.read_entry("ESIXRUN", 1));
   {
      DEFINE_PARAMETER(Lambda12);
      slha_io.read_block("ESIXLAMBDA", Lambda12);
      model.set_Lambda12(Lambda12);
   }
   model.set_MuPr(slha_io.read_entry("ESIXRUN", 0));


   model.set_scale(read_scale());
}

template <class T>
void CNE6SSMSusy_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const CNE6SSMSusy_slha<T>& model,
   const QedQcd& qedqcd, const CNE6SSMSusy_scales& scales)
{
   CNE6SSMSusy_slha_io slha_io;
   const CNE6SSMSusy_input_parameters<T>& input = model.get_input();
   const Problems<CNE6SSMSusy_info::NUMBER_OF_PARTICLES>& problems
      = model.get_problems();
   const bool error = problems.have_problem();

   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(qedqcd);
   slha_io.set_minpar(input);
   slha_io.set_extpar(input);
   if (!error) {
      slha_io.set_spectrum(model);
      slha_io.set_extra(model, scales);
   }

   slhaea = slha_io.get_slha_io().get_data();
}

template <class T>
SLHAea::Coll CNE6SSMSusy_slha_io::fill_slhaea(
   const CNE6SSMSusy_slha<T>& model, const QedQcd& qedqcd)
{
   CNE6SSMSusy_scales scales;

   return fill_slhaea(model, qedqcd, scales);
}

template <class T>
SLHAea::Coll CNE6SSMSusy_slha_io::fill_slhaea(
   const CNE6SSMSusy_slha<T>& model, const QedQcd& qedqcd,
   const CNE6SSMSusy_scales& scales)
{
   SLHAea::Coll slhaea;
   CNE6SSMSusy_slha_io::fill_slhaea(slhaea, model, qedqcd, scales);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class T>
void CNE6SSMSusy_slha_io::set_model_parameters(const CNE6SSMSusy_slha<T>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "gY")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
            << FORMAT_ELEMENT(4, (MODELPARAMETER(g1p)), "g1p")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block("Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block("Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(30, (MODELPARAMETER(XiF)), "XiF")
            << FORMAT_ELEMENT(31, (MODELPARAMETER(MuPhi)), "MuPhi")
            << FORMAT_ELEMENT(32, (MODELPARAMETER(KappaPr)), "KappaPr")
            << FORMAT_ELEMENT(33, (MODELPARAMETER(Sigmax)), "Sigmax")
            << FORMAT_ELEMENT(102, (MODELPARAMETER(vd)), "vd")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(vu)), "vu")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("ESIXHEYUK", MODELPARAMETER(hE), "hE", model.get_scale());
   {
      std::ostringstream block;
      block << "Block ESIXRUN Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(42, (MODELPARAMETER(SigmaL)), "SigmaL")
            << FORMAT_ELEMENT(11, (MODELPARAMETER(vs)), "vs")
            << FORMAT_ELEMENT(12, (MODELPARAMETER(vsb)), "vsb")
            << FORMAT_ELEMENT(13, (MODELPARAMETER(vphi)), "vphi")
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Lambdax)), "Lambdax")
            << FORMAT_ELEMENT(0, (MODELPARAMETER(MuPr)), "MuPr")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("ESIXGDYUK", MODELPARAMETER(gD), "gD", model.get_scale());
   slha_io.set_block("ESIXFUYUK", MODELPARAMETER(fu), "fu", model.get_scale());
   slha_io.set_block("ESIXFDYUK", MODELPARAMETER(fd), "fd", model.get_scale());
   slha_io.set_block("ESIXKAPPA", MODELPARAMETER(Kappa), "Kappa", model.get_scale());
   slha_io.set_block("ESIXLAMBDA", MODELPARAMETER(Lambda12), "Lambda12", model.get_scale());

}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 */
template <class T>
void CNE6SSMSusy_slha_io::set_extra(
   const CNE6SSMSusy_slha<T>& model, const CNE6SSMSusy_scales& scales)
{

   {
      std::ostringstream block;
      block << "Block FlexibleSUSYOutput Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(0, (SCALES(HighScale)), "HighScale")
            << FORMAT_ELEMENT(2, (SCALES(LowScale)), "LowScale")
      ;
      slha_io.set_block(block);
   }

}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void CNE6SSMSusy_slha_io::set_spectrum(const CNE6SSMSusy<T>& model)
{
   const CNE6SSMSusy_slha<T> model_slha(model);
   set_spectrum(model_slha);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class T>
void CNE6SSMSusy_slha_io::set_spectrum(const CNE6SSMSusy_slha<T>& model)
{
   set_model_parameters(model);
}

} // namespace flexiblesusy

#undef Pole
#undef MODELPARAMETER
#undef DEFINE_PARAMETER
#undef LowEnergyConstant
#undef SCALES

#endif
