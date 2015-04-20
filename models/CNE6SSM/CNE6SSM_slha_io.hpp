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

#ifndef CNE6SSM_SLHA_IO_H
#define CNE6SSM_SLHA_IO_H

#include "CNE6SSM_two_scale_model_slha.hpp"
#include "CNE6SSM_info.hpp"
#include "CNE6SSM_physical.hpp"
#include "slha_io.hpp"
#include "ckm.hpp"
#include "ew_input.hpp"

#include <Eigen/Core>
#include <string>
#include <utility>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) model.get_##p()
#define DEFINE_PARAMETER(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(MODELPARAMETER(p))>::type>::type p;
#define DEFINE_POLE_MASS(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(PHYSICAL(p))>::type>::type p;
#define SM(p) Electroweak_constants::p
#define SCALES(p) scales.p

namespace flexiblesusy {

struct CNE6SSM_input_parameters;
class Spectrum_generator_settings;

struct CNE6SSM_scales {
   CNE6SSM_scales() : HighScale(0.), SUSYScale(0.), LowScale(0.) {}
   double HighScale, SUSYScale, LowScale;
};

class CNE6SSM_slha_io {
public:
   CNE6SSM_slha_io();
   ~CNE6SSM_slha_io() {}

   void clear();

   void fill(QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(CNE6SSM_input_parameters&) const;
   template <class T> void fill(CNE6SSM_slha<T>&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void set_extpar(const CNE6SSM_input_parameters&);
   template <class T> void set_extra(const CNE6SSM_slha<T>&, const CNE6SSM_scales&);
   void set_minpar(const CNE6SSM_input_parameters&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class T> void set_spectrum(const CNE6SSM_slha<T>&);
   template <class T> void set_spectrum(const CNE6SSM<T>&);
   void set_spinfo(const Problems<CNE6SSM_info::NUMBER_OF_PARTICLES>&);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& ostr = std::cout) { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(CNE6SSM_input_parameters&, int, double);
   static void fill_extpar_tuple(CNE6SSM_input_parameters&, int, double);
   static void fill_flexiblesusy_tuple(Spectrum_generator_settings&, int, double);

   template <class T>
   static void fill_slhaea(SLHAea::Coll&, const CNE6SSM_slha<T>&, const QedQcd&, const CNE6SSM_scales&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const CNE6SSM_slha<T>&, const QedQcd&);

   template <class T>
   static SLHAea::Coll fill_slhaea(const CNE6SSM_slha<T>&, const QedQcd&, const CNE6SSM_scales&);

private:
   SLHA_io slha_io; ///< SLHA io class
   static unsigned const NUMBER_OF_DRBAR_BLOCKS = 29;
   static char const * const drbar_blocks[NUMBER_OF_DRBAR_BLOCKS];

   void set_mass(const CNE6SSM_physical&, bool);
   void set_mixing_matrices(const CNE6SSM_physical&, bool);
   template <class T> void set_model_parameters(const CNE6SSM_slha<T>&);
   void set_ckm(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   void set_pmns(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   double read_scale() const;
   template <class T> void fill_drbar_parameters(CNE6SSM_slha<T>&) const;
   template <class T> void fill_physical(CNE6SSM_slha<T>&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class T>
void CNE6SSM_slha_io::fill(CNE6SSM_slha<T>& model) const
{
   fill_drbar_parameters(model);
   fill_physical(model);
}

/**
 * Reads DR-bar parameters from a SLHA output file.
 */
template <class T>
void CNE6SSM_slha_io::fill_drbar_parameters(CNE6SSM_slha<T>& model) const
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
   {
      DEFINE_PARAMETER(TYe);
      slha_io.read_block("Te", TYe);
      model.set_TYe(TYe);
   }
   {
      DEFINE_PARAMETER(TYd);
      slha_io.read_block("Td", TYd);
      model.set_TYd(TYd);
   }
   {
      DEFINE_PARAMETER(TYu);
      slha_io.read_block("Tu", TYu);
      model.set_TYu(TYu);
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
   model.set_TKappaPr(slha_io.read_entry("ESIXRUN", 28));
   model.set_TSigmax(slha_io.read_entry("ESIXRUN", 29));
   {
      DEFINE_PARAMETER(ThE);
      slha_io.read_block("ESIXTHETRI", ThE);
      model.set_ThE(ThE);
   }
   model.set_TSigmaL(slha_io.read_entry("ESIXRUN", 43));
   {
      DEFINE_PARAMETER(TgD);
      slha_io.read_block("ESIXTGDTRI", TgD);
      model.set_TgD(TgD);
   }
   {
      DEFINE_PARAMETER(Tfu);
      slha_io.read_block("ESIXTFUTRI", Tfu);
      model.set_Tfu(Tfu);
   }
   {
      DEFINE_PARAMETER(Tfd);
      slha_io.read_block("ESIXTFDTRI", Tfd);
      model.set_Tfd(Tfd);
   }
   model.set_BMuPhi(slha_io.read_entry("ESIXRUN", 30));
   model.set_LXiF(slha_io.read_entry("HMIX", 34));
   {
      DEFINE_PARAMETER(mq2);
      slha_io.read_block("MSQ2", mq2);
      model.set_mq2(mq2);
   }
   {
      DEFINE_PARAMETER(me2);
      slha_io.read_block("MSE2", me2);
      model.set_me2(me2);
   }
   {
      DEFINE_PARAMETER(ml2);
      slha_io.read_block("MSL2", ml2);
      model.set_ml2(ml2);
   }
   {
      DEFINE_PARAMETER(mu2);
      slha_io.read_block("MSU2", mu2);
      model.set_mu2(mu2);
   }
   {
      DEFINE_PARAMETER(md2);
      slha_io.read_block("MSD2", md2);
      model.set_md2(md2);
   }
   model.set_mHd2(slha_io.read_entry("MSOFT", 21));
   model.set_mHu2(slha_io.read_entry("MSOFT", 22));
   {
      DEFINE_PARAMETER(mDx2);
      slha_io.read_block("mX2", mDx2);
      model.set_mDx2(mDx2);
   }
   {
      DEFINE_PARAMETER(mDxbar2);
      slha_io.read_block("mXBar2", mDxbar2);
      model.set_mDxbar2(mDxbar2);
   }
   model.set_ms2(slha_io.read_entry("MSOFT", 23));
   model.set_msbar2(slha_io.read_entry("MSOFT", 24));
   model.set_mphi2(slha_io.read_entry("MSOFT", 25));
   model.set_mHp2(slha_io.read_entry("MSOFT", 26));
   model.set_mHpbar2(slha_io.read_entry("MSOFT", 27));
   model.set_MassB(slha_io.read_entry("MSOFT", 1));
   model.set_MassWB(slha_io.read_entry("MSOFT", 2));
   model.set_MassG(slha_io.read_entry("MSOFT", 3));
   model.set_MassBp(slha_io.read_entry("MSOFT", 4));
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
   {
      DEFINE_PARAMETER(TKappa);
      slha_io.read_block("ESIXTKAPPA", TKappa);
      model.set_TKappa(TKappa);
   }
   model.set_Lambdax(slha_io.read_entry("ESIXRUN", 1));
   model.set_TLambdax(slha_io.read_entry("ESIXRUN", 2));
   {
      DEFINE_PARAMETER(Lambda12);
      slha_io.read_block("ESIXLAMBDA", Lambda12);
      model.set_Lambda12(Lambda12);
   }
   {
      DEFINE_PARAMETER(TLambda12);
      slha_io.read_block("ESIXTLAMBDA", TLambda12);
      model.set_TLambda12(TLambda12);
   }
   model.set_MuPr(slha_io.read_entry("ESIXRUN", 0));
   model.set_BMuPr(slha_io.read_entry("ESIXRUN", 101));


   model.set_scale(read_scale());
}

/**
 * Reads pole masses and mixing matrices from a SLHA output file.
 */
template <class T>
void CNE6SSM_slha_io::fill_physical(CNE6SSM_slha<T>& model) const
{
   {
      DEFINE_PARAMETER(ZD);
      slha_io.read_block("DSQMIX", ZD);
      PHYSICAL(ZD) = ZD;
   }
   {
      DEFINE_PARAMETER(ZU);
      slha_io.read_block("USQMIX", ZU);
      PHYSICAL(ZU) = ZU;
   }
   {
      DEFINE_PARAMETER(ZE);
      slha_io.read_block("SELMIX", ZE);
      PHYSICAL(ZE) = ZE;
   }
   {
      DEFINE_PARAMETER(ZV);
      slha_io.read_block("SNUMIX", ZV);
      PHYSICAL(ZV) = ZV;
   }
   {
      DEFINE_PARAMETER(ZDX);
      slha_io.read_block("ESIXZDX", ZDX);
      PHYSICAL(ZDX) = ZDX;
   }
   {
      DEFINE_PARAMETER(ZH);
      slha_io.read_block("NMHMIX", ZH);
      PHYSICAL(ZH) = ZH;
   }
   {
      DEFINE_PARAMETER(ZA);
      slha_io.read_block("NMAMIX", ZA);
      PHYSICAL(ZA) = ZA;
   }
   {
      DEFINE_PARAMETER(ZP);
      slha_io.read_block("CHARGEMIX", ZP);
      PHYSICAL(ZP) = ZP;
   }
   {
      DEFINE_PARAMETER(ZN);
      slha_io.read_block("NMNMIX", ZN);
      PHYSICAL(ZN) = ZN;
   }
   {
      DEFINE_PARAMETER(ZNp);
      slha_io.read_block("ZNPMIX", ZNp);
      PHYSICAL(ZNp) = ZNp;
   }
   {
      DEFINE_PARAMETER(ZNI);
      slha_io.read_block("ZNIMIX", ZNI);
      PHYSICAL(ZNI) = ZNI;
   }
   {
      DEFINE_PARAMETER(ZMI);
      slha_io.read_block("ESIXZMI", ZMI);
      PHYSICAL(ZMI) = ZMI;
   }
   {
      DEFINE_PARAMETER(ZPI);
      slha_io.read_block("ESIXZPI", ZPI);
      PHYSICAL(ZPI) = ZPI;
   }
   {
      DEFINE_PARAMETER(UP);
      slha_io.read_block("VMIX", UP);
      PHYSICAL(UP) = UP;
   }
   {
      DEFINE_PARAMETER(UM);
      slha_io.read_block("UMIX", UM);
      PHYSICAL(UM) = UM;
   }
   {
      DEFINE_PARAMETER(ZEL);
      slha_io.read_block("UELMIX", ZEL);
      PHYSICAL(ZEL) = ZEL;
   }
   {
      DEFINE_PARAMETER(ZER);
      slha_io.read_block("UERMIX", ZER);
      PHYSICAL(ZER) = ZER;
   }
   {
      DEFINE_PARAMETER(ZDL);
      slha_io.read_block("UDLMIX", ZDL);
      PHYSICAL(ZDL) = ZDL;
   }
   {
      DEFINE_PARAMETER(ZDR);
      slha_io.read_block("UDRMIX", ZDR);
      PHYSICAL(ZDR) = ZDR;
   }
   {
      DEFINE_PARAMETER(ZUL);
      slha_io.read_block("UULMIX", ZUL);
      PHYSICAL(ZUL) = ZUL;
   }
   {
      DEFINE_PARAMETER(ZUR);
      slha_io.read_block("UURMIX", ZUR);
      PHYSICAL(ZUR) = ZUR;
   }
   {
      DEFINE_PARAMETER(ZDXL);
      slha_io.read_block("ESIXZXL", ZDXL);
      PHYSICAL(ZDXL) = ZDXL;
   }
   {
      DEFINE_PARAMETER(ZDXR);
      slha_io.read_block("ESIXZXR", ZDXR);
      PHYSICAL(ZDXR) = ZDXR;
   }
   {
      DEFINE_PARAMETER(UHp0);
      slha_io.read_block("UHNPMIX", UHp0);
      PHYSICAL(UHp0) = UHp0;
   }
   {
      DEFINE_PARAMETER(UHpp);
      slha_io.read_block("UHPPMIX", UHpp);
      PHYSICAL(UHpp) = UHpp;
   }
   {
      DEFINE_PARAMETER(UHI0);
      slha_io.read_block("UHNIMIX", UHI0);
      PHYSICAL(UHI0) = UHI0;
   }
   {
      DEFINE_PARAMETER(UHIPM);
      slha_io.read_block("UHPPMIX", UHIPM);
      PHYSICAL(UHIPM) = UHIPM;
   }

   PHYSICAL(MVG) = slha_io.read_entry("MASS", 21);
   PHYSICAL(MGlu) = slha_io.read_entry("MASS", 1000021);
   PHYSICAL(MFv)(0) = slha_io.read_entry("MASS", 12);
   PHYSICAL(MFv)(1) = slha_io.read_entry("MASS", 14);
   PHYSICAL(MFv)(2) = slha_io.read_entry("MASS", 16);
   PHYSICAL(MChaP) = slha_io.read_entry("MASS", 1000091);
   PHYSICAL(MVP) = slha_io.read_entry("MASS", 22);
   PHYSICAL(MVZ) = slha_io.read_entry("MASS", 23);
   PHYSICAL(MVZp) = slha_io.read_entry("MASS", 31);
   PHYSICAL(MSd)(0) = slha_io.read_entry("MASS", 1000001);
   PHYSICAL(MSd)(1) = slha_io.read_entry("MASS", 1000003);
   PHYSICAL(MSd)(2) = slha_io.read_entry("MASS", 1000005);
   PHYSICAL(MSd)(3) = slha_io.read_entry("MASS", 2000001);
   PHYSICAL(MSd)(4) = slha_io.read_entry("MASS", 2000003);
   PHYSICAL(MSd)(5) = slha_io.read_entry("MASS", 2000005);
   PHYSICAL(MSv)(0) = slha_io.read_entry("MASS", 1000012);
   PHYSICAL(MSv)(1) = slha_io.read_entry("MASS", 1000014);
   PHYSICAL(MSv)(2) = slha_io.read_entry("MASS", 1000016);
   PHYSICAL(MSu)(0) = slha_io.read_entry("MASS", 1000002);
   PHYSICAL(MSu)(1) = slha_io.read_entry("MASS", 1000004);
   PHYSICAL(MSu)(2) = slha_io.read_entry("MASS", 1000006);
   PHYSICAL(MSu)(3) = slha_io.read_entry("MASS", 2000002);
   PHYSICAL(MSu)(4) = slha_io.read_entry("MASS", 2000004);
   PHYSICAL(MSu)(5) = slha_io.read_entry("MASS", 2000006);
   PHYSICAL(MSe)(0) = slha_io.read_entry("MASS", 1000011);
   PHYSICAL(MSe)(1) = slha_io.read_entry("MASS", 1000013);
   PHYSICAL(MSe)(2) = slha_io.read_entry("MASS", 1000015);
   PHYSICAL(MSe)(3) = slha_io.read_entry("MASS", 2000011);
   PHYSICAL(MSe)(4) = slha_io.read_entry("MASS", 2000013);
   PHYSICAL(MSe)(5) = slha_io.read_entry("MASS", 2000015);
   PHYSICAL(MSDX)(0) = slha_io.read_entry("MASS", 1000051);
   PHYSICAL(MSDX)(1) = slha_io.read_entry("MASS", 2000051);
   PHYSICAL(MSDX)(2) = slha_io.read_entry("MASS", 1000052);
   PHYSICAL(MSDX)(3) = slha_io.read_entry("MASS", 2000052);
   PHYSICAL(MSDX)(4) = slha_io.read_entry("MASS", 1000053);
   PHYSICAL(MSDX)(5) = slha_io.read_entry("MASS", 2000053);
   PHYSICAL(Mhh)(0) = slha_io.read_entry("MASS", 25);
   PHYSICAL(Mhh)(1) = slha_io.read_entry("MASS", 35);
   PHYSICAL(Mhh)(2) = slha_io.read_entry("MASS", 45);
   PHYSICAL(Mhh)(3) = slha_io.read_entry("MASS", 55);
   PHYSICAL(Mhh)(4) = slha_io.read_entry("MASS", 65);
   PHYSICAL(MAh)(2) = slha_io.read_entry("MASS", 91191138);
   PHYSICAL(MAh)(3) = slha_io.read_entry("MASS", 36);
   PHYSICAL(MAh)(4) = slha_io.read_entry("MASS", 91191137);
   PHYSICAL(MHpm)(1) = slha_io.read_entry("MASS", 37);
   PHYSICAL(MChi)(0) = slha_io.read_entry("MASS", 1000022);
   PHYSICAL(MChi)(1) = slha_io.read_entry("MASS", 1000023);
   PHYSICAL(MChi)(2) = slha_io.read_entry("MASS", 1000025);
   PHYSICAL(MChi)(3) = slha_io.read_entry("MASS", 1000035);
   PHYSICAL(MChi)(4) = slha_io.read_entry("MASS", 1000045);
   PHYSICAL(MChi)(5) = slha_io.read_entry("MASS", 1000055);
   PHYSICAL(MChi)(6) = slha_io.read_entry("MASS", 1000065);
   PHYSICAL(MChi)(7) = slha_io.read_entry("MASS", 1000075);
   PHYSICAL(MCha)(0) = slha_io.read_entry("MASS", 1000024);
   PHYSICAL(MCha)(1) = slha_io.read_entry("MASS", 1000037);
   PHYSICAL(MFe)(0) = slha_io.read_entry("MASS", 11);
   PHYSICAL(MFe)(1) = slha_io.read_entry("MASS", 13);
   PHYSICAL(MFe)(2) = slha_io.read_entry("MASS", 15);
   PHYSICAL(MFd)(0) = slha_io.read_entry("MASS", 1);
   PHYSICAL(MFd)(1) = slha_io.read_entry("MASS", 3);
   PHYSICAL(MFd)(2) = slha_io.read_entry("MASS", 5);
   PHYSICAL(MFu)(0) = slha_io.read_entry("MASS", 2);
   PHYSICAL(MFu)(1) = slha_io.read_entry("MASS", 4);
   PHYSICAL(MFu)(2) = slha_io.read_entry("MASS", 6);
   PHYSICAL(MFDX)(0) = slha_io.read_entry("MASS", 51);
   PHYSICAL(MFDX)(1) = slha_io.read_entry("MASS", 52);
   PHYSICAL(MFDX)(2) = slha_io.read_entry("MASS", 53);
   PHYSICAL(MSHI0)(0) = slha_io.read_entry("MASS", 82);
   PHYSICAL(MSHI0)(1) = slha_io.read_entry("MASS", 86);
   PHYSICAL(MSHI0)(2) = slha_io.read_entry("MASS", 84);
   PHYSICAL(MSHI0)(3) = slha_io.read_entry("MASS", 88);
   PHYSICAL(MSHI0)(4) = slha_io.read_entry("MASS", 9994453);
   PHYSICAL(MSHI0)(5) = slha_io.read_entry("MASS", 9994454);
   PHYSICAL(MSHI0)(6) = slha_io.read_entry("MASS", 9994455);
   PHYSICAL(MSHIPM)(0) = slha_io.read_entry("MASS", 81);
   PHYSICAL(MSHIPM)(1) = slha_io.read_entry("MASS", 85);
   PHYSICAL(MSHIPM)(2) = slha_io.read_entry("MASS", 83);
   PHYSICAL(MSHIPM)(3) = slha_io.read_entry("MASS", 87);
   PHYSICAL(MChaI)(0) = slha_io.read_entry("MASS", 1000088);
   PHYSICAL(MChaI)(1) = slha_io.read_entry("MASS", 1000089);
   PHYSICAL(MChiI)(0) = slha_io.read_entry("MASS", 1000081);
   PHYSICAL(MChiI)(1) = slha_io.read_entry("MASS", 1000082);
   PHYSICAL(MChiI)(2) = slha_io.read_entry("MASS", 1000083);
   PHYSICAL(MChiI)(3) = slha_io.read_entry("MASS", 1000084);
   PHYSICAL(MChiI)(4) = slha_io.read_entry("MASS", 1000085);
   PHYSICAL(MChiI)(5) = slha_io.read_entry("MASS", 1000086);
   PHYSICAL(MChiI)(6) = slha_io.read_entry("MASS", 1000087);
   PHYSICAL(MSHp0)(0) = slha_io.read_entry("MASS", 92);
   PHYSICAL(MSHp0)(1) = slha_io.read_entry("MASS", 94);
   PHYSICAL(MSHpp)(0) = slha_io.read_entry("MASS", 91);
   PHYSICAL(MSHpp)(1) = slha_io.read_entry("MASS", 93);
   PHYSICAL(MChiP)(0) = slha_io.read_entry("MASS", 1000092);
   PHYSICAL(MChiP)(1) = slha_io.read_entry("MASS", 1000094);
   PHYSICAL(MVWm) = slha_io.read_entry("MASS", 24);

}

template <class T>
void CNE6SSM_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const CNE6SSM_slha<T>& model,
   const QedQcd& qedqcd, const CNE6SSM_scales& scales)
{
   CNE6SSM_slha_io slha_io;
   const CNE6SSM_input_parameters& input = model.get_input();
   const Problems<CNE6SSM_info::NUMBER_OF_PARTICLES>& problems
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
SLHAea::Coll CNE6SSM_slha_io::fill_slhaea(
   const CNE6SSM_slha<T>& model, const QedQcd& qedqcd)
{
   CNE6SSM_scales scales;

   return fill_slhaea(model, qedqcd, scales);
}

template <class T>
SLHAea::Coll CNE6SSM_slha_io::fill_slhaea(
   const CNE6SSM_slha<T>& model, const QedQcd& qedqcd,
   const CNE6SSM_scales& scales)
{
   SLHAea::Coll slhaea;
   CNE6SSM_slha_io::fill_slhaea(slhaea, model, qedqcd, scales);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class T>
void CNE6SSM_slha_io::set_model_parameters(const CNE6SSM_slha<T>& model)
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
   slha_io.set_block("Te", MODELPARAMETER(TYe_slha), "TYe", model.get_scale());
   slha_io.set_block("Td", MODELPARAMETER(TYd_slha), "TYd", model.get_scale());
   slha_io.set_block("Tu", MODELPARAMETER(TYu_slha), "TYu", model.get_scale());
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(30, (MODELPARAMETER(XiF)), "XiF")
            << FORMAT_ELEMENT(31, (MODELPARAMETER(MuPhi)), "MuPhi")
            << FORMAT_ELEMENT(32, (MODELPARAMETER(KappaPr)), "KappaPr")
            << FORMAT_ELEMENT(33, (MODELPARAMETER(Sigmax)), "Sigmax")
            << FORMAT_ELEMENT(34, (MODELPARAMETER(LXiF)), "LXiF")
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
            << FORMAT_ELEMENT(28, (MODELPARAMETER(TKappaPr)), "TKappaPr")
            << FORMAT_ELEMENT(29, (MODELPARAMETER(TSigmax)), "TSigmax")
            << FORMAT_ELEMENT(43, (MODELPARAMETER(TSigmaL)), "TSigmaL")
            << FORMAT_ELEMENT(30, (MODELPARAMETER(BMuPhi)), "BMuPhi")
            << FORMAT_ELEMENT(11, (MODELPARAMETER(vs)), "vs")
            << FORMAT_ELEMENT(12, (MODELPARAMETER(vsb)), "vsb")
            << FORMAT_ELEMENT(13, (MODELPARAMETER(vphi)), "vphi")
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Lambdax)), "Lambdax")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(TLambdax)), "TLambdax")
            << FORMAT_ELEMENT(0, (MODELPARAMETER(MuPr)), "MuPr")
            << FORMAT_ELEMENT(101, (MODELPARAMETER(BMuPr)), "BMuPr")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("ESIXGDYUK", MODELPARAMETER(gD), "gD", model.get_scale());
   slha_io.set_block("ESIXFUYUK", MODELPARAMETER(fu), "fu", model.get_scale());
   slha_io.set_block("ESIXFDYUK", MODELPARAMETER(fd), "fd", model.get_scale());
   slha_io.set_block("ESIXTHETRI", MODELPARAMETER(ThE), "ThE", model.get_scale());
   slha_io.set_block("ESIXTGDTRI", MODELPARAMETER(TgD), "TgD", model.get_scale());
   slha_io.set_block("ESIXTFUTRI", MODELPARAMETER(Tfu), "Tfu", model.get_scale());
   slha_io.set_block("ESIXTFDTRI", MODELPARAMETER(Tfd), "Tfd", model.get_scale());
   slha_io.set_block("MSQ2", MODELPARAMETER(mq2_slha), "mq2", model.get_scale());
   slha_io.set_block("MSE2", MODELPARAMETER(me2_slha), "me2", model.get_scale());
   slha_io.set_block("MSL2", MODELPARAMETER(ml2_slha), "ml2", model.get_scale());
   slha_io.set_block("MSU2", MODELPARAMETER(mu2_slha), "mu2", model.get_scale());
   slha_io.set_block("MSD2", MODELPARAMETER(md2_slha), "md2", model.get_scale());
   {
      std::ostringstream block;
      block << "Block MSOFT Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(21, (MODELPARAMETER(mHd2)), "mHd2")
            << FORMAT_ELEMENT(22, (MODELPARAMETER(mHu2)), "mHu2")
            << FORMAT_ELEMENT(23, (MODELPARAMETER(ms2)), "ms2")
            << FORMAT_ELEMENT(24, (MODELPARAMETER(msbar2)), "msbar2")
            << FORMAT_ELEMENT(25, (MODELPARAMETER(mphi2)), "mphi2")
            << FORMAT_ELEMENT(26, (MODELPARAMETER(mHp2)), "mHp2")
            << FORMAT_ELEMENT(27, (MODELPARAMETER(mHpbar2)), "mHpbar2")
            << FORMAT_ELEMENT(1, (MODELPARAMETER(MassB)), "MassB")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(MassWB)), "MassWB")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(MassG)), "MassG")
            << FORMAT_ELEMENT(4, (MODELPARAMETER(MassBp)), "MassBp")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("mX2", MODELPARAMETER(mDx2), "mDx2", model.get_scale());
   slha_io.set_block("mXBar2", MODELPARAMETER(mDxbar2), "mDxbar2", model.get_scale());
   slha_io.set_block("ESIXKAPPA", MODELPARAMETER(Kappa), "Kappa", model.get_scale());
   slha_io.set_block("ESIXTKAPPA", MODELPARAMETER(TKappa), "TKappa", model.get_scale());
   slha_io.set_block("ESIXLAMBDA", MODELPARAMETER(Lambda12), "Lambda12", model.get_scale());
   slha_io.set_block("ESIXTLAMBDA", MODELPARAMETER(TLambda12), "TLambda12", model.get_scale());

}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 */
template <class T>
void CNE6SSM_slha_io::set_extra(
   const CNE6SSM_slha<T>& model, const CNE6SSM_scales& scales)
{
   const CNE6SSM_physical physical(model.get_physical_slha());


}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void CNE6SSM_slha_io::set_spectrum(const CNE6SSM<T>& model)
{
   const CNE6SSM_slha<T> model_slha(model);
   set_spectrum(model_slha);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class T>
void CNE6SSM_slha_io::set_spectrum(const CNE6SSM_slha<T>& model)
{
   const CNE6SSM_physical physical(model.get_physical_slha());
   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);

   if (slha_io.get_modsel().quark_flavour_violated)
      set_ckm(model.get_ckm_matrix(), model.get_scale());

   if (slha_io.get_modsel().lepton_flavour_violated)
      set_pmns(model.get_pmns_matrix(), model.get_scale());
}

} // namespace flexiblesusy

#endif
