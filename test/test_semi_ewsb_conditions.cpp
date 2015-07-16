// ====================================================================
// Test suite for implementation of semianalytic EWSB conditions
// ====================================================================

#include <iomanip>
#include <limits>
#include <random>
#include <chrono>
#include <Eigen/Core>

#include "wrappers.hpp"
#include "CNE6SSM_semi_two_scale_model.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_semi_ewsb_conditions

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace std;
using namespace flexiblesusy;

void set_test_model_parameters(CNE6SSM_semianalytic<Two_scale>& model)
{
   model.set_Yd(0, 0, 0.);
   model.set_Yd(1, 1, 0.);
   model.set_Yd(2, 2, 0.1);

   model.set_Ye(0, 0, 0.);
   model.set_Ye(1, 1, 0.);
   model.set_Ye(2, 2, 0.1);

   model.set_SigmaL(0.47);
   model.set_KappaPr(0.016);
   model.set_Sigmax(0.074);

   model.set_Kappa(0, 0, 0.42);
   model.set_Kappa(1, 1, 0.42);
   model.set_Kappa(2, 2, 0.42);

   model.set_Lambda12(0, 0, 0.57);
   model.set_Lambda12(1, 1, 0.57);

   model.set_Lambdax(0.11);

   model.set_Yu(0, 0, 0.);
   model.set_Yu(1, 1, 0.);
   model.set_Yu(2, 2, 0.87);

   model.set_MuPr(17000.);
   model.set_MuPhi(0.);
   model.set_BMuPr(10000.);
   model.set_BMuPhi(0.);

   model.set_g1(0.463);
   model.set_g2(0.648);
   model.set_g3(1.161);
   model.set_g1p(0.442);

   model.set_vd(23.63);
   model.set_vu(236.28);

   model.set_TSigmax(513.2);
   model.set_TLambdax(124.2);
   model.set_TKappaPr(320.1);

   model.set_QS(model.get_input().QSInput);
}

// automatically generated versions of EWSB conditions
double get_automatic_ewsb_eq_hh_1(const CNE6SSM_mass_eigenstates& model)
{
   const double Sigmax = model.get_Sigmax();
   const double Lambdax = model.get_Lambdax();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double TLambdax = model.get_TLambdax();
   const double mHd2 = model.get_mHd2();

   double result = Re(mHd2*vd - 0.35355339059327373*vs*vu*Conj(TLambdax) + 0.25
      *vphi*vsb*vu*Conj(Sigmax)*Lambdax + 0.25*vphi*vsb*vu*Conj(Lambdax)*Sigmax +
      0.075*Power(vd,3)*Sqr(g1) + 0.1125*Power(vd,3)*Sqr(g1p) + 0.125*Power(vd,3)*
      Sqr(g2) + 0.5*vd*AbsSqr(Lambdax)*Sqr(vs) - 0.0375*QS*vd*Sqr(g1p)*Sqr(vs) +
      0.0375*QS*vd*Sqr(g1p)*Sqr(vsb) + 0.5*vd*AbsSqr(Lambdax)*Sqr(vu) - 0.075*vd*
      Sqr(g1)*Sqr(vu) + 0.075*vd*Sqr(g1p)*Sqr(vu) - 0.125*vd*Sqr(g2)*Sqr(vu) -
      0.35355339059327373*vs*vu*TLambdax);

   return result;
}

double get_automatic_ewsb_eq_hh_2(const CNE6SSM_mass_eigenstates& model)
{
   const double Sigmax = model.get_Sigmax();
   const double Lambdax = model.get_Lambdax();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double TLambdax = model.get_TLambdax();
   const double mHu2 = model.get_mHu2();

   double result = Re(mHu2*vu - 0.35355339059327373*vd*vs*Conj(TLambdax) + 0.25
      *vd*vphi*vsb*Conj(Sigmax)*Lambdax + 0.25*vd*vphi*vsb*Conj(Lambdax)*Sigmax +
      0.075*Power(vu,3)*Sqr(g1) + 0.05*Power(vu,3)*Sqr(g1p) + 0.125*Power(vu,3)*
      Sqr(g2) + 0.5*vu*AbsSqr(Lambdax)*Sqr(vd) - 0.075*vu*Sqr(g1)*Sqr(vd) + 0.075*
      vu*Sqr(g1p)*Sqr(vd) - 0.125*vu*Sqr(g2)*Sqr(vd) + 0.5*vu*AbsSqr(Lambdax)*Sqr(
      vs) - 0.025*QS*vu*Sqr(g1p)*Sqr(vs) + 0.025*QS*vu*Sqr(g1p)*Sqr(vsb) -
      0.35355339059327373*vd*vs*TLambdax);

   return result;
}

double get_automatic_ewsb_eq_hh_3(const CNE6SSM_mass_eigenstates& model)
{
   const double Sigmax = model.get_Sigmax();
   const double KappaPr = model.get_KappaPr();
   const double MuPhi = model.get_MuPhi();
   const double Lambdax = model.get_Lambdax();
   const double XiF = model.get_XiF();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double TSigmax = model.get_TSigmax();
   const double TLambdax = model.get_TLambdax();
   const double ms2 = model.get_ms2();

   double result = Re(ms2*vs - 0.35355339059327373*MuPhi*vphi*vsb*Conj(Sigmax)
      - 0.35355339059327373*vd*vu*Conj(TLambdax) - 0.35355339059327373*vphi*vsb*
      Conj(TSigmax) - 0.5*vsb*Conj(Sigmax)*XiF - 0.35355339059327373*vphi*vsb*Conj
      (MuPhi)*Sigmax - 0.5*vsb*Conj(XiF)*Sigmax + 0.0125*Power(vs,3)*Sqr(g1p)*Sqr(
      QS) + 0.5*vs*AbsSqr(Lambdax)*Sqr(vd) - 0.0375*QS*vs*Sqr(g1p)*Sqr(vd) + 0.5*
      vs*AbsSqr(Sigmax)*Sqr(vphi) - 0.25*vsb*Conj(Sigmax)*KappaPr*Sqr(vphi) - 0.25
      *vsb*Conj(KappaPr)*Sigmax*Sqr(vphi) + 0.5*vs*AbsSqr(Sigmax)*Sqr(vsb) -
      0.0125*vs*Sqr(g1p)*Sqr(QS)*Sqr(vsb) + 0.5*vs*AbsSqr(Lambdax)*Sqr(vu) - 0.025
      *QS*vs*Sqr(g1p)*Sqr(vu) - 0.35355339059327373*vd*vu*TLambdax -
      0.35355339059327373*vphi*vsb*TSigmax);

   return result;
}

double get_automatic_ewsb_eq_hh_4(const CNE6SSM_mass_eigenstates& model)
{
   const double Sigmax = model.get_Sigmax();
   const double KappaPr = model.get_KappaPr();
   const double MuPhi = model.get_MuPhi();
   const double Lambdax = model.get_Lambdax();
   const double XiF = model.get_XiF();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double TSigmax = model.get_TSigmax();
   const double msbar2 = model.get_msbar2();

   double result = Re(msbar2*vsb - 0.35355339059327373*MuPhi*vphi*vs*Conj(
      Sigmax) - 0.35355339059327373*vphi*vs*Conj(TSigmax) + 0.25*vd*vphi*vu*Conj(
      Sigmax)*Lambdax - 0.5*vs*Conj(Sigmax)*XiF - 0.35355339059327373*vphi*vs*Conj
      (MuPhi)*Sigmax + 0.25*vd*vphi*vu*Conj(Lambdax)*Sigmax - 0.5*vs*Conj(XiF)*
      Sigmax + 0.0125*Power(vsb,3)*Sqr(g1p)*Sqr(QS) + 0.0375*QS*vsb*Sqr(g1p)*Sqr(
      vd) + 0.5*vsb*AbsSqr(Sigmax)*Sqr(vphi) - 0.25*vs*Conj(Sigmax)*KappaPr*Sqr(
      vphi) - 0.25*vs*Conj(KappaPr)*Sigmax*Sqr(vphi) + 0.5*vsb*AbsSqr(Sigmax)*Sqr(
      vs) - 0.0125*vsb*Sqr(g1p)*Sqr(QS)*Sqr(vs) + 0.025*QS*vsb*Sqr(g1p)*Sqr(vu) -
      0.35355339059327373*vphi*vs*TSigmax);

   return result;
}

double get_automatic_ewsb_eq_hh_5(const CNE6SSM_mass_eigenstates& model)
{
   const double Sigmax = model.get_Sigmax();
   const double KappaPr = model.get_KappaPr();
   const double MuPhi = model.get_MuPhi();
   const double Lambdax = model.get_Lambdax();
   const double XiF = model.get_XiF();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double TSigmax = model.get_TSigmax();
   const double TKappaPr = model.get_TKappaPr();
   const double BMuPhi = model.get_BMuPhi();
   const double LXiF = model.get_LXiF();
   const double mphi2 = model.get_mphi2();

   double result = Re(mphi2*vphi + vphi*AbsSqr(MuPhi) + Power(vphi,3)*AbsSqr(
      KappaPr) + 0.5*vphi*BMuPhi + 0.5*vphi*Conj(BMuPhi) + 0.7071067811865475*
      MuPhi*Conj(XiF) - 0.35355339059327373*MuPhi*vs*vsb*Conj(Sigmax) +
      0.7071067811865475*Conj(LXiF) - 0.35355339059327373*vs*vsb*Conj(TSigmax) +
      vphi*Conj(XiF)*KappaPr - 0.5*vphi*vs*vsb*Conj(Sigmax)*KappaPr + 0.25*vd*vsb*
      vu*Conj(Sigmax)*Lambdax + 0.7071067811865475*Conj(MuPhi)*XiF + vphi*Conj(
      KappaPr)*XiF - 0.35355339059327373*vs*vsb*Conj(MuPhi)*Sigmax - 0.5*vphi*vs*
      vsb*Conj(KappaPr)*Sigmax + 0.25*vd*vsb*vu*Conj(Lambdax)*Sigmax +
      0.7071067811865475*LXiF + 1.0606601717798212*MuPhi*Conj(KappaPr)*Sqr(vphi) +
      0.35355339059327373*Conj(TKappaPr)*Sqr(vphi) + 1.0606601717798212*Conj(
      MuPhi)*KappaPr*Sqr(vphi) + 0.5*vphi*AbsSqr(Sigmax)*Sqr(vs) + 0.5*vphi*AbsSqr
      (Sigmax)*Sqr(vsb) + 0.35355339059327373*Sqr(vphi)*TKappaPr -
      0.35355339059327373*vs*vsb*TSigmax);

   return result;
}

Eigen::Array<double,5,1> get_automatic_ewsb_tree_level_soft_masses(const CNE6SSM_mass_eigenstates& model)
{
   const double Sigmax = model.get_Sigmax();
   const double KappaPr = model.get_KappaPr();
   const double Lambdax = model.get_Lambdax();
   const double MuPhi = model.get_MuPhi();
   const double XiF = model.get_XiF();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double TSigmax = model.get_TSigmax();
   const double TKappaPr = model.get_TKappaPr();
   const double TLambdax = model.get_TLambdax();
   const double BMuPhi = model.get_BMuPhi();
   const double LXiF = model.get_LXiF();

   const double mHd2 = Re((0.0125*(28.284271247461902*vs*vu*Conj(TLambdax)
      - 20*vphi*vsb*vu*Conj(Sigmax)*Lambdax - 20*vphi*vsb*vu*Conj(Lambdax)*Sigmax
      - 6*Power(vd,3)*Sqr(g1) - 9*Power(vd,3)*Sqr(g1p) - 10*Power(vd,3)*Sqr(g2) -
      40*vd*AbsSqr(Lambdax)*Sqr(vs) + 3*QS*vd*Sqr(g1p)*Sqr(vs) - 3*QS*vd*Sqr(g1p)*
      Sqr(vsb) - 40*vd*AbsSqr(Lambdax)*Sqr(vu) + 6*vd*Sqr(g1)*Sqr(vu) - 6*vd*Sqr(
      g1p)*Sqr(vu) + 10*vd*Sqr(g2)*Sqr(vu) + 28.284271247461902*vs*vu*TLambdax))
      /vd);

   const double mHu2 = Re((0.025*(14.142135623730951*vd*vs*Conj(TLambdax) -
      10*vd*vphi*vsb*Conj(Sigmax)*Lambdax - 10*vd*vphi*vsb*Conj(Lambdax)*Sigmax -
      3*Power(vu,3)*Sqr(g1) - 2*Power(vu,3)*Sqr(g1p) - 5*Power(vu,3)*Sqr(g2) - 20
      *vu*AbsSqr(Lambdax)*Sqr(vd) + 3*vu*Sqr(g1)*Sqr(vd) - 3*vu*Sqr(g1p)*Sqr(vd) +
      5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(Lambdax)*Sqr(vs) + QS*vu*Sqr(g1p)*Sqr(
      vs) - QS*vu*Sqr(g1p)*Sqr(vsb) + 14.142135623730951*vd*vs*TLambdax))/vu);

   const double ms2 = Re((0.0125*(28.284271247461902*MuPhi*vphi*vsb*Conj(
      Sigmax) + 28.284271247461902*vd*vu*Conj(TLambdax) + 28.284271247461902*vphi*
      vsb*Conj(TSigmax) + 40*vsb*Conj(Sigmax)*XiF + 28.284271247461902*vphi*vsb*
      Conj(MuPhi)*Sigmax + 40*vsb*Conj(XiF)*Sigmax - Power(vs,3)*Sqr(g1p)*Sqr(QS)
      - 40*vs*AbsSqr(Lambdax)*Sqr(vd) + 3*QS*vs*Sqr(g1p)*Sqr(vd) - 40*vs*AbsSqr(
      Sigmax)*Sqr(vphi) + 20*vsb*Conj(Sigmax)*KappaPr*Sqr(vphi) + 20*vsb*Conj(
      KappaPr)*Sigmax*Sqr(vphi) - 40*vs*AbsSqr(Sigmax)*Sqr(vsb) + vs*Sqr(g1p)*Sqr(
      QS)*Sqr(vsb) - 40*vs*AbsSqr(Lambdax)*Sqr(vu) + 2*QS*vs*Sqr(g1p)*Sqr(vu) +
      28.284271247461902*vd*vu*TLambdax + 28.284271247461902*vphi*vsb*TSigmax))/vs
      );

   const double msbar2 = Re((0.0125*(28.284271247461902*MuPhi*vphi*vs*Conj(
      Sigmax) + 28.284271247461902*vphi*vs*Conj(TSigmax) - 20*vd*vphi*vu*Conj(
      Sigmax)*Lambdax + 40*vs*Conj(Sigmax)*XiF + 28.284271247461902*vphi*vs*Conj(
      MuPhi)*Sigmax - 20*vd*vphi*vu*Conj(Lambdax)*Sigmax + 40*vs*Conj(XiF)*Sigmax
      - Power(vsb,3)*Sqr(g1p)*Sqr(QS) - 3*QS*vsb*Sqr(g1p)*Sqr(vd) - 40*vsb*AbsSqr(
      Sigmax)*Sqr(vphi) + 20*vs*Conj(Sigmax)*KappaPr*Sqr(vphi) + 20*vs*Conj(
      KappaPr)*Sigmax*Sqr(vphi) - 40*vsb*AbsSqr(Sigmax)*Sqr(vs) + vsb*Sqr(g1p)*Sqr
      (QS)*Sqr(vs) - 2*QS*vsb*Sqr(g1p)*Sqr(vu) + 28.284271247461902*vphi*vs*
      TSigmax))/vsb);

   const double mphi2 = Re((0.25*(-4*vphi*AbsSqr(MuPhi) - 4*Power(vphi,3)*
      AbsSqr(KappaPr) - 2*vphi*BMuPhi - 2*vphi*Conj(BMuPhi) - 2.8284271247461903*
      MuPhi*Conj(XiF) + 1.4142135623730951*MuPhi*vs*vsb*Conj(Sigmax) -
      2.8284271247461903*Conj(LXiF) + 1.4142135623730951*vs*vsb*Conj(TSigmax) - 4*
      vphi*Conj(XiF)*KappaPr + 2*vphi*vs*vsb*Conj(Sigmax)*KappaPr - vd*vsb*vu*Conj
      (Sigmax)*Lambdax - 2.8284271247461903*Conj(MuPhi)*XiF - 4*vphi*Conj(KappaPr)
      *XiF + 1.4142135623730951*vs*vsb*Conj(MuPhi)*Sigmax + 2*vphi*vs*vsb*Conj(
      KappaPr)*Sigmax - vd*vsb*vu*Conj(Lambdax)*Sigmax - 2.8284271247461903*LXiF -
      4.242640687119286*MuPhi*Conj(KappaPr)*Sqr(vphi) - 1.4142135623730951*Conj(
      TKappaPr)*Sqr(vphi) - 4.242640687119286*Conj(MuPhi)*KappaPr*Sqr(vphi) - 2*
      vphi*AbsSqr(Sigmax)*Sqr(vs) - 2*vphi*AbsSqr(Sigmax)*Sqr(vsb) -
      1.4142135623730951*Sqr(vphi)*TKappaPr + 1.4142135623730951*vs*vsb*TSigmax))
      /vphi);

   Eigen::Array<double,5,1> masses;
   masses << mHd2, mHu2, ms2, msbar2, mphi2;

   return masses;
}

BOOST_AUTO_TEST_CASE( test_soft_higgs_masses )
{
   CNE6SSM_semianalytic_input_parameters<Two_scale> input;

   input.sInput = 650000.0;
   input.QSInput = 0.;
   input.TanBeta = 10.;
   input.m12 = 1600;
   input.Azero = 0.;

   CNE6SSM_semianalytic<Two_scale> model(input);

   set_test_model_parameters(model);
   model.set_scale(1000.);

   const double TanTheta = 0.95;
   model.set_vs(input.sInput * Cos(ArcTan(TanTheta)));
   model.set_vsb(input.sInput * Sin(ArcTan(TanTheta)));
   model.set_vphi(23013.5);
   model.set_XiF(1.0e9);
   model.set_LXiF(1.0e10);

   Eigen::Array<double,5,1> tree_level_soft_masses
      = model.get_ewsb_tree_level_soft_masses();

   Eigen::Array<double,5,1> tree_level_soft_masses_old
      = get_automatic_ewsb_tree_level_soft_masses(model);

   const double prec = 1.0e-10;

   BOOST_CHECK_CLOSE(tree_level_soft_masses(0),
                     tree_level_soft_masses_old(0), prec);
   BOOST_CHECK_CLOSE(tree_level_soft_masses(1),
                     tree_level_soft_masses_old(1), prec);
   BOOST_CHECK_CLOSE(tree_level_soft_masses(2),
                     tree_level_soft_masses_old(2), prec);
   BOOST_CHECK_CLOSE(tree_level_soft_masses(3),
                     tree_level_soft_masses_old(3), prec);
   BOOST_CHECK_CLOSE(tree_level_soft_masses(4),
                     tree_level_soft_masses_old(4), prec);
}

// tests that the rearranged EWSB conditions are solved
// by the same parameters as the automatically generated expressions
BOOST_AUTO_TEST_CASE( test_rearranged_ewsb_eqs )
{
   CNE6SSM_semianalytic_input_parameters<Two_scale> input;

   input.sInput = 650000.0;
   input.QSInput = 0.;
   input.TanBeta = 10.;
   input.m12 = 1600;
   input.Azero = 0.;

   CNE6SSM_semianalytic<Two_scale> model(input);

   set_test_model_parameters(model);
   model.set_scale(1000.);

   const double TanTheta = 0.95;
   model.set_vs(input.sInput * Cos(ArcTan(TanTheta)));
   model.set_vsb(input.sInput * Sin(ArcTan(TanTheta)));
   model.set_vphi(23013.5);
   model.set_XiF(1.0e9);
   model.set_LXiF(1.0e10);

   Eigen::Array<double,5,1> tree_level_soft_masses
      = get_automatic_ewsb_tree_level_soft_masses(model);

   model.set_mHd2(tree_level_soft_masses(0));
   model.set_mHu2(tree_level_soft_masses(1));
   model.set_ms2(tree_level_soft_masses(2));
   model.set_msbar2(tree_level_soft_masses(3));
   model.set_mphi2(tree_level_soft_masses(4));

   const double ewsb_eq_hh_1 = model.get_ewsb_eq_hh_1();
   const double ewsb_eq_hh_2 = model.get_ewsb_eq_hh_2();
   const double ewsb_eq_hh_3 = model.get_ewsb_eq_hh_3();
   const double ewsb_eq_hh_4 = model.get_ewsb_eq_hh_4();
   const double ewsb_eq_hh_5 = model.get_ewsb_eq_hh_5();

   const double prec = 1.0e-5;

   BOOST_CHECK_LT(Abs(ewsb_eq_hh_1), prec);
   BOOST_CHECK_LT(Abs(ewsb_eq_hh_2), prec);
   BOOST_CHECK_LT(Abs(ewsb_eq_hh_3), prec);
   BOOST_CHECK_LT(Abs(ewsb_eq_hh_4), prec);
   BOOST_CHECK_LT(Abs(ewsb_eq_hh_5), prec);
}

// tests that the parameters that zero the ordinary
// EWSB conditions satisfy the modified ones
// Note: assumes the rearranged versions match the
// generated versions, i.e. this is only valid if
// all of the tests above pass
BOOST_AUTO_TEST_CASE( test_tree_level_ewsb_eqs_consistent )
{
   const unsigned ewsb_loop_order = 0;
   CNE6SSM_semianalytic_input_parameters<Two_scale> input;

   input.sInput = 650000.0;
   input.QSInput = 0.;
   input.TanBeta = 10.;
   input.m12 = 1600;
   input.Azero = 0.;

   CNE6SSM_semianalytic<Two_scale> model(input);

   set_test_model_parameters(model);
   model.set_scale(1000.);
   model.set_ewsb_loop_order(ewsb_loop_order);

   const double TanTheta = 0.95;
   model.set_vs(input.sInput * Cos(ArcTan(TanTheta)));
   model.set_vsb(input.sInput * Sin(ArcTan(TanTheta)));
   model.set_vphi(23013.5);
   model.set_XiF(1.0e9);
   model.set_LXiF(1.0e10);

   Eigen::Array<double,5,1> tree_level_soft_masses
      = get_automatic_ewsb_tree_level_soft_masses(model);

   model.set_mHd2(tree_level_soft_masses(0));
   model.set_mHu2(tree_level_soft_masses(1));
   model.set_ms2(tree_level_soft_masses(2));
   model.set_msbar2(tree_level_soft_masses(3));
   model.set_mphi2(tree_level_soft_masses(4));

   BOOST_REQUIRE(Abs(model.get_ewsb_eq_hh_1()) < 1.0e-3);
   BOOST_REQUIRE(Abs(model.get_ewsb_eq_hh_2()) < 1.0e-3);
   BOOST_REQUIRE(Abs(model.get_ewsb_eq_hh_3()) < 1.0e-3);
   BOOST_REQUIRE(Abs(model.get_ewsb_eq_hh_4()) < 1.0e-3);
   BOOST_REQUIRE(Abs(model.get_ewsb_eq_hh_5()) < 1.0e-3);

   BOOST_CHECK_LT(Abs(model.get_tree_level_ewsb_eq_hh_1()), 1.0e-3);
   BOOST_CHECK_LT(Abs(model.get_tree_level_ewsb_eq_hh_2()), 1.0e-3);
   BOOST_CHECK_LT(Abs(model.get_tree_level_ewsb_eq_hh_3()), 1.0e-3);
   BOOST_CHECK_LT(Abs(model.get_tree_level_ewsb_eq_hh_4()), 1.0e-3);
   BOOST_CHECK_LT(Abs(model.get_tree_level_ewsb_eq_hh_5()), 1.0e-3);
}

BOOST_AUTO_TEST_CASE( test_semi_ewsb_tree_level )
{
   CNE6SSM_semianalytic_input_parameters<Two_scale> input;

   input.sInput = 40000.;
   input.QSInput = 5.;
   input.TanBeta = 10.;
   input.m12 = 2000.;
   input.Azero = 0.;

   CNE6SSM_semianalytic<Two_scale> model(input);

   set_test_model_parameters(model);
   model.set_scale(1000.);

   const double ewsb_iteration_precision = 1.0e-3;
   model.set_ewsb_iteration_precision(ewsb_iteration_precision);

   const double beta_loops = 2;
   model.set_loops(beta_loops);
   const unsigned ewsb_loops = 0;
   model.set_ewsb_loop_order(ewsb_loops);

   const double high_scale_guess = 2.0e16;
   model.calculate_coefficients(high_scale_guess);

   model.solve_ewsb();

   double ewsb_eqs[model.number_of_tadpole_equations];

   model.tadpole_equations(ewsb_eqs);

   double residual = Abs(ewsb_eqs[0]);

   for (unsigned i = 1; i < model.number_of_tadpole_equations; ++i)
      residual += Abs(ewsb_eqs[i]);

   BOOST_CHECK_LT(residual, ewsb_iteration_precision);

}

BOOST_AUTO_TEST_CASE( test_semi_ewsb_one_loop )
{
   CNE6SSM_semianalytic_input_parameters<Two_scale> input;

   input.sInput = 40000.;
   input.QSInput = 5.;
   input.TanBeta = 10.;
   input.m12 = 5000.;
   input.Azero = 0.;

   CNE6SSM_semianalytic<Two_scale> model(input);

   set_test_model_parameters(model);
   model.set_scale(1000.);

   const double ewsb_iteration_precision = 1.0e-3;
   model.set_ewsb_iteration_precision(ewsb_iteration_precision);

   const double beta_loops = 2;
   model.set_loops(beta_loops);
   const unsigned ewsb_loops = 1;
   model.set_ewsb_loop_order(ewsb_loops);

   const double high_scale_guess = 2.0e16;
   model.calculate_coefficients(high_scale_guess);

   model.solve_ewsb();

   double ewsb_eqs[model.number_of_tadpole_equations];

   model.tadpole_equations(ewsb_eqs);

   double residual = Abs(ewsb_eqs[0]);

   for (unsigned i = 1; i < model.number_of_tadpole_equations; ++i)
      residual += Abs(ewsb_eqs[i]);

   BOOST_CHECK_LT(residual, ewsb_iteration_precision);
}

BOOST_AUTO_TEST_CASE( test_semi_ewsb_two_loop )
{
   CNE6SSM_semianalytic_input_parameters<Two_scale> input;

   input.sInput = 40000.;
   input.QSInput = 5.;
   input.TanBeta = 10.;
   input.m12 = 2000.;
   input.Azero = 0.;

   CNE6SSM_semianalytic<Two_scale> model(input);

   set_test_model_parameters(model);
   model.set_scale(1000.);

   const double ewsb_iteration_precision = 1.0e-3;
   model.set_ewsb_iteration_precision(ewsb_iteration_precision);

   const double beta_loops = 2;
   model.set_loops(beta_loops);
   const unsigned ewsb_loops = 2;
   model.set_ewsb_loop_order(ewsb_loops);

   const double high_scale_guess = 2.0e16;
   model.calculate_coefficients(high_scale_guess);

   model.solve_ewsb();

   double ewsb_eqs[model.number_of_tadpole_equations];

   model.tadpole_equations(ewsb_eqs);

   double residual = Abs(ewsb_eqs[0]);

   for (unsigned i = 1; i < model.number_of_tadpole_equations; ++i)
      residual += Abs(ewsb_eqs[i]);

   BOOST_CHECK_LT(residual, ewsb_iteration_precision);
}
