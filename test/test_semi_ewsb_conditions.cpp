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

// tests that the rearranged EWSB conditions are equal to the 
// automatically generated expressions
BOOST_AUTO_TEST_CASE( test_rearranged_ewsb_eqs )
{
   CNE6SSM_semianalytic_input_parameters<Two_scale> input;

   // for s too large lack of numerical precision
   // causes this test to fail.
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

   const double automatic_ewsb_eq_hh_1 = get_automatic_ewsb_eq_hh_1(model);
   const double automatic_ewsb_eq_hh_2 = get_automatic_ewsb_eq_hh_2(model);
   const double automatic_ewsb_eq_hh_3 = get_automatic_ewsb_eq_hh_3(model);
   const double automatic_ewsb_eq_hh_4 = get_automatic_ewsb_eq_hh_4(model);
   const double automatic_ewsb_eq_hh_5 = get_automatic_ewsb_eq_hh_5(model);

   const double ewsb_eq_hh_1 = model.get_ewsb_eq_hh_1();
   const double ewsb_eq_hh_2 = model.get_ewsb_eq_hh_2();
   const double ewsb_eq_hh_3 = model.get_ewsb_eq_hh_3();
   const double ewsb_eq_hh_4 = model.get_ewsb_eq_hh_4();
   const double ewsb_eq_hh_5 = model.get_ewsb_eq_hh_5();

   const double prec = 1.0e-14;

   BOOST_CHECK_CLOSE(automatic_ewsb_eq_hh_1, ewsb_eq_hh_1, prec);
   BOOST_CHECK_CLOSE(automatic_ewsb_eq_hh_2, ewsb_eq_hh_2, prec);
   BOOST_CHECK_CLOSE(automatic_ewsb_eq_hh_3, ewsb_eq_hh_3, prec);
   BOOST_CHECK_CLOSE(automatic_ewsb_eq_hh_4, ewsb_eq_hh_4, prec);
   BOOST_CHECK_CLOSE(automatic_ewsb_eq_hh_5, ewsb_eq_hh_5, prec);
}

// tests that the parameters that zero the ordinary
// EWSB conditions satisfy the modified ones
BOOST_AUTO_TEST_CASE( test_tree_level_ewsb_eqs_consistent )
{
   const unsigned ewsb_loop_order = 0;
   CNE6SSM_semianalytic_input_parameters<Two_scale> input;

   // for s too large lack of numerical precision
   // causes this test to fail.
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
      = model.get_ewsb_tree_level_soft_masses();

   model.set_mHd2(tree_level_soft_masses(0));
   model.set_mHu2(tree_level_soft_masses(1));
   model.set_ms2(tree_level_soft_masses(2));
   model.set_msbar2(tree_level_soft_masses(3));
   model.set_mphi2(tree_level_soft_masses(4));

   // note
   std::cout << "ewsb_eq_hh_1 = " << model.get_ewsb_eq_hh_1() << '\n';
   std::cout << "ewsb_eq_hh_2 = " << model.get_ewsb_eq_hh_2() << '\n';
   std::cout << "ewsb_eq_hh_3 = " << model.get_ewsb_eq_hh_3() << '\n';
   std::cout << "ewsb_eq_hh_4 = " << model.get_ewsb_eq_hh_4() << '\n';
   std::cout << "ewsb_eq_hh_5 = " << model.get_ewsb_eq_hh_5() << '\n';

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
      std::cout << "Abs(ewsb_eqs[0]) = " << residual << '\n';
   for (unsigned i = 1; i < model.number_of_tadpole_equations; ++i) {
      std::cout << "Abs(ewsb_eqs[" << i << "]) = " << Abs(ewsb_eqs[i]) << '\n';
      residual += Abs(ewsb_eqs[i]);
   }
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
   std::cout << "Abs(ewsb_eqs[0]) = " << residual << '\n';
   for (unsigned i = 1; i < model.number_of_tadpole_equations; ++i) {
      std::cout << "Abs(ewsb_eqs[" << i << "]) = " << Abs(ewsb_eqs[i]) << '\n';
      residual += Abs(ewsb_eqs[i]);
   }
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
   std::cout << "Abs(ewsb_eqs[0]) = " << residual << '\n';
   for (unsigned i = 1; i < model.number_of_tadpole_equations; ++i) {
      std::cout << "Abs(ewsb_eqs[" << i << "]) = " << Abs(ewsb_eqs[i]) << '\n';
      residual += Abs(ewsb_eqs[i]);
   }
   BOOST_CHECK_LT(residual, ewsb_iteration_precision);
}
