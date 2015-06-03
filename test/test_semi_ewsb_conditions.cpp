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

// tests that the parameters that zero the ordinary
// EWSB conditions satisfy the modified ones
BOOST_AUTO_TEST_CASE( test_tree_level_ewsb_eqs_consistent )
{
   const unsigned ewsb_loop_order = 0;
   CNE6SSM_semianalytic_input_parameters<Two_scale> input;

   // for s too large lack of numerical precision
   // causes this test to fail.
   input.sInput = 5000.0;
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
