// ====================================================================
// Test suite for implementation of rotated CP-even Higgs mass matrix
// ====================================================================

#include <iomanip>
#include <limits>
#include <random>
#include <chrono>
#include <Eigen/Core>

#include "wrappers.hpp"
#include "CNE6SSM_semi_two_scale_model.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_rotated_mass_matrix_hh

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
   model.set_MuPhi(125.3);
   model.set_BMuPr(10000.);
   model.set_BMuPhi(543.2);

   model.set_XiF(1893.2);
   model.set_LXiF(4319.6);

   model.set_g1(0.463);
   model.set_g2(0.648);
   model.set_g3(1.161);
   model.set_g1p(0.442);

   model.set_vd(23.63);
   model.set_vu(236.28);
   model.set_vs(14235.3);
   model.set_vsb(14673.7);
   model.set_vphi(7231.3);

   model.set_TLambdax(5.4193);
   model.set_TSigmax(56.3994);
   model.set_TKappaPr(32.391);

   model.set_QS(model.get_input().QSInput);
}

// tests that the resulting mass matrix matches
// the expected value
BOOST_AUTO_TEST_CASE( test_rotated_mass_matrix_hh )
{
   CNE6SSM_semianalytic_input_parameters<Two_scale> input;

   input.sInput = 5000.0;
   input.QSInput = 5.;
   input.TanBeta = 10.;
   input.m12 = 1600;
   input.Azero = 0.;

   CNE6SSM_semianalytic<Two_scale> model(input);
   set_test_model_parameters(model);
   model.set_scale(1000.);

   Eigen::Matrix<double,5,5> mass_matrix_hh
      = model.get_rotated_mass_matrix_hh();

   const double expected_Mh00 = 5.0624234833915815e7;
   const double expected_Mh01 = -12427.79297188903;
   const double expected_Mh02 = -24564.271402040395;
   const double expected_Mh03 = 17276.694757626472;
   const double expected_Mh04 = -210917.8462479296;
   const double expected_Mh11 = 1.14399600387621e6;
   const double expected_Mh12 = -315197.5112848233;
   const double expected_Mh13 = -4129.728288077248;
   const double expected_Mh14 = 30191.35836227309;
   const double expected_Mh22 = 1.5644404007465127e6;
   const double expected_Mh23 = -13900.62758618708;
   const double expected_Mh24 = 2808.45009751649;
   const double expected_Mh33 = -3.810296946643915e6;
   const double expected_Mh34 = 1379.9336124121455;
   const double expected_Mh44 = 8555.027848943877;

   BOOST_CHECK_CLOSE(mass_matrix_hh(0,0), expected_Mh00, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(0,1), expected_Mh01, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(0,2), expected_Mh02, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(0,3), expected_Mh03, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(0,4), expected_Mh04, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(1,1), expected_Mh11, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(1,2), expected_Mh12, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(1,3), expected_Mh13, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(1,4), expected_Mh14, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(2,2), expected_Mh22, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(2,3), expected_Mh23, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(2,4), expected_Mh24, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(3,3), expected_Mh33, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(3,4), expected_Mh34, 1.0e-10);
   BOOST_CHECK_CLOSE(mass_matrix_hh(4,4), expected_Mh44, 1.0e-10);

   BOOST_CHECK_EQUAL(mass_matrix_hh(1,0), mass_matrix_hh(0,1));
   BOOST_CHECK_EQUAL(mass_matrix_hh(2,0), mass_matrix_hh(0,2));
   BOOST_CHECK_EQUAL(mass_matrix_hh(3,0), mass_matrix_hh(0,3));
   BOOST_CHECK_EQUAL(mass_matrix_hh(4,0), mass_matrix_hh(0,4));
   BOOST_CHECK_EQUAL(mass_matrix_hh(2,1), mass_matrix_hh(1,2));
   BOOST_CHECK_EQUAL(mass_matrix_hh(3,1), mass_matrix_hh(1,3));
   BOOST_CHECK_EQUAL(mass_matrix_hh(4,1), mass_matrix_hh(1,4));
   BOOST_CHECK_EQUAL(mass_matrix_hh(3,2), mass_matrix_hh(2,3));
   BOOST_CHECK_EQUAL(mass_matrix_hh(4,2), mass_matrix_hh(2,4));
   BOOST_CHECK_EQUAL(mass_matrix_hh(4,3), mass_matrix_hh(3,4));
}
