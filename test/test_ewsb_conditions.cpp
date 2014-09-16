// ====================================================================
// Test suite for implementation of EWSB conditions
// ====================================================================

#include "CNE6SSM_two_scale_model.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_ewsb_conditions

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace std;
using namespace flexiblesusy;

// Test that the default EWSB conditions give
// expected values for a set of parameters
BOOST_AUTO_TEST_CASE( test_default_ewsb_tree_level )
{
   CNE6SSM<Two_scale> default_model = CNE6SSM<Two_scale>();

   BOOST_CHECK_EQUAL(default_model.get_ewsb_eq_hh_1(), 0.);
   BOOST_CHECK_EQUAL(default_model.get_ewsb_eq_hh_2(), 0.);
   BOOST_CHECK_EQUAL(default_model.get_ewsb_eq_hh_3(), 0.);
   BOOST_CHECK_EQUAL(default_model.get_ewsb_eq_hh_4(), 0.);
   BOOST_CHECK_EQUAL(default_model.get_ewsb_eq_hh_5(), 0.);

   CNE6SSM_input_parameters inputs;
   inputs.QS = 5.;
   default_model.set_input_parameters(inputs);

   default_model.set_vd(24.4779);
   default_model.set_vu(244.779);
   default_model.set_vs(5381.382);
   default_model.set_vsb(5919.520587);
   default_model.set_vphi(8000.);
   default_model.set_g1(0.46);
   default_model.set_g1p(0.47);
   default_model.set_g2(0.61);
   default_model.set_Lambdax(0.01);
   default_model.set_Sigmax(0.12);
   default_model.set_TLambdax(1100.);
   default_model.set_TSigmax(540.);
   default_model.set_MuPhi(100.);
   default_model.set_BMuPhi(150.);
   default_model.set_KappaPr(0.5);
   default_model.set_TKappaPr(3000.);
   default_model.set_XiF(0.4);
   default_model.set_LXiF(10.);
   default_model.set_mHd2(2.3e6);
   default_model.set_mHu2(-1.3e6);
   default_model.set_ms2(-1.9e6);
   default_model.set_msbar2(2.1e6);
   default_model.set_mphi2(1.2e5);

   BOOST_CHECK_CLOSE(default_model.get_ewsb_eq_hh_1(), -9.551900152e8, 1.0e-7);
   BOOST_CHECK_CLOSE(default_model.get_ewsb_eq_hh_2(), -3.774462196e8, 1.0e-7);
   BOOST_CHECK_CLOSE(default_model.get_ewsb_eq_hh_3(), -3.850971629e10, 1.0e-7);
   BOOST_CHECK_CLOSE(default_model.get_ewsb_eq_hh_4(), -8.248055767e9, 1.0e-7);
   BOOST_CHECK_CLOSE(default_model.get_ewsb_eq_hh_5(), 2.475560387e11, 1.0e-7);

}

// Test that the custom EWSB conditions give
// expected values for a set of parameters
BOOST_AUTO_TEST_CASE( test_custom_ewsb_tree_level )
{
   CNE6SSM<Two_scale> custom_model = CNE6SSM<Two_scale>();
   
   //BOOST_CHECK_EQUAL(custom_model.get_alternate_ewsb_eq_hh_1(), 0.);
   //BOOST_CHECK_EQUAL(custom_model.get_alternate_ewsb_eq_hh_2(), 0.);
   //BOOST_CHECK_EQUAL(custom_model.get_alternate_ewsb_eq_hh_3(), 0.);
   //BOOST_CHECK_EQUAL(custom_model.get_alternate_ewsb_eq_hh_4(), 0.);
   //BOOST_CHECK_EQUAL(custom_model.get_alternate_ewsb_eq_hh_5(), 0.);


}
