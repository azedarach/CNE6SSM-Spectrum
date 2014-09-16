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
BOOST_AUTO_TEST_CASE( test_default_ewsb )
{
   CNE6SSM<Two_scale> model = CNE6SSM<Two_scale>();

   BOOST_CHECK_EQUAL(model.get_ewsb_eq_hh_1(), 0.0);
   BOOST_CHECK_EQUAL(model.get_ewsb_eq_hh_2(), 0.0);
   BOOST_CHECK_EQUAL(model.get_ewsb_eq_hh_3(), 0.0);
   BOOST_CHECK_EQUAL(model.get_ewsb_eq_hh_4(), 0.0);
   BOOST_CHECK_EQUAL(model.get_ewsb_eq_hh_5(), 0.0);
}

// Test that the custom EWSB conditions give
// expected values for a set of parameters
BOOST_AUTO_TEST_CASE( test_custom_ewsb )
{

}
