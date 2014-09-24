// ====================================================================
// Test suite for implementation of fixed point iteration class
// ====================================================================

#include <gsl/gsl_vector.h>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_fixed_point_iterator

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "fixed_point_iterator.hpp"

using namespace flexiblesusy;

// Solution of \f$x^2 - 3x / 8 - 1/16 = 0\f$, solutions
// are 1/2 and -1/8, of which -1/8 is an attractive
// fixed point while 1/2 is unstable.
BOOST_AUTO_TEST_CASE( test_absolute_fixedpoint_1dim )
{
   struct FixedPoint1D {
      static int func(const gsl_vector* xold, 
                      void*, gsl_vector* xnew)
         {
            double x = gsl_vector_get(xold, 0);
            gsl_vector_set(xnew, 0, (16.*x*x - 1.) / 6.);
            return GSL_SUCCESS;
         }
   };

   const int max_iters = 100;
   const double tolerance = 1.0e-5;
   bool use_absolute_check = true;

   Fixed_point_iterator<1> fp_iter_1d(FixedPoint1D::func, NULL, max_iters, 
                                      tolerance, use_absolute_check);

   const double start[1] = { 0.499 };

   const int status = fp_iter_1d.find_fixed_point(start);

   const double soln = -0.125;

   BOOST_CHECK_EQUAL(status, GSL_SUCCESS);
   // N.B. CHECK_CLOSE_FRACTION(a,b,tol) == |a-b|/|b| < tol in this case,
   // which doesn't quite match the convergence criterion.
   BOOST_CHECK_LE(Abs(soln - fp_iter_1d.get_fixed_point(0)), tolerance);
}

BOOST_AUTO_TEST_CASE( test_relative_fixedpoint_1dim )
{
   struct FixedPoint1D {
      static int func(const gsl_vector* xold, 
                      void*, gsl_vector* xnew)
         {
            double x = gsl_vector_get(xold, 0);
            gsl_vector_set(xnew, 0, (16.*x*x - 1.) / 6.);
            return GSL_SUCCESS;
         }
   };

   const int max_iters = 100;
   const double tolerance = 1.0e-5;
   bool use_absolute_check = false;

   Fixed_point_iterator<1> fp_iter_1d(FixedPoint1D::func, NULL, max_iters, 
                                      tolerance, use_absolute_check);

   const double start[1] = { 0.499 };

   const int status = fp_iter_1d.find_fixed_point(start);

   const double soln = -0.125;

   BOOST_CHECK_EQUAL(status, GSL_SUCCESS);
   // N.B. in this case for x < y, the relative convergence criterion in the 
   // iteration is (|y| - |x|) / |y| \leq |y - x| / |y| < tolerance.
   BOOST_CHECK_CLOSE_FRACTION(fp_iter_1d.get_fixed_point(0), soln, tolerance);
}

BOOST_AUTO_TEST_CASE( test_fixedpoint_2dim )
{
   struct FixedPoint2D {
      static int func(const gsl_vector* xold, 
                      void*, gsl_vector* xnew)
         {
            double x = gsl_vector_get(xold, 0);
            double y = gsl_vector_get(xold, 1);

            return GSL_SUCCESS;
         }
   };

}

BOOST_AUTO_TEST_CASE( test_fixedpoint_3dim )
{

}

BOOST_AUTO_TEST_CASE( test_number_of_calls )
{

}

BOOST_AUTO_TEST_CASE( test_max_iters_termination )
{

}

BOOST_AUTO_TEST_CASE( test_runaway_iteration )
{

   struct Repeller {
      static int repel(const gsl_vector* xold, 
                       void*, gsl_vector* xnew)
         {
            gsl_vector_set(xnew, 0, Sqr(gsl_vector_get(xold, 0)) 
                           + 2.*gsl_vector_get(xold, 1));
            gsl_vector_set(xnew, 1, Sqr(gsl_vector_get(xold, 1)) 
                           - gsl_vector_get(xold, 0));
            return GSL_SUCCESS;
         }
   };
   
   Fixed_point_iterator<2> fp_iter(Repeller::repel, NULL, 100, 1.0e-5, false);

   const double start[2] = {1., 0.5};

   const int status = fp_iter.find_fixed_point(start);

   BOOST_CHECK_EQUAL(status, GSL_EBADFUNC);

}
