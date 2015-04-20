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

   const std::size_t max_iters = 100;
   const double tolerance = 1.0e-5;

   fixed_point_iterator::Convergence_tester_absolute convergence_tester(tolerance);

   Fixed_point_iterator<1, fixed_point_iterator::Convergence_tester_absolute> 
      fp_iter_1d(FixedPoint1D::func, NULL, max_iters, convergence_tester);

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

   const std::size_t max_iters = 100;
   const double tolerance = 1.0e-5;

   fixed_point_iterator::Convergence_tester_relative convergence_tester(tolerance);

   Fixed_point_iterator<1, fixed_point_iterator::Convergence_tester_relative> 
      fp_iter_1d(FixedPoint1D::func, NULL, max_iters, convergence_tester);

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

   const std::size_t max_iters = 100;
   const double tolerance = 1.0e-5;

   fixed_point_iterator::Convergence_tester_relative convergence_tester(tolerance);   

   Fixed_point_iterator<2> fp_iter(Repeller::repel, NULL, max_iters, convergence_tester);

   const double start[2] = {1., 0.5};

   const int status = fp_iter.find_fixed_point(start);

   BOOST_CHECK_EQUAL(status, GSL_EBADFUNC);

}

class Parabola {
public:
   static void reset() { number_of_calls = 0; }
   static unsigned get_number_of_calls() { return number_of_calls; }

   /**
    * Finding root of f(x,y) = ((x-5)^2, (y-1)^2),
    *
    * => Update steps
    *
    * (x,y) = (-25/(x-10), -1/(y-2))
    *
    * @param x touple (x,y)
    *
    * @return fixed point iteration update steps
    */
   static int func(const gsl_vector* x, void*, gsl_vector* f)
      {
         const double y = gsl_vector_get(x, 0);
         const double z = gsl_vector_get(x, 1);
         //gsl_vector_set(f, 0, -25./(y - 10.));
         //gsl_vector_set(f, 1, -1./(z - 2.));
         gsl_vector_set(f, 0, 0.1 * (25.0 + y * y));
         gsl_vector_set(f, 1, 0.5 * (1.0 + z * z));
         number_of_calls++;
         return GSL_SUCCESS;
      }

private:
   static unsigned number_of_calls;
};

unsigned Parabola::number_of_calls = 0;

BOOST_AUTO_TEST_CASE( test_parabola_2dim )
{
   const double precision = 1.0e-4;
   const std::size_t max_iters = 1000;
   const double start[2] = {4.95, 1.01};
   fixed_point_iterator::Convergence_tester_relative convergence_tester(precision);
   Fixed_point_iterator<2> fpi(Parabola::func, NULL, max_iters, convergence_tester);
   int status = GSL_SUCCESS;

   Parabola::reset();

   status = fpi.find_fixed_point(start);

   const double residual_1 = MaxRelDiff(5.0, fpi.get_fixed_point(0));
   const double residual_2 = MaxRelDiff(1.0, fpi.get_fixed_point(1));

   // Note: The convergence criterion 
   // MaxRelDiff(x_{n+1}, x_{n}) < precision
   // is not very good: the method converges slowly. This means
   // subsequent steps are very close to each other, but x_n
   // might not be close to the true fixed point.

   BOOST_REQUIRE(status == GSL_SUCCESS);
   BOOST_CHECK_LT(residual_1, 100*precision);
   BOOST_CHECK_LT(residual_2, 100*precision);
   BOOST_MESSAGE("fixed point iterator used " << Parabola::get_number_of_calls() << " calls");
}
