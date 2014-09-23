// ====================================================================
// Test suite for implementation of fixed point iteration class
// ====================================================================
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_fixed_point_iterator

#include <boost/test/unit_test.hpp>

#include "fixed_point_iterator.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_fixedpoint_1dim )
{

}

BOOST_AUTO_TEST_CASE( test_fixedpoint_2dim )
{

}

BOOST_AUTO_TEST_CASE( test_fixedpoint_3dim )
{

}

BOOST_AUTO_TEST_CASE( test_number_of_calls )
{

}

struct Repeller {
   static unsigned number_of_iters;
   static int repel(const std::vector<double>& xold, 
                    void*, std::vector<double>& xnew)
      {
         number_of_iters++;
         xnew[0] = Sqr(xold[0]) + 2.*xold[1];
         xnew[1] = Sqr(xold[1]) - xold[0];
         return GSL_SUCCESS;
      }
};

unsigned Repeller::number_of_iters = 0;

BOOST_AUTO_TEST_CASE( test_max_iters_termination )
{
   Fixed_point_iterator<2> fp_iter(Repeller::repel, NULL, 100, 1.0e-5, false);

   const std::vector<double> start = {1., 0.5};

   const int status = fp_iter.find_fixed_point(start);

   std::cout << "n = " << Repeller::number_of_iters << std::endl;
}

