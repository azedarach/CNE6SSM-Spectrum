// ====================================================================
// Class to do fixed point iteration. Uses std::vector instead of
// gsl_vector at the moment, as the syntax is slightly nicer.
//
// TODO:
//   - implement check for no progress towards solution
// ====================================================================

#ifndef FIXED_POINT_ITERATOR_H
#define FIXED_POINT_ITERATOR_H

#include <iostream>
#include <cassert>
#include <vector>
#include <gsl/gsl_errno.h>

#include "wrappers.hpp"
#include "logger.hpp"
#include "error.hpp"

namespace flexiblesusy {

/**
 * @class Fixed_point_iterator
 * @brief Does fixed point iteration
 *
 * The user has to provide the function (of which a fixed point 
 * should be found) of the type Function_t. This function gets as
 * arguments a std::vector<double> of length 'dimension', a 
 * pointer to the parameters (of type void*) and a std::vector<double>
 * where the next point must be stored.
 *
 * Example:
 * @code
 * struct MyExample {
 *    static int func(const std::vector<double>& xold, void*, 
 *                    std::vector<double>& xnew)
 *       {
 *         return GSL_SUCCESS;
 *       }
 * };
 * 
 * Fixed_point_iterator<2> fp_iter(MyExample::func, NULL, 100, 1.0e-5);
 * const std::vector<double> start = {a, b};
 * const int status = fp_iter.find_fixed_point(start);
 * @endcode
 */
template <std::size_t dimension>
class Fixed_point_iterator {
public:
   typedef int (*Function_t)(const std::vector<double> &, void*,
                             std::vector<double>&);

   Fixed_point_iterator();
   Fixed_point_iterator(Function_t, void*, std::size_t, double, bool);

   double get_fixed_point(std::size_t) const;
   void set_function(Function_t f) { function = f; }
   void set_parameters(void* m) { parameters = m; }
   void set_precision(double p) { precision = p; }
   void set_max_iterations(std::size_t n) { max_iterations = n; }
   void set_test_absolute_errors(bool ae) { test_on_absolute = ae; } 
   int find_fixed_point(const std::vector<double> &);

private:
   std::size_t max_iterations;       ///< maximum number of iterations
   double precision;                 ///< precision goal
   bool test_on_absolute;            ///< use absolute convergence criterion
   std::vector<double> xn;           ///< current iteration point
   std::vector<double> fixed_point;  ///< vector of fixed point estimate
   void* parameters;                 ///< pointer to parameters
   Function_t function;              ///< function defining fixed point

   int fixed_point_iterator_iterate();
   int fixed_point_iterator_test_relative() const;
   int fixed_point_iterator_test_absolute() const;
   void print_state(std::size_t) const;
};

/**
 * Default constructor
 */
template <std::size_t dimension>
Fixed_point_iterator<dimension>::Fixed_point_iterator()
   : max_iterations(100)
   , precision(1.0e-2)
   , test_on_absolute(false)
   , xn(dimension)
   , fixed_point(dimension)
   , parameters(NULL)
   , function(NULL)
{}

/**
 * Constructor
 *
 * @param function_ pointer to the function to find fixed point for
 * @param parameters_ pointer to the parameters (for example the model)
 * @param max_iterations_ maximum number of iterations
 * @param precision_ precision goal
 * @param absolute_ use absolute convergence test
 */
template <std::size_t dimension>
Fixed_point_iterator<dimension>::Fixed_point_iterator(Function_t function_, 
                                                      void* parameters_, 
                                                      std::size_t max_iterations_, 
                                                      double precision_, bool absolute_)
   : max_iterations(max_iterations_)
   , precision(precision_)
   , test_on_absolute(absolute_)
   , xn(dimension)
   , fixed_point(dimension)
   , parameters(parameters_)
   , function(function_)
{}

/**
 * Start the iteration
 *
 * @param start starting point
 *
 * @return GSL error code (GSL_SUCCESS if fixed point found)
 */
template <std::size_t dimension>
int Fixed_point_iterator<dimension>::find_fixed_point(const std::vector<double> & start)
{
   assert(function && "Fixed_point_iterator<dimension>::find_fixed_point: function pointer"
          " must not be zero!");
   assert(start.size() == dimension && "Fixed_point_iterator<dimension>::find_fixed_point: start vector wrong size!");

   int status;
   std::size_t iter = 0;

#ifndef ENABLE_DEBUG
   gsl_set_error_handler_off();
#endif

   fixed_point = start;

#ifdef ENABLE_VERBOSE
   print_state(iter);
#endif

   do {
      iter++;
      status = fixed_point_iterator_iterate();

#ifdef ENABLE_VERBOSE
      print_state(iter);
#endif

      if (status)   // check if iterator has problems
         break;

      status = (test_on_absolute ? 
                fixed_point_iterator_test_absolute() 
                : fixed_point_iterator_test_relative());

   } while (status == GSL_CONTINUE && iter < max_iterations);

#ifdef ENABLE_VERBOSE
   std::cout << "\tFixed_point_iterator status = " 
             << gsl_strerror(status) << "\n";
#endif

   return status;
}

/**
 * Perform a single step of the fixed point iteration
 *
 * @return GSL error code
 */
template <std::size_t dimension>
int Fixed_point_iterator<dimension>::fixed_point_iterator_iterate()
{
   xn = fixed_point;
   int status = function(xn, parameters, fixed_point);

   if (status != GSL_SUCCESS) {
      return GSL_EBADFUNC;
   }

   return GSL_SUCCESS;
}

/**
 * Test whether the relative difference is less
 * than the set precision. The relative difference test used
 * here is identical to that defined in wrappers.hpp 
 * (see MaxRelDiff defined there), applied to each element of 
 * the vector.
 *
 * @return GSL error code (GSL_SUCCESS or GSL_CONTINUE)
 */
template <std::size_t dimension>
int Fixed_point_iterator<dimension>::fixed_point_iterator_test_relative() const
{
   double rel_diff = 0.;

   if (precision < 0.) {
      GSL_ERROR("relative tolerance is negative", GSL_EBADTOL);
   }

   for (std::size_t i =0; i < dimension; ++i) {
      rel_diff = MaxRelDiff(xn[i], fixed_point[i]);

      if (rel_diff > precision) {
         return GSL_CONTINUE;
      }
   }

   return GSL_SUCCESS;
}

/**
 * Test whether the absolute value of the residual, defined by
 * |x_{n+1}-x_n| = \f$\sqrt{\sum_i (x_{n+1}(i) - x_n(i))^2}\f$, 
 * is less than the set precision.
 *
 * @return GSL error code (GSL_SUCCESS or GSL_CONTINUE)
 */
template <std::size_t dimension>
int Fixed_point_iterator<dimension>::fixed_point_iterator_test_absolute() const
{
   double residual = 0;

   if (precision < 0.) {
      GSL_ERROR("absolute tolerance is negative", GSL_EBADTOL);
   }

   for (std::size_t i = 0; i < dimension; ++i) {
      residual += Sqr(fixed_point[i] - xn[i]);
   }

   residual = Sqrt(residual);

   return (residual < precision ? GSL_SUCCESS : GSL_CONTINUE);
}

/**
 * Print state of the fixed point iterator
 *
 * @param iteration iteration number
 */
template <std::size_t dimension>
void Fixed_point_iterator<dimension>::print_state(std::size_t iteration) const
{
   std::cout << "\tIteration n = " << iteration << ": x(n) =";
   for (std::size_t i = 0; i < dimension; ++i) {
      std::cout << " " << xn[i];
   }
   std::cout << ", x(n+1) =";
   for (std::size_t i = 0; i < dimension; ++i) {
      std::cout << " " << fixed_point[i];
   }
   std::cout << "\n";
}

template <std::size_t dimension>
double Fixed_point_iterator<dimension>::get_fixed_point(std::size_t i) const
{
   assert(i < dimension && "Fixed_point_iterator<>::get_fixed_point: index out"
          " of bounds");
   return fixed_point[i];
}

} // namespace flexiblesusy

#endif
