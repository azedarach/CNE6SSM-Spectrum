// ====================================================================
// Class to do fixed point iteration.
// ====================================================================

#ifndef FIXED_POINT_ITERATOR_H
#define FIXED_POINT_ITERATOR_H

#include <iostream>
#include <cassert>
#include <vector>
#include <gsl/gsl_errno.h>

#include "logger.hpp"
#include "error.hpp"

namespace flexiblesusy {

/**
 * @class Fixed_point_iterator
 * @brief Does fixed point iteration
 *
 * Use std::vector instead of gsl_vector
 */
template <std::size_t dimension>
class Fixed_point_iterator {
public:

   Fixed_point_iterator();

   double get_fixed_point(std::size_t) const;
   int do_iteration(const std::vector<double> &);

private:
   std::size_t max_iterations;       ///< maximum number of iterations
   double precision;                 ///< precision goal
   std::vector<double> xn;           ///< current iteration point
   std::vector<double> fixed_point;  ///< vector of fixed point estimate
   void* parameters;                 ///< pointer to parameters

   void print_state(std::size_t) const;
};

/**
 * Default constructor
 */
template <std::size_t dimension>
Fixed_point_iterator<dimension>::Fixed_point_iterator()
   : max_iterations(100)
   , precision(1.0e-2)
   , parameters(NULL)
{}

/**
 * Start the iteration
 *
 * @param start starting point
 *
 * @return GSL error code (GSL_SUCCESS if fixed point found)
 */
template <std::size_t dimension>
int Fixed_point_iterator<dimension>::do_iteration(const std::vector<double> & start)
{

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
   for (std::vector<double>::const_iterator it = xn.begin(),
           end = xn.end(); it != end; ++it) {
      std::cout << " " << *it;
   }
   std::cout << ", x(n+1) =";
   for (std::vector<double>::const_iterator it = fixed_point.begin(),
           end = fixed_point.end(); it != end; ++it) {
      std::cout << " " << *it;
   }
   std::cout << "\n";
}

template <std::size_t dimension>
double Fixed_point_iterator<dimension>::get_fixed_point(std::size_t i) const
{
   assert(i < dimension && "Fixed_point_iterator<>::get_fixed_point: index out"
          " of bounds");
   return fixed_point.at(i);
}

} // namespace flexiblesusy

#endif
