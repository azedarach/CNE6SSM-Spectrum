// ====================================================================
// Implementation of running precision calculator class to be used
// with semianalytic solver
// ====================================================================

#include "semianalytic_running_precision.hpp"
#include "logger.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace flexiblesusy {

Semianalytic_constant_precision::Semianalytic_constant_precision(double precision_)
   : Semianalytic_running_precision()
   , precision(precision_)
{
}

Semianalytic_constant_precision::~Semianalytic_constant_precision()
{
}

double Semianalytic_constant_precision::get_precision(unsigned)
{
   return precision;
}

Semianalytic_increasing_precision::Semianalytic_increasing_precision(double decreasing_factor_, double minimum_precision_)
   : Semianalytic_running_precision
} // namespace flexiblesusy
