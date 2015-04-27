// ====================================================================
// Implementation of convergence tester class to be used with 
// semianalytic solver
// ====================================================================

#include "semianalytic_convergence_tester.hpp"

#include "numerics.hpp"

namespace flexiblesusy {

Convergence_tester<Semianalytic>::~Convergence_tester()
{
}

bool Convergence_tester<Semianalytic>::is_equal(double a, double b) const
{
   return flexiblesusy::is_equal(a, b);
}

bool Convergence_tester<Semianalytic>::is_zero(double a) const
{
   return flexiblesusy::is_zero(a);
}

} // namespace flexiblesusy
