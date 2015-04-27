// ====================================================================
// Convergence tester class to be used with semianalytic solver
// ====================================================================

#ifndef SEMIANALYTIC_CONVERGENCE_TESTER_H
#define SEMIANALYTIC_CONVERGENCE_TESTER_H

#include "convergence_tester.hpp"

namespace flexiblesusy {

class Semianalytic;

template<>
class Convergence_tester<Semianalytic> {
public:
   virtual ~Convergence_tester();
   virtual bool inner_accuracy_goal_reached() = 0;
   virtual bool outer_accuracy_goal_reached() = 0;
   virtual unsigned int max_iterations() const = 0;
   virtual void reset_inner_iteration_count() = 0;        ///< reset number of inner iterations to zero

protected:
   bool is_equal(double, double) const;
   bool is_zero(double) const;
};

}

#endif
