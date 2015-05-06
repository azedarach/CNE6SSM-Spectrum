// ====================================================================
// Implementation of convergence tester for CNE6SSM semianalytic solver
// ====================================================================

#ifndef CNE6SSM_SEMIANALYTIC_CONVERGENCE_TESTER_H
#define CNE6SSM_SEMIANALYTIC_CONVERGENCE_TESTER_H

#include "CNE6SSM_convergence_tester.hpp"
#include "CNE6SSM_semianalytic_model.hpp"
#include "semianalytic_convergence_tester_drbar.hpp"

namespace flexiblesusy {

class Semianalytic;

template<>
class CNE6SSM_convergence_tester<Semianalytic> : public Convergence_tester_DRbar<CNE6SSM<Semianalytic> > {
public:
   CNE6SSM_convergence_tester(CNE6SSM<Semianalytic>*, double);
   virtual ~CNE6SSM_convergence_tester();

protected:
   virtual double max_rel_diff_inner() const;
   virtual double max_rel_diff_outer() const;
};

} // namespace flexiblesusy

#endif
