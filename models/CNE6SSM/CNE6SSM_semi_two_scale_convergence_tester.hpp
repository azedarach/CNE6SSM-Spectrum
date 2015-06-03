// ====================================================================
// Implementation of convergence tester for semianalytic version of
// the two-scale algorithm, checking the DRbar masses for convergence
// ====================================================================

#ifndef CNE6SSM_SEMI_TWO_SCALE_CONVERGENCE_TESTER_H
#define CNE6SSM_SEMI_TWO_SCALE_CONVERGENCE_TESTER_H

#include "CNE6SSM_semi_convergence_tester.hpp"
#include "CNE6SSM_semi_two_scale_model.hpp"
#include "two_scale_convergence_tester_drbar.hpp"

namespace flexiblesusy {

class Two_scale;

template<>
class CNE6SSM_semianalytic_convergence_tester<Two_scale> : public Convergence_tester_DRbar<CNE6SSM_semianalytic<Two_scale> > {
public:
   CNE6SSM_semianalytic_convergence_tester(CNE6SSM_semianalytic<Two_scale>*, double);
   virtual ~CNE6SSM_semianalytic_convergence_tester();

protected:
   virtual double max_rel_diff() const;
};

} // namespace flexiblesusy

#endif
