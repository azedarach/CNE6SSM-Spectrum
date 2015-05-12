// ====================================================================
// Specialisation of a constraint handler for setting up the solver
// when using the semianalytic algorithm.
// ====================================================================

#ifndef CNE6SSM_SEMIANALYTIC_CONSTRAINT_HANDLER_H
#define CNE6SSM_SEMIANALYTIC_CONSTRAINT_HANDLER_H

#include "CNE6SSM_constraint_handler.hpp"
#include "CNE6SSM_semianalytic_high_scale_constraint.hpp"
#include "CNE6SSM_semianalytic_susy_scale_constraint.hpp"
#include "CNE6SSM_semianalytic_low_scale_constraint.hpp"
#include "CNE6SSM_semianalytic_initial_guesser.hpp"

#include "semianalytic_solver.hpp"

namespace flexiblesusy {

class Semianalytic;
class Semianalytic_model;

template<>
class CNE6SSM_constraint_handler<Semianalytic> {
public:
   CNE6SSM_constraint_handler();
   virtual ~CNE6SSM_constraint_handler();

   double get_highest_scale() const;
   double get_susy_scale() const;
   double get_lowest_scale() const;

   void initialize_constraints(Semianalytic_model*, const QedQcd&);
   CNE6SSM_initial_guesser<Semianalytic> get_initial_guesser(CNE6SSM<Semianalytic>*,
                                                             const QedQcd&) const;
   void add_constraints_to_solver(Semianalytic_model*, RGFlow<Semianalytic>&);
   void set_sm_parameters(const QedQcd& oneset);
   const QedQcd& get_sm_parameters() const;

private:
   CNE6SSM_high_scale_constraint<Semianalytic> high_scale_constraint;
   CNE6SSM_susy_scale_constraint<Semianalytic> susy_scale_constraint;
   CNE6SSM_low_scale_constraint<Semianalytic> low_scale_constraint;
};

} // namespace flexiblesusy

#endif
