// ====================================================================
// Specialisation of a constraint handler for setting up the solver
// when using the two_scale algorithm.
// ====================================================================

#ifndef CNE6SSM_TWO_SCALE_CONSTRAINT_HANDLER_H
#define CNE6SSM_TWO_SCALE_CONSTRAINT_HANDLER_H

#include "CNE6SSM_constraint_handler.hpp" 
#include "CNE6SSM_two_scale_high_scale_constraint.hpp"
#include "CNE6SSM_two_scale_susy_scale_constraint.hpp"
#include "CNE6SSM_two_scale_low_scale_constraint.hpp"
#include "CNE6SSM_two_scale_initial_guesser.hpp"

#include "two_scale_solver.hpp"

namespace flexiblesusy {

class Two_scale;
class Two_scale_model;

template<>
class CNE6SSM_constraint_handler<Two_scale> {
public:
   CNE6SSM_constraint_handler();
   virtual ~CNE6SSM_constraint_handler();

   double get_highest_scale() const;
   double get_susy_scale() const;
   double get_lowest_scale() const;

   void initialize_constraints(Two_scale_model*, const QedQcd&);
   CNE6SSM_initial_guesser<Two_scale> get_initial_guesser(CNE6SSM<Two_scale>*,
                                                          const QedQcd&) const;
   void add_constraints_to_solver(Two_scale_model*, RGFlow<Two_scale>&);

   void set_sm_parameters(const QedQcd& oneset);
   const QedQcd& get_sm_parameters() const;

private:
   CNE6SSM_high_scale_constraint<Two_scale> high_scale_constraint;
   CNE6SSM_susy_scale_constraint<Two_scale> susy_scale_constraint;
   CNE6SSM_low_scale_constraint<Two_scale> low_scale_constraint;
};

} // namespace flexiblesusy

#endif
