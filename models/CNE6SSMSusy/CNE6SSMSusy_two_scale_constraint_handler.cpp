// ====================================================================
// Specialisation of a constraint handler for setting up the solver
// when using the two_scale algorithm.
// ====================================================================

#include "CNE6SSMSusy_two_scale_constraint_handler.hpp"

namespace flexiblesusy {

CNE6SSMSusy_constraint_handler<Two_scale>::CNE6SSMSusy_constraint_handler()
   : high_scale_constraint()
   , low_scale_constraint()
{
}

CNE6SSMSusy_constraint_handler<Two_scale>::~CNE6SSMSusy_constraint_handler()
{
}


double CNE6SSMSusy_constraint_handler<Two_scale>::get_highest_scale() const
{
   return high_scale_constraint.get_scale();
}

double CNE6SSMSusy_constraint_handler<Two_scale>::get_lowest_scale() const
{
   return low_scale_constraint.get_scale();
}

void CNE6SSMSusy_constraint_handler<Two_scale>::initialize_constraints(Two_scale_model* model, const QedQcd& oneset)
{
   high_scale_constraint.clear();
   low_scale_constraint.clear();

   high_scale_constraint.set_model(model);
   low_scale_constraint.set_model(model);

   low_scale_constraint.set_sm_parameters(oneset);

   high_scale_constraint.initialize();
   low_scale_constraint.initialize();
}

CNE6SSMSusy_initial_guesser<Two_scale> CNE6SSMSusy_constraint_handler<Two_scale>::get_initial_guesser(CNE6SSMSusy<Two_scale>* model, const QedQcd& oneset) const
{
   CNE6SSMSusy_initial_guesser<Two_scale> initial_guesser(model, oneset,
                                                          low_scale_constraint,
                                                          high_scale_constraint);

   return initial_guesser;
}

void CNE6SSMSusy_constraint_handler<Two_scale>::add_constraints_to_solver(Two_scale_model* model, RGFlow<Two_scale>& solver)
{
   std::vector<Constraint<Two_scale>*> upward_constraints(2);
   upward_constraints[0] = &low_scale_constraint;
   upward_constraints[1] = &high_scale_constraint;

   std::vector<Constraint<Two_scale>*> downward_constraints(3);
   downward_constraints[0] = &high_scale_constraint;
   downward_constraints[2] = &low_scale_constraint;

   solver.add_model(model, upward_constraints, downward_constraints);
}

void CNE6SSMSusy_constraint_handler<Two_scale>::set_sm_parameters(const QedQcd& oneset)
{
   low_scale_constraint.set_sm_parameters(oneset);
}

const QedQcd& CNE6SSMSusy_constraint_handler<Two_scale>::get_sm_parameters() const
{
   return low_scale_constraint.get_sm_parameters();
}

} // namespace flexiblesusy
