// ====================================================================
// Specialisation of a constraint handler for setting up the solver
// when using the semianalytic algorithm.
// ====================================================================

#include "CNE6SSM_semianalytic_constraint_handler.hpp"

namespace flexiblesusy {

CNE6SSM_constraint_handler<Semianalytic>::CNE6SSM_constraint_handler()
   : high_scale_constraint()
   , susy_scale_constraint()
   , low_scale_constraint()
{
}

CNE6SSM_constraint_handler<Semianalytic>::~CNE6SSM_constraint_handler()
{
}

double CNE6SSM_constraint_handler<Semianalytic>::get_highest_scale() const
{
   return high_scale_constraint.get_scale();
}

double CNE6SSM_constraint_handler<Semianalytic>::get_susy_scale() const
{
   return susy_scale_constraint.get_scale();
}

double CNE6SSM_constraint_handler<Semianalytic>::get_lowest_scale() const
{
   return low_scale_constraint.get_scale();
}

void CNE6SSM_constraint_handler<Semianalytic>::initialize_constraints(Semianalytic_model* model, const QedQcd& oneset)
{
   high_scale_constraint.clear();
   susy_scale_constraint.clear();
   low_scale_constraint.clear();

   high_scale_constraint.set_model(model);
   susy_scale_constraint.set_model(model);
   low_scale_constraint.set_model(model);

   low_scale_constraint.set_sm_parameters(oneset);

   susy_scale_constraint.set_input_scale_constraint(&high_scale_constraint);

   high_scale_constraint.initialize();
   susy_scale_constraint.initialize();
   low_scale_constraint.initialize();
}

CNE6SSM_initial_guesser<Semianalytic> CNE6SSM_constraint_handler<Semianalytic>::get_initial_guesser(CNE6SSM<Semianalytic>* model, const QedQcd& oneset) const
{
   CNE6SSM_initial_guesser<Semianalytic> initial_guesser(model, oneset,
                                                         low_scale_constraint,
                                                         susy_scale_constraint,
                                                         high_scale_constraint);

   return initial_guesser;
}

void CNE6SSM_constraint_handler<Semianalytic>::add_constraints_to_solver(Semianalytic_model* model, RGFlow<Semianalytic>& solver)
{
   std::vector<Constraint<Semianalytic>*> inner_upward_constraints(2);
   inner_upward_constraints[0] = &low_scale_constraint;
   inner_upward_constraints[1] = &high_scale_constraint;

   std::vector<Constraint<Semianalytic>*> inner_downward_constraints(2);
   inner_downward_constraints[0] = &high_scale_constraint;
   inner_downward_constraints[1] = &low_scale_constraint;

   std::vector<Constraint<Semianalytic>*> outer_upward_constraints(1);
   outer_upward_constraints[0] = &susy_scale_constraint;

   std::vector<Constraint<Semianalytic>*> outer_downward_constraints(2);
   outer_downward_constraints[0] = &susy_scale_constraint;
   outer_downward_constraints[1] = &low_scale_constraint;

   solver.add_model(model, inner_upward_constraints, inner_downward_constraints,
                    outer_upward_constraints, outer_downward_constraints);
}

void CNE6SSM_constraint_handler<Semianalytic>::set_sm_parameters(const QedQcd& oneset)
{
   low_scale_constraint.set_sm_parameters(oneset);
}

const QedQcd& CNE6SSM_constraint_handler<Semianalytic>::get_sm_parameters() const
{
   return low_scale_constraint.get_sm_parameters();
}

} // namespace flexiblesusy
