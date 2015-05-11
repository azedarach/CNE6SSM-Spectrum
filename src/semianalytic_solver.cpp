// ====================================================================
// Test implementation of semi-analytic solver routines
// ====================================================================

#include "semianalytic_solver.hpp"
#include "semianalytic_constraint.hpp"
#include "semianalytic_convergence_tester.hpp"
#include "semianalytic_initial_guesser.hpp"
#include "semianalytic_matching.hpp"
#include "semianalytic_model.hpp"
#include "two_scale_running_precision.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "functors.hpp"

#include <cmath>
#include <algorithm>
#include <iterator>
#include <limits>
#include <cassert>

/**
 * @file semianalytic_solver.cpp
 * @brief contains the implementation of the RGFlow<Semianalytic> class members
 */

namespace flexiblesusy {

RGFlow<Semianalytic>::RGFlow()
   : models()
   , inner_iteration(0)
   , outer_iteration(0)
   , convergence_tester(NULL)
   , initial_guesser(NULL)
   , running_precision_calculator(NULL)
   , inner_running_precision(1.0e-3)
   , outer_running_precision(1.0e-3)
   , model_at_this_scale(NULL)
   , only_solve_inner_iteration(false)
{
}

RGFlow<Semianalytic>::~RGFlow()
{
   delete_models();
}

/**
 * @brief Solves the boundary value problem.
 *
 * At first the initial_guess() is called.  Afterwards, the semianalytic
 * solution is applied: first, the function
 * iteratively runs the tower up and down and imposes the boundary
 * conditions.  The iteration stops if either the maximum number of
 * iterations is reached or the precision goal is achieved (defined by
 * the convergence_tester for the inner two-scale iteration). Then
 * the next, semianalytic step is applied. This is iterated until
 * either the maximum number of iterations is reached or the precision goal
 * is achieved (defined by the convergence_tester for the outer iteration).
 */
void RGFlow<Semianalytic>::solve()
{
   check_setup();

   const unsigned max_iterations = get_max_iterations();
   if (models.empty() || max_iterations == 0)
      return;

   initial_guess();

   outer_iteration = 0;
   bool outer_accuracy_reached = false;

   if (!only_solve_inner_iteration) {
      while (outer_iteration < max_iterations && !outer_accuracy_reached) {
         update_running_precision(Iteration_stage::Outer);

         // do inner two scale iteration
         solve_inner_iteration();
         
         // apply semianalytic step
         clear_problems();
         run_up(Iteration_stage::Outer);
         run_down(Iteration_stage::Outer);
         
         outer_accuracy_reached = accuracy_goal_reached(Iteration_stage::Outer);
         ++outer_iteration;
      }

      apply_lowest_constraint(Iteration_stage::Outer);
      
      if (!outer_accuracy_reached)
         throw NoConvergenceError(max_iterations);
   } else {
      solve_inner_iteration();
   }

   VERBOSE_MSG("convergence reached after " << iteration << " iterations");
}

void RGFlow<Semianalytic>::solve_inner_iteration()
{
   const unsigned max_iterations = get_max_iterations();

   reset_inner_iteration_convergence_test();
   inner_iteration = 0;
   bool inner_accuracy_reached = false;
   while (inner_iteration < max_iterations && !inner_accuracy_reached) {
      update_running_precision(Iteration_stage::Inner);
      clear_problems();
      run_up(Iteration_stage::Inner);
      run_down(Iteration_stage::Inner);
      inner_accuracy_reached = accuracy_goal_reached(Iteration_stage::Inner);
      ++inner_iteration;
   }

   apply_lowest_constraint(Iteration_stage::Inner);

   if (!inner_accuracy_reached)
      throw NoConvergenceError(max_iterations);

   VERBOSE_MSG("inner convergence reached after " << inner_iteration << " iterations");
}

/*
 * Sanity checks the models and boundary conditions
 */
void RGFlow<Semianalytic>::check_setup() const
{
   for (size_t m = 0; m < models.size(); ++m) {
      const TModel* model = models[m];
      if (!model->model) {
         std::stringstream message;
         message << "RGFlow<Semianalytic>::Error: model pointer ["
                 << m << "] is NULL";
         throw SetupError(message.str());
      }

      // check whether last model has a non-zero matching condition
      if (m + 1 == models.size()) {
         if (model->matching_condition)
            WARNING("the matching condition of the " << model->model->name()
                    << " is non-zero but will not be used");
      } else {
         if (model->matching_condition == NULL) {
            std::stringstream message;
            message << "RGFlow<Semianalytic>::Error: matching condition "
                    << "of the " << model->model->name() << " to the "
                    << models[m + 1]->model->name() << " is NULL";
            throw SetupError(message.str());
         }
      }
   }

   if (!convergence_tester) {
      throw SetupError("RGFlow<Semianalytic>::Error: convergence tester must "
                       "not be NULL");
   }
}

void RGFlow<Semianalytic>::clear_problems()
{
   VERBOSE_MSG("> clearing problems ...");

   const size_t number_of_models = models.size();
   for (size_t m = 0; m < number_of_models; ++m) {
      TModel* model = models[m];
      model->model->clear_problems();
   }
}

void RGFlow<Semianalytic>::delete_models()
{
   for_each(models.begin(), models.end(), Delete_object());
}

/**
 * Does the initial guess by calling the guess() method of the initial
 * guesser (if given).
 */
void RGFlow<Semianalytic>::initial_guess()
{
   if (initial_guesser)
      initial_guesser->guess();
}

/**
 * Run the model tower to the highest scale as part of the inner/outer iteration.
 * Thereby all upwards constraints are imposed (by calling the apply() functions) and the
 * matching conditions are applied (by calling
 * match_low_to_high_scale_model() functions).
 */
void RGFlow<Semianalytic>::run_up(Iteration_stage stage)
{
   VERBOSE_MSG("> running tower up (iteration " 
               << (stage == Iteration_stage::Inner ? 
                   inner_iteration : outer_iteration)  << ") ...");
   const size_t number_of_models = models.size();
   for (size_t m = 0; m < number_of_models; ++m) {
      TModel* model = models[m];
      model->model->set_precision(get_precision(stage));
      VERBOSE_MSG("> \tselecting model " << model->model->name());
      // apply all constraints
      const size_t n_upwards_constraints = (stage == Iteration_stage::Inner ?
                                            model->inner_upwards_constraints.size()
                                            : model->outer_upwards_constraints.size());
      for (size_t c = 0; c < n_upwards_constraints; ++c) {
         Constraint<Semianalytic>* constraint = 
            (stage == Iteration_stage::Inner ? model->inner_upwards_constraints[c]
             : model->outer_upwards_constraints[c]);
         const double scale = constraint->get_scale();
         VERBOSE_MSG("> \t\tselecting constraint " << c << " at scale " << scale);
         VERBOSE_MSG("> \t\t\trunning model to scale " << scale);
         model->model->run_to(scale);
         VERBOSE_MSG("> \t\t\tapplying constraint");
         constraint->apply();
      }
      // apply matching condition if this is not the last model
      if (m != number_of_models - 1) {
         VERBOSE_MSG("> \tmatching to model " << models[m + 1]->model->name());
         Matching<Semianalytic>* mc = model->matching_condition;
         const double scale = mc->get_scale();
         VERBOSE_MSG("> \t\t\trunning model to scale " << scale);
	 model->model->run_to(scale);
         VERBOSE_MSG("> \t\t\tapplying matching condition");
         mc->match_low_to_high_scale_model();
      }
   }
   VERBOSE_MSG("> running up finished");
}

/**
 * Run the model tower to the lowest scale as part of inner/outer iteration.
 * Thereby all downwards constraints are imposed (by calling the apply()
 * functions) and the matching conditions are applied (by calling
 * match_high_to_low_scale_model() functions).
 */
void RGFlow<Semianalytic>::run_down(Iteration_stage stage)
{
   assert(models.size() > 0 && "model size must not be zero");
   VERBOSE_MSG("< running tower down ...");
   const size_t number_of_models = models.size();
   for (long m = number_of_models - 1; m >= 0; --m) {
      TModel* model = models[m];
      VERBOSE_MSG("< \tselecting model " << model->model->name());
      // apply all constraints:
      // If m is the last model, do not apply the highest constraint,
      // because it was already applied when we ran up.
      const size_t c_begin = (m + 1 == (long)number_of_models ? 1 : 0);
      const size_t c_end = (stage == Iteration_stage::Inner ?
                            model->inner_downwards_constraints.size()
                            : model->outer_downwards_constraints.size());
      for (size_t c = c_begin; c < c_end; ++c) {
         Constraint<Semianalytic>* constraint = 
            (stage == Iteration_stage::Inner ? model->inner_downwards_constraints[c]
             : model->outer_downwards_constraints[c]);
         const double scale = constraint->get_scale();
         VERBOSE_MSG("< \t\tselecting constraint " << c << " at scale " << scale);
         VERBOSE_MSG("< \t\t\trunning model to scale " << scale);
         model->model->run_to(scale);
         // If m is the lowest energy model, do not apply the lowest
         // constraint, because it will be applied when we run up next
         // time.
         if (m != 0 || c + 1 != c_end) {
            VERBOSE_MSG("< \t\t\tapplying constraint");
            constraint->apply();
         }
      }
      // apply matching condition if this is not the first model
      if (m > 0) {
         Matching<Semianalytic>* mc = models[m - 1]->matching_condition;
         VERBOSE_MSG("< \tmatching to model " << models[m - 1]->model->name());
         const double scale = mc->get_scale();
         VERBOSE_MSG("< \t\t\trunning model to scale " << scale);
         model->model->run_to(scale);
         VERBOSE_MSG("> \t\t\tapplying matching condition");
         mc->match_high_to_low_scale_model();
      }
   }
   VERBOSE_MSG("< running down finished");
}

/**
 * Impose the constraint at the lowest energy scale (by calling the
 * apply() function).
 */
void RGFlow<Semianalytic>::apply_lowest_constraint(Iteration_stage stage)
{
   if (models.empty())
      return;

   TModel* model = models[0];
   model_at_this_scale = model->model;

   Constraint_container constraints;

   switch (stage) {
   case Iteration_stage::Inner: {
      constraints = model->inner_downwards_constraints;
      break;
   }
   case Iteration_stage::Outer: {
      constraints = model->outer_downwards_constraints;
      break;
   }
   }

   if (constraints.empty())
      return;

   Constraint<Semianalytic>* constraint = constraints.back();
   const double scale = constraint->get_scale();
   VERBOSE_MSG("| selecting constraint 0 at scale " << scale);
   VERBOSE_MSG("| \trunning model " << model->model->name() << " to scale " << scale);
   model->model->run_to(scale);
   VERBOSE_MSG("| \tapplying constraint");
   constraint->apply();
}

/**
 * Returns the precision of the RG running.
 *
 * @return RG running precision
 */
double RGFlow<Semianalytic>::get_precision(Iteration_stage stage)
{
   return (stage == Iteration_stage::Inner ? inner_running_precision
           : outer_running_precision);
}

/**
 * Recalculates the precision of the RG running using the user defined
 * Semianalytic_running_precision_calculator class.
 */
void RGFlow<Semianalytic>::update_running_precision(Iteration_stage stage)
{
   const unsigned int iteration = (stage == Iteration_stage::Inner ? inner_iteration
                                   : outer_iteration);

   if (running_precision_calculator) {
      const double precision = running_precision_calculator->get_precision(iteration);

      (stage == Iteration_stage::Inner ? inner_running_precision
       : outer_running_precision) = precision;
   }
}

void RGFlow<Semianalytic>::reset_inner_iteration_convergence_test()
{
   convergence_tester->reset_inner_iteration_count();
}

/**
 * Add a model and the corresponding model constraints.  With this
 * function the user can use different contraints for the up and down
 * running of the model tower.
 *
 * @param model model
 * @param inner_upwards_constraints model constraints for running up in inner iteration
 * @param inner_downwards_constraints model constraints for running down in inner iteration
 * @param outer_upwards_constraints model constraints for running up in outer iteration
 * @param outer_downwards_constraints model constraints for running down in outer iteration
 */
void RGFlow<Semianalytic>::add_model(Semianalytic_model* model,
                                     const std::vector<Constraint<Semianalytic>*>& inner_upwards_constraints,
                                     const std::vector<Constraint<Semianalytic>*>& inner_downwards_constraints,
                                     const std::vector<Constraint<Semianalytic>*>& outer_upwards_constraints,
                                     const std::vector<Constraint<Semianalytic>*>& outer_downwards_constraints)
{
   add_model(model, NULL, inner_upwards_constraints, inner_downwards_constraints,
             outer_upwards_constraints, outer_downwards_constraints);
}

/**
 * Add a model, the corresponding model constraints and the matching
 * condition to the next higher model.  With this function the user
 * can use different contraints for the up and down running of the
 * model tower.
 *
 * @param model model
 * @param mc matching condition to the next higher model
 * @param inner_upwards_constraints model constraints for running up in inner iteration
 * @param inner_downwards_constraints model constraints for running down in inner iteration
 * @param outer_upwards_constraints model constraints for running up in outer iteration
 * @param outer_downwards_constraints model constraints for running down in outer iteration
 */
void RGFlow<Semianalytic>::add_model(Semianalytic_model* model, 
                                     Matching<Semianalytic>* mc,
                                     const std::vector<Constraint<Semianalytic>*>& inner_upwards_constraints,
                                     const std::vector<Constraint<Semianalytic>*>& inner_downwards_constraints,
                                     const std::vector<Constraint<Semianalytic>*>& outer_upwards_constraints,
                                     const std::vector<Constraint<Semianalytic>*>& outer_downwards_constraints)
{
   TModel* new_model = new TModel(model, inner_upwards_constraints, inner_downwards_constraints,
                                  outer_upwards_constraints, outer_downwards_constraints, mc);

   for (Constraint_container::iterator it = new_model->inner_upwards_constraints.begin(),
           end = new_model->inner_upwards_constraints.end(); it != end; ++it)
      (*it)->set_model(model);

   for (Constraint_container::iterator it = new_model->inner_downwards_constraints.begin(),
           end = new_model->inner_downwards_constraints.end(); it != end; ++it)
      (*it)->set_model(model);

   for (Constraint_container::iterator it = new_model->outer_upwards_constraints.begin(),
           end = new_model->outer_upwards_constraints.end(); it != end; ++it)
      (*it)->set_model(model);

   for (Constraint_container::iterator it = new_model->outer_downwards_constraints.begin(),
           end = new_model->outer_downwards_constraints.end(); it != end; ++it)
      (*it)->set_model(model);

   models.push_back(new_model);
}

/**
 * Returns the value returned by the accuracy_goal_reached() method of
 * the convergence tester.
 */
bool RGFlow<Semianalytic>::accuracy_goal_reached(Iteration_stage stage) const
{
   return (stage == Iteration_stage::Inner ?
           convergence_tester->inner_accuracy_goal_reached()
           : convergence_tester->outer_accuracy_goal_reached());
}

/**
 * Sets a flag indicating whether to only solve the inner iteration
 *
 * @param flag flag determining if only inner iteration should be solved
 */
void RGFlow<Semianalytic>::solve_inner_iteration_only(bool flag)
{
   only_solve_inner_iteration = flag;
}

/**
 * Set the convergence tester to be used during the iteration.
 *
 * @param convergence_tester_ the convergence tester to be used
 */
void RGFlow<Semianalytic>::set_convergence_tester(Convergence_tester<Semianalytic>* convergence_tester_)
{
   convergence_tester = convergence_tester_;
}

/**
 * Set the initial guesser to be used to start the iteration.
 *
 * @param initial_guesser_ the initial guesser to be used
 */
void RGFlow<Semianalytic>::set_initial_guesser(Initial_guesser<Semianalytic>* initial_guesser_)
{
   initial_guesser = initial_guesser_;
}

/**
 * Set RG running precision calculator.
 *
 * @param rp running precision calculator
 */
void RGFlow<Semianalytic>::set_running_precision(Two_scale_running_precision* rp)
{
   running_precision_calculator = rp;
}

/**
 * Returns the number of performed iterations
 * @return number of performed iterations
 */
unsigned int RGFlow<Semianalytic>::number_of_iterations_done() const
{
   return outer_iteration;
}

/**
 * Returns the maximum number of iterations set in the convergence
 * tester.
 */
unsigned int RGFlow<Semianalytic>::get_max_iterations() const
{
   return convergence_tester->max_iterations();
}

/**
 * @brief resets the solver to the initial condition
 *
 * The pointers to the models, matching conditions, convergence
 * tester, initial guesser, and running precision calculator are set
 * to zero.  The running precision is set to the default value 0.001.
 */
void RGFlow<Semianalytic>::reset()
{
   delete_models();
   models.clear();

   outer_iteration = 0;
   inner_iteration = 0;
   convergence_tester = NULL;
   initial_guesser = NULL;
   running_precision_calculator = NULL;
   outer_running_precision = 1.0e-3;
   inner_running_precision = 1.0e-3;
   model_at_this_scale = NULL;
}

/**
 * Run the model tower to the given scale.
 *
 * @param scale scale to run to
 */
void RGFlow<Semianalytic>::run_to(double scale)
{
   // find model which is defined at `scale'
   model_at_this_scale = NULL;
   const size_t number_of_models = models.size();

   for (size_t m = 0; m < models.size(); ++m) {
      const TModel* model = models[m];
      double highest_scale, lowest_scale;

      if (!model) {
         std::ostringstream msg;
         msg << "RGFlow<Semianalytic>::run_to: pointer to model " << m
             << " is zero";
         throw SetupError(msg.str());
      }

      if (m != number_of_models - 1) {
         // if this is not the last model, the matching condition is
         // the highest scale
         const Matching<Semianalytic>* mc = model->matching_condition;
         if (!mc) {
            std::ostringstream msg;
            msg << "RGFlow<Semianalytic>::run_to: pointer to matching condition"
               " of model " << m << " is zero";
            throw SetupError(msg.str());
         }
         highest_scale = mc->get_scale();
      } else {
         // otherwise the last constraint is at the highest scale
         double inner_highest_scale, outer_highest_scale;

         if (model->inner_upwards_constraints.empty())
            inner_highest_scale = std::numeric_limits<double>::max();
         else 
            inner_highest_scale = model->inner_upwards_constraints.back()->get_scale();

         if (model->outer_upwards_constraints.empty())
            outer_highest_scale = std::numeric_limits<double>::max();
         else 
            outer_highest_scale = model->outer_upwards_constraints.back()->get_scale();

         highest_scale = (inner_highest_scale > outer_highest_scale ? inner_highest_scale
                          : outer_highest_scale);
      }

      if (m > 0) {
         // if this is not the first model, the previous matching
         // condition is the lowest scale
         lowest_scale = models[m-1]->matching_condition->get_scale();
      } else {
         // otherwise the first constraint is at the lowest scale
         double inner_lowest_scale, outer_lowest_scale;

         if (model->inner_upwards_constraints.empty())
            inner_lowest_scale = 0.;
         else
            inner_lowest_scale = model->inner_upwards_constraints[0]->get_scale();

         if (model->outer_upwards_constraints.empty())
            outer_lowest_scale = 0.;
         else
            outer_lowest_scale = model->outer_upwards_constraints[0]->get_scale();

         lowest_scale = (inner_lowest_scale < outer_lowest_scale ? inner_lowest_scale
                         : outer_lowest_scale);
      }

      if (lowest_scale <= scale && scale <= highest_scale) {
         model_at_this_scale = model->model;
         break;
      }
   }

   if (model_at_this_scale)
      model_at_this_scale->run_to(scale);
   else {
      std::ostringstream msg;
      msg << "No model at scale " << scale << " found!";
      throw SetupError(msg.str());
   }
}

/**
 * Returns the pointer to the model at the current scale.
 */
Semianalytic_model* RGFlow<Semianalytic>::get_model() const
{
   return model_at_this_scale;
}

} // namespace flexiblesusy
