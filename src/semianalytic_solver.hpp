// ====================================================================
// Test implementation of semi-analytic solver routines
// ====================================================================

#ifndef SEMIANALYTIC_SOLVER_H
#define SEMIANALYTIC_SOLVER_H

#include "rg_flow.hpp"

#include <vector>
#include <string>
#include <sstream>

/**
 * @file semianalytic_solver.hpp
 * @brief contains the definition of the RGFlow<Semianalytic> class
 */

namespace flexiblesusy {

template <class T> class Constraint;
template <class T> class Matching;
template <class T> class Convergence_tester;
template <class T> class Initial_guesser;
class Semianalytic;
class Semianalytic_model;
class Two_scale_running_precision;

/**
 * @class RGFlow<Semianalytic>
 * @brief Boundary condition solver (semianalytic algorithm)
 *
 * This boundary condition solver uses a semianalytic algorithm
 * to solve the boundary value problem. A two-scale iteration is
 * used to determine an initial set of parameters, followed by
 * a second step in which these parameters are used to to derive
 * the remaining parameters using a semi-analytic approach provided
 * by the user. The whole process is iterated until convergence is
 * reached.
 */

template<>
class RGFlow<Semianalytic> {
public:
   RGFlow();
   ~RGFlow();

   /// add models and constraints
   void add_model(Semianalytic_model*,
                  const std::vector<Constraint<Semianalytic>*>&,
                  const std::vector<Constraint<Semianalytic>*>&);

   /// add models and constraints
   void add_model(Semianalytic_model*,
                  const std::vector<Constraint<Semianalytic>*>&,
                  const std::vector<Constraint<Semianalytic>*>&,
                  const std::vector<Constraint<Semianalytic>*>&,
                  const std::vector<Constraint<Semianalytic>*>&);

   /// add model, constraints and matching condition
   void add_model(Semianalytic_model*,
                  Matching<Semianalytic>* m = NULL,
                  const std::vector<Constraint<Semianalytic>*>&
                  inner_constraints = std::vector<Constraint<Semianalytic>*>(),
                  const std::vector<Constraint<Semianalytic>*>&
                  outer_constraints = std::vector<Constraint<Semianalytic>*>());

   /// add models, constraints and matching condition
   void add_model(Semianalytic_model*,
                  Matching<Semianalytic>*,
                  const std::vector<Constraint<Semianalytic>*>&,
                  const std::vector<Constraint<Semianalytic>*>&,
                  const std::vector<Constraint<Semianalytic>*>&,
                  const std::vector<Constraint<Semianalytic>*>&);

   /// get model at current scale
   Semianalytic_model* get_model() const;
   /// get number of used iterations
   unsigned int number_of_iterations_done() const;
   /// clear all internal data
   void reset();
   /// pick valid model and run it to the given scale
   void run_to(double);
   /// set flag to only solve inner iteration
   void solve_inner_iteration_only(bool);
   /// set convergence tester
   void set_convergence_tester(Convergence_tester<Semianalytic>*);
   /// set running precision calculator
   void set_running_precision(Two_scale_running_precision*);
   /// set initial guesser
   void set_initial_guesser(Initial_guesser<Semianalytic>*);
   /// solves the boundary value problem
   void solve();

private:
   typedef std::vector<Constraint<Semianalytic>*> Constraint_container;

   /**
    * @class TModel
    * @brief contains model, constraints and matching condition
    *
    * This class lumps together the model, its constraints and the
    * matching condition to the next higher model.
    */
   struct TModel {
      Semianalytic_model* model;                          ///< the model
      Constraint_container inner_upwards_constraints; ///< upwards two-scale constraints
      Constraint_container inner_downwards_constraints; ///< downwards two-scale constraints
      Constraint_container outer_upwards_constraints; ///< upwards semianalytic constraints
      Constraint_container outer_downwards_constraints; ///< downwards semianalytic constraints
      Matching<Semianalytic>* matching_condition;         ///< matching condition

      TModel(Semianalytic_model* m,
             const Constraint_container& iuc,
             const Constraint_container& idc,
             const Constraint_container& ouc,
             const Constraint_container& odc,
             Matching<Semianalytic>* mc)
         : model(m)
         , inner_upwards_constraints(iuc)
         , inner_downwards_constraints(idc)
         , outer_upwards_constraints(ouc)
         , outer_downwards_constraints(odc)
         , matching_condition(mc)
         {}
   };  

   enum Iteration_stage : unsigned {Inner, Outer};

   std::vector<TModel*> models; ///< tower of models (from low to high scale)
   unsigned int inner_iteration; ///< inner iteration number (starting at 0)
   unsigned int outer_iteration; ///< outer iteration number (starting at 0)
   Convergence_tester<Semianalytic>* convergence_tester; ///< the convergence tester
   Initial_guesser<Semianalytic>* initial_guesser; ///< does initial guess
   Two_scale_running_precision* running_precision_calculator; ///< RG running precision calculator
   double inner_running_precision; ///< RG running precision used in inner solver
   double outer_running_precision; ///< RG running precision
   Semianalytic_model* model_at_this_scale; ///< model at current scale
   bool only_solve_inner_iteration; ///< flag to only solve inner two-scale iteration

   bool accuracy_goal_reached(Iteration_stage) const; ///< check if accuracy goal is reached
   void check_setup() const; ///< check the setup
   void clear_problems(); ///< clear model problems
   void delete_models(); ///< delete all models
   unsigned int get_max_iterations() const; ///< returns max. number of iterations
   void initial_guess(); ///< initial guess
   void reset_inner_iteration_convergence_test();
   void solve_inner_iteration(); ///< solves inner iteration
   void run_up(Iteration_stage); ///< run all models up during two-scale iteration
   void run_down(Iteration_stage); ///< run all models down during two-scale iteration
   void apply_lowest_constraint(Iteration_stage); ///< apply lowest constraint
   double get_precision(Iteration_stage); ///< returns running precision
   void update_running_precision(Iteration_stage); ///< update the RG running precision
};

} // namespace flexiblesusy

#endif
