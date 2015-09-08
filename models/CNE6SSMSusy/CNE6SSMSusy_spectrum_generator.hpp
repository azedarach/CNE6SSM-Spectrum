// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Sun 19 Apr 2015 20:37:13

#ifndef CNE6SSMSusy_SPECTRUM_GENERATOR_H
#define CNE6SSMSusy_SPECTRUM_GENERATOR_H

#include "CNE6SSMSusy_two_scale_model.hpp"
#include "CNE6SSMSusy_two_scale_convergence_tester.hpp"
#include "CNE6SSMSusy_two_scale_initial_guesser.hpp"
#include "CNE6SSMSusy_two_scale_constraint_handler.hpp"

#include "CNE6SSMSusy_utilities.hpp"

#include "coupling_monitor.hpp"
#include "error.hpp"
#include "two_loop_corrections.hpp"
#include "numerics2.hpp"
#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"

namespace flexiblesusy {

template <class T>
class CNE6SSMSusy_spectrum_generator {
public:
   CNE6SSMSusy_spectrum_generator()
      : solver(), model()
      , constraint_handler()
      , high_scale(0.)
      , low_scale(0.)
      , parameter_output_scale(0.)
      , precision_goal(1.0e-4)
      , max_iterations(0)
      , beta_loop_order(2)
   {}
   ~CNE6SSMSusy_spectrum_generator() {}

   double get_high_scale() const { return high_scale; }
   double get_low_scale()  const { return low_scale;  }
   const CNE6SSMSusy<T>& get_model() const { return model; }
   const Problems<CNE6SSMSusy_info::NUMBER_OF_PARTICLES>& get_problems() const {
      return model.get_problems();
   }
   int get_exit_code() const { return get_problems().have_problem(); }
   void set_parameter_output_scale(double s) { parameter_output_scale = s; }
   void set_precision_goal(double precision_goal_) { precision_goal = precision_goal_; }
   void set_beta_loop_order(unsigned l) { beta_loop_order = l; }
   void set_max_iterations(unsigned n) { max_iterations = n; }

   void run(const QedQcd& oneset, const CNE6SSMSusy_input_parameters<T>& input);
   void write_running_couplings(const std::string& filename = "CNE6SSMSusy_rge_running.dat") const;

private:
   RGFlow<T> solver;
   CNE6SSMSusy<T> model;
   CNE6SSMSusy_constraint_handler<T> constraint_handler;
   double high_scale, low_scale;
   double parameter_output_scale; ///< output scale for running parameters
   double precision_goal; ///< precision goal
   unsigned max_iterations; ///< maximum number of iterations
   unsigned beta_loop_order; ///< beta-function loop order
};

/**
 * @brief Run's the RG solver with the given input parameters
 *
 * This function sets up the RG solver using a high-scale, susy-scale
 * and low-scale constraint.  Afterwards the solver is run until
 * convergence is reached or an error occours.  Finally the particle
 * spectrum (pole masses) is calculated.
 *
 * @param oneset Standard Model input parameters
 * @param input model input parameters
 */
template <class T>
void CNE6SSMSusy_spectrum_generator<T>::run(const QedQcd& oneset,
                                            const CNE6SSMSusy_input_parameters<T>& input)
{
   model.clear();
   model.set_input_parameters(input);
   model.set_loops(beta_loop_order);

   constraint_handler.initialize_constraints(&model, oneset);

   CNE6SSMSusy_convergence_tester<T> convergence_tester(&model, precision_goal);
   if (max_iterations > 0)
      convergence_tester.set_max_iterations(max_iterations);

   CNE6SSMSusy_initial_guesser<T> initial_guesser(constraint_handler.get_initial_guesser(&model, oneset));

   Two_scale_increasing_precision precision(10.0, precision_goal);

   solver.reset();
   solver.set_convergence_tester(&convergence_tester);
   solver.set_running_precision(&precision);
   solver.set_initial_guesser(&initial_guesser);
   constraint_handler.add_constraints_to_solver(&model, solver);

   high_scale = low_scale = 0.;

   try {
      solver.solve();
      high_scale = constraint_handler.get_highest_scale();
      low_scale = constraint_handler.get_lowest_scale();

      // run to output scale (if scale > 0)
      if (!is_zero(parameter_output_scale)) {
         model.run_to(parameter_output_scale);
      }
   } catch (const NoConvergenceError&) {
      model.get_problems().flag_no_convergence();
   } catch (const NonPerturbativeRunningError&) {
      model.get_problems().flag_no_perturbative();
   } catch (const NoRhoConvergenceError&) {
      model.get_problems().flag_no_rho_convergence();
   } catch (const Error& error) {
      model.get_problems().flag_thrown(error.what());
   } catch (const std::string& str) {
      model.get_problems().flag_thrown(str);
   } catch (const char* str) {
      model.get_problems().flag_thrown(str);
   } catch (const std::exception& error) {
      model.get_problems().flag_thrown(error.what());
   }
}

/**
 * Create a text file which contains the values of all model
 * parameters at all scales between the low-scale and the high-scale.
 *
 * @param filename name of output file
 */
template <class T>
void CNE6SSMSusy_spectrum_generator<T>::write_running_couplings(const std::string& filename) const
{
   CNE6SSMSusy<T> tmp_model(model);
   tmp_model.run_to(low_scale);

   CNE6SSMSusy_parameter_getter parameter_getter;
   Coupling_monitor<CNE6SSMSusy<T>, CNE6SSMSusy_parameter_getter>
      coupling_monitor(tmp_model, parameter_getter);

   coupling_monitor.run(low_scale, high_scale, 100, true);
   coupling_monitor.write_to_file(filename);
}

} // namespace flexiblesusy

#endif
