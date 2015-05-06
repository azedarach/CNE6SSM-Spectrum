// ====================================================================
// Test implementation for a class to solve EWSB and do spectrum
// calculation using semianalytic algorithm
// ====================================================================

/**
 * @file CNE6SSM_semianalytic_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the semianalytic solver by solving EWSB
 *        and determine the pole masses and mixings
 */

#ifndef CNE6SSM_SEMIANALYTIC_MODEL
#define CNE6SSM_SEMIANALYTIC_MODEL

#include "CNE6SSM_model.hpp"
#include "CNE6SSM_mass_eigenstates.hpp"
#include "CNE6SSM_info.hpp"
#include "CNE6SSM_semianalytic_input_parameters.hpp"
#include "semianalytic_model.hpp"
#include "config.h"

#include <iosfwd>
#include <string>

#include <gsl/gsl_vector.h>
#include <Eigen/Core>

namespace flexiblesusy {

class EWSB_solver;
class Semianalytic;
/**
 * @class CNE6SSM<Semianalytic>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class CNE6SSM<Semianalytic> : public Semianalytic_model, public CNE6SSM_mass_eigenstates {
public:
   explicit CNE6SSM(const CNE6SSM_input_parameters<Semianalytic>& input_ = CNE6SSM_input_parameters<Semianalytic>());
   virtual ~CNE6SSM();

   virtual void clear();
   void set_ewsb_iteration_precision(double);
   void set_ewsb_loop_order(unsigned);
   void set_number_of_ewsb_iterations(std::size_t);
   const CNE6SSM_input_parameters<Semianalytic>& get_input() const;
   void set_input_parameters(const CNE6SSM_input_parameters<Semianalytic>&);
   double get_ewsb_iteration_precision() const;
   double get_ewsb_loop_order() const;
   int solve_ewsb_tree_level();
   int solve_ewsb_one_loop();
   int solve_ewsb();            ///< solve EWSB at ewsb_loop_order level

   // interface functions
   virtual void calculate_spectrum();
   virtual void clear_problems();
   virtual std::string name() const;
   virtual void run_to(double scale, double eps = -1.0);
   virtual void print(std::ostream&) const;
   virtual void set_precision(double);

   // model interface
   double get_parameter(unsigned) const;
   void set_parameter(unsigned, double);

private:
   struct EWSB_args {
      CNE6SSM<Semianalytic>* model;
      unsigned ewsb_loop_order;
   };

   // input parameters
   CNE6SSM_input_parameters<Semianalytic> input;

   std::size_t number_of_ewsb_iterations;
   unsigned ewsb_loop_order;
   double precision;              ///< RG running precision
   double ewsb_iteration_precision;
   static const std::size_t number_of_ewsb_equations = 5;

   int solve_ewsb_iteratively();
   int solve_ewsb_iteratively(unsigned);
   int solve_ewsb_iteratively_with(EWSB_solver*, const double[number_of_ewsb_equations]);
   int check_ewsb_solution(double);
   void ewsb_initial_guess(double[number_of_ewsb_equations]);
   int ewsb_step(double[number_of_ewsb_equations]);
   static int ewsb_step(const gsl_vector*, void*, gsl_vector*);
   void tadpole_equations(double[number_of_ewsb_equations]) const;
   static int tadpole_equations(const gsl_vector*, void*, gsl_vector*);

};

   std::ostream& operator<<(std::ostream&, const CNE6SSM<Semianalytic>&);

} // namespace flexiblesusy

#endif
