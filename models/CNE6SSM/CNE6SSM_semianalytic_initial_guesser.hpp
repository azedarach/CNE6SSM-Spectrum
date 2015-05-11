// ====================================================================
// Test implementation of initial guesser for semianalytic
// algorithm
// ====================================================================

#ifndef CNE6SSM_SEMIANALYTIC_INITIAL_GUESSER_H
#define CNE6SSM_SEMIANALYTIC_INITIAL_GUESSER_H

#include "CNE6SSM_initial_guesser.hpp"
#include "CNE6SSM_semianalytic_low_scale_constraint.hpp"
#include "CNE6SSM_semianalytic_susy_scale_constraint.hpp"
#include "CNE6SSM_semianalytic_high_scale_constraint.hpp"
#include "semianalytic_initial_guesser.hpp"
#include "two_scale_initial_guesser.hpp"

#include <sstream>

namespace flexiblesusy {

template <class T>
class CNE6SSM;

class Two_scale;
class Semianalytic;

/**
 * @class CNE6SSM_initial_guesser<Semianalytic>
 * @brief initial guesser for CNE6SSM using semianalytic solver
 */

template<>
class CNE6SSM_initial_guesser<Semianalytic> : public Initial_guesser<Semianalytic> {
public:
   CNE6SSM_initial_guesser(CNE6SSM<Semianalytic>*,
                               const QedQcd&,
                               const CNE6SSM_low_scale_constraint<Semianalytic>&,
                               const CNE6SSM_susy_scale_constraint<Semianalytic>&,
                               const CNE6SSM_high_scale_constraint<Semianalytic>&);
   virtual ~CNE6SSM_initial_guesser();
   virtual void guess(); ///< initial guess

   void set_running_precision(double p) { running_precision = p; }
   void set_max_iterations(unsigned iterations) { max_iterations = iterations; }

private:
   struct Susy_parameters_guesser : Initial_guesser<Semianalytic> {
      CNE6SSM_initial_guesser<Semianalytic>* initial_guesser;
      void guess()
         {
            initial_guesser->guess_low_scale_parameters();
            initial_guesser->guess_high_scale_parameters();
         }
   };

   CNE6SSM<Semianalytic>* model; ///< pointer to model class
   QedQcd oneset;   ///< Standard Model low-energy date
   double mu_guess; ///< guessed DR-bar mass of up-quark
   double mc_guess; ///< guessed DR-bar mass of charm-quark
   double mt_guess; ///< guessed DR-bar mass of top-quark
   double md_guess; ///< guessed DR-bar mass of down-quark
   double ms_guess; ///< guessed DR-bar mass of strange-quark
   double mb_guess; ///< guessed DR-bar mass of bottom-quark
   double me_guess; ///< guessed DR-bar mass of electron
   double mm_guess; ///< guessed DR-bar mass of muon
   double mtau_guess; ///< guessed DR-bar mass of tau
   double running_precision; ///< Runge-Kutta RG running precision
   unsigned max_iterations; ///< max. iterations in tree level iteration
   CNE6SSM_low_scale_constraint<Semianalytic> low_constraint;
   CNE6SSM_susy_scale_constraint<Semianalytic> susy_constraint;
   CNE6SSM_high_scale_constraint<Semianalytic> high_constraint;

   void guess_low_scale_parameters();
   void guess_high_scale_parameters();
   void guess_susy_parameters();
   void guess_soft_parameters();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();

};

} // namespace flexiblesusy

#endif
