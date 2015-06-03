// ====================================================================
// Test implementation of initial guesser to be used in semianalytic
// version of the two-scale algorithm
// ====================================================================

#ifndef CNE6SSM_SEMI_TWO_SCALE_INITIAL_GUESSER_H
#define CNE6SSM_SEMI_TWO_SCALE_INITIAL_GUESSER_H

#include "CNE6SSM_semi_initial_guesser.hpp"
#include "CNE6SSM_semi_two_scale_low_scale_constraint.hpp"
#include "CNE6SSM_semi_two_scale_susy_scale_constraint.hpp"
#include "CNE6SSM_semi_two_scale_high_scale_constraint.hpp"
#include "two_scale_initial_guesser.hpp"

#include <sstream>

namespace flexiblesusy {

template <class T>
class CNE6SSM_semianalytic;

class Two_scale;

/**
 * @class CNE6SSM_semianalytic_initial_guesser<Two_scale>
 * @brief initial guesser for the CNE6SSM
 */

template<>
class CNE6SSM_semianalytic_initial_guesser<Two_scale> : public Initial_guesser<Two_scale> {
public:
   CNE6SSM_semianalytic_initial_guesser(CNE6SSM_semianalytic<Two_scale>*,
                                        const QedQcd&,
                                        const CNE6SSM_semianalytic_low_scale_constraint<Two_scale>&,
                                        const CNE6SSM_semianalytic_susy_scale_constraint<Two_scale>&,
                                        const CNE6SSM_semianalytic_high_scale_constraint<Two_scale>&,
                                        CNE6SSM_semianalytic_high_scale_constraint<Two_scale>*);
   virtual ~CNE6SSM_semianalytic_initial_guesser();
   virtual void guess(); ///< initial guess

   void set_running_precision(double p) { running_precision = p; }

private:
   CNE6SSM_semianalytic<Two_scale>* model; ///< pointer to model class
   QedQcd oneset;   ///< Standard Model low-energy data
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
   CNE6SSM_semianalytic_low_scale_constraint<Two_scale> low_constraint;
   CNE6SSM_semianalytic_susy_scale_constraint<Two_scale> susy_constraint;
   CNE6SSM_semianalytic_high_scale_constraint<Two_scale> high_constraint;
   // Following is a pointer to the high-scale constraint for the
   // full iteration, to try to allow setting the initial guess to be
   // the first approximation for M_GUT with no thresholds
   CNE6SSM_semianalytic_high_scale_constraint<Two_scale>* main_iteration_high_constraint;

   static void guess_low_scale_parameters(CNE6SSM_semianalytic<Two_scale>*, const QedQcd&,
                                          double, double, double);
   void guess_high_scale_parameters();
   void guess_susy_parameters();
   void guess_soft_parameters();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();

   struct Initial_guess_low_scale_constraint : public Constraint<Two_scale> {
      CNE6SSM_semianalytic<Two_scale>* model;
      QedQcd oneset;
      double g1;
      double g2;
      double g3;
      double scale;

      Initial_guess_low_scale_constraint(CNE6SSM_semianalytic<Two_scale>* model_,
                                         const QedQcd& oneset_, double g1_, double g2_,
                                         double g3_, double scale_)
         : Constraint<Two_scale>()
         , model(model_)
         , oneset(oneset_)
         , g1(g1_)
         , g2(g2_)
         , g3(g3_)
         , scale(scale_)
         {}

      virtual void apply() { guess_low_scale_parameters(model, oneset, g1, g2, g3); }
      virtual double get_scale() const { return scale; }
      virtual void set_model(Two_scale_model* model_);
   };

};

} // namespace flexiblesusy

#endif
