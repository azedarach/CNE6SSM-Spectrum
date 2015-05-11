// ====================================================================
// Test implementation of initial guesser for semianalytic
// algorithm
// ====================================================================

#include "CNE6SSM_semianalytic_initial_guesser.hpp"
#include "CNE6SSM_semianalytic_model.hpp"
#include "CNE6SSM_semianalytic_convergence_tester.hpp"
#include "lowe.h"
#include "error.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"
#include "two_scale_running_precision.hpp"
#include "semianalytic_solver.hpp"
#include "config.h"

#include <Eigen/Core>
#include <cassert>
#include <vector>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define MODEL model

CNE6SSM_initial_guesser<Semianalytic>::CNE6SSM_initial_guesser(
   CNE6SSM<Semianalytic>* model_,
   const QedQcd& oneset_,
   const CNE6SSM_low_scale_constraint<Semianalytic>& low_constraint_,
   const CNE6SSM_susy_scale_constraint<Semianalytic>& susy_constraint_,
   const CNE6SSM_high_scale_constraint<Semianalytic>& high_constraint_
)
   : Initial_guesser<Semianalytic>()
   , model(model_)
   , oneset(oneset_)
   , mu_guess(0.)
   , mc_guess(0.)
   , mt_guess(0.)
   , md_guess(0.)
   , ms_guess(0.)
   , mb_guess(0.)
   , me_guess(0.)
   , mm_guess(0.)
   , mtau_guess(0.)
   , running_precision(1.0e-3)
   , max_iterations(0)
   , low_constraint(low_constraint_)
   , susy_constraint(susy_constraint_)
   , high_constraint(high_constraint_)
{
   assert(model && "CNE6SSM_initial_guesser: Error: pointer to model"
          " CNE6SSM<Semianalytic> must not be zero");
}

CNE6SSM_initial_guesser<Semianalytic>::~CNE6SSM_initial_guesser()
{
}

/**
 * Guesses the DR-bar model parameters by calling
 * guess_susy_parameters() and guess_soft_parameters().
 */
void CNE6SSM_initial_guesser<Semianalytic>::guess()
{
   guess_susy_parameters();
   guess_soft_parameters();
}

/**
 * Guesses the SUSY parameters (gauge, Yukawa couplings) at
 * \f$m_\text{top}^\text{pole}\f$ from the Standard Model gauge
 * couplings and fermion masses. Threshold corrections are ignored.
 * The user-defined initial guess at the low-scale
 * (InitialGuessAtLowScale) is applied here:
 *
 * \code{.cpp}
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto ssumInput = INPUTPARAMETER(ssumInput);

   MODEL->set_vd(LowEnergyConstant(vev)/Sqrt(1 + Sqr(TanBeta)));
   MODEL->set_vu((TanBeta*LowEnergyConstant(vev))/Sqrt(1 + Sqr(TanBeta)));
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
   MODEL->set_XiF(Power(LowEnergyConstant(MZ),3));
   MODEL->set_LXiF(Power(LowEnergyConstant(MZ),3));
   MODEL->set_vs(0.7071067811865475*sInput);
   MODEL->set_vsb(0.7071067811865475*sInput);
   MODEL->set_vphi(0.5*sInput);
   MODEL->set_QS(QSInput);

 * \endcode
 */
void CNE6SSM_initial_guesser<Semianalytic>::guess_low_scale_parameters()
{
   QedQcd leAtMt(oneset);
   const double MZ = Electroweak_constants::MZ;
   const double MW = Electroweak_constants::MW;
   const double sinThetaW2 = 1.0 - Sqr(MW / MZ);
   const double mtpole = leAtMt.displayPoleMt();

   mu_guess = leAtMt.displayMass(mUp);
   mc_guess = leAtMt.displayMass(mCharm);
   mt_guess = leAtMt.displayMass(mTop) - 30.0;
   md_guess = leAtMt.displayMass(mDown);
   ms_guess = leAtMt.displayMass(mStrange);
   mb_guess = leAtMt.displayMass(mBottom);
   me_guess = leAtMt.displayMass(mElectron);
   mm_guess = leAtMt.displayMass(mMuon);
   mtau_guess = leAtMt.displayMass(mTau);

   // guess gauge couplings at mt
   const DoubleVector alpha_sm(leAtMt.getGaugeMu(mtpole, sinThetaW2));

   model->set_g1(sqrt(4.0 * M_PI * alpha_sm(1)));
   model->set_g2(sqrt(4.0 * M_PI * alpha_sm(2)));
   model->set_g3(sqrt(4.0 * M_PI * alpha_sm(3)));
   model->set_scale(mtpole);

   // apply user-defined initial guess at the low scale
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto sInput = INPUTPARAMETER(sInput);
   const auto QSInput = INPUTPARAMETER(QSInput);

   MODEL->set_vd(LowEnergyConstant(vev)/Sqrt(1 + Sqr(TanBeta)));
   MODEL->set_vu((TanBeta*LowEnergyConstant(vev))/Sqrt(1 + Sqr(TanBeta)));
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
   MODEL->set_XiF(Power(LowEnergyConstant(MZ),3));
   MODEL->set_LXiF(Power(LowEnergyConstant(MZ),3));
   MODEL->set_vs(0.7071067811865475 * sInput);
   MODEL->set_vsb(0.7071067811865475 * sInput);
   MODEL->set_vphi(0.5 * sInput);
   MODEL->set_QS(QSInput);

}

void CNE6SSM_initial_guesser<Semianalytic>::calculate_DRbar_yukawa_couplings()
{
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

/**
 * Calculates the Yukawa couplings Yu of the up-type quarks
 * from the Standard Model up-type quark masses (ignoring threshold
 * corrections).
 */
void CNE6SSM_initial_guesser<Semianalytic>::calculate_Yu_DRbar()
{
   Eigen::Matrix<double,3,3> topDRbar(Eigen::Matrix<double,3,3>::Zero());
   topDRbar(0,0) = mu_guess;
   topDRbar(1,1) = mc_guess;
   topDRbar(2,2) = mt_guess;

   const auto vu = MODELPARAMETER(vu);
   MODEL->set_Yu(Diag((1.4142135623730951*topDRbar)/vu));

}

/**
 * Calculates the Yukawa couplings Yd of the down-type
 * quarks from the Standard Model down-type quark masses (ignoring
 * threshold corrections).
 */
void CNE6SSM_initial_guesser<Semianalytic>::calculate_Yd_DRbar()
{
   Eigen::Matrix<double,3,3> bottomDRbar(Eigen::Matrix<double,3,3>::Zero());
   bottomDRbar(0,0) = md_guess;
   bottomDRbar(1,1) = ms_guess;
   bottomDRbar(2,2) = mb_guess;

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Yd(Diag((1.4142135623730951*bottomDRbar)/vd));

}

/**
 * Calculates the Yukawa couplings Ye of the leptons
 * from the Standard Model down-type lepton masses (ignoring threshold
 * corrections).
 */
void CNE6SSM_initial_guesser<Semianalytic>::calculate_Ye_DRbar()
{
   Eigen::Matrix<double,3,3> electronDRbar(Eigen::Matrix<double,3,3>::Zero());
   electronDRbar(0,0) = me_guess;
   electronDRbar(1,1) = mm_guess;
   electronDRbar(2,2) = mtau_guess;

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Ye(Diag((1.4142135623730951*electronDRbar)/vd));

}

/**
 * Guesses the SUSY parameters at the high-scale. At first it runs
 * to the guess of the high-scale (HighScaleFirstGuess) and imposes the
 * high-scale constraint (HighScaleInput):
 *
 * \code{.cpp}

 * \endcode
 *
 * It then performs a two-scale iteration at tree level to obtain an
 * estimate for the high-scale parameters at the low scale.
 */
void CNE6SSM_initial_guesser<Semianalytic>::guess_high_scale_parameters()
{
   const double low_scale_guess = low_constraint.get_initial_scale_guess();
   const double high_scale_guess = high_constraint.get_initial_scale_guess();

   model->run_to(high_scale_guess, running_precision);

   // apply high-scale constraint
   high_constraint.set_model(model);
   high_constraint.apply();

   model->run_to(low_scale_guess, running_precision);
}

/**
 * Guesses the SUSY parameters at the low- and high-scales. A 
 * two-scale iteration is performed at tree level to find
 * a consistent set of SUSY parameters at each scale.
 */
void CNE6SSM_initial_guesser<Semianalytic>::guess_susy_parameters()
{
   // ensure tree level iteration, with no thresholds
   const unsigned threshold_loop_order = model->get_thresholds();

   model->set_thresholds(0);

   high_constraint.clear();
   low_constraint.clear();

   high_constraint.set_model(model);
   low_constraint.set_model(model);

   low_constraint.set_sm_parameters(oneset);

   high_constraint.initialize();
   low_constraint.initialize();

   std::vector<Constraint<Semianalytic>*> inner_upward_constraints(2);
   inner_upward_constraints[0] = &low_constraint;
   inner_upward_constraints[1] = &high_constraint;

   std::vector<Constraint<Semianalytic>*> inner_downward_constraints(2);
   inner_downward_constraints[0] = &high_constraint;
   inner_downward_constraints[1] = &low_constraint;

   std::vector<Constraint<Semianalytic>*> outer_upward_constraints(1);
   outer_upward_constraints[0] = &susy_constraint;

   std::vector<Constraint<Semianalytic>*> outer_downward_constraints(2);
   outer_downward_constraints[0] = &susy_constraint;
   outer_downward_constraints[1] = &low_constraint;

   CNE6SSM_convergence_tester<Semianalytic> convergence_tester(model, running_precision);
   if (max_iterations > 0)
      convergence_tester.set_max_iterations(max_iterations);

   Susy_parameters_guesser susy_guesser;
   susy_guesser.initial_guesser = this;

   Two_scale_constant_precision precision(running_precision);

   RGFlow<Semianalytic> solver;
   solver.set_convergence_tester(&convergence_tester);
   solver.set_running_precision(&precision);
   solver.set_initial_guesser(&susy_guesser);
   solver.add_model(model, inner_upward_constraints, inner_downward_constraints,
                    outer_upward_constraints, outer_downward_constraints);
   solver.solve_inner_iteration_only(true);

   try {
      solver.solve();
   } catch (const NoConvergenceError&) {
#ifdef ENABLE_VERBOSE
      WARNING("CNE6SSM_initial_guesser<Semianalytic>::guess_susy_parameters(): "
              "tree level iteration has not converged!");
#endif
      // do initial guess with no iteration
      high_constraint.initialize();
      guess_low_scale_parameters();
      guess_high_scale_parameters();
   } catch (const NonPerturbativeRunningError&) {
#ifdef ENABLE_VERBOSE
      WARNING("CNE6SSM_initial_guesser<Semianalytic>::guess_susy_parameters(): "
              "non-perturbative running in tree level iteration!");
#endif
      // do initial guess with no iteration
      high_constraint.initialize();
      guess_low_scale_parameters();
      guess_high_scale_parameters();
   } catch (const NoRhoConvergenceError&) {
#ifdef ENABLE_VERBOSE
      WARNING("CNE6SSM_initial_guesser<Semianalytic>::guess_susy_parameters(): "
              "no rho convergence in tree level iteration!");
#endif
      // do initial guess with no iteration
      high_constraint.initialize();
      guess_low_scale_parameters();
      guess_high_scale_parameters();
   } catch (const Error& error) {
#ifdef ENABLE_VERBOSE
      WARNING("CNE6SSM_initial_guesser<Semianalytic>::guess_susy_parameters(): "
              "error encountered in tree level iteration!");
#endif
      // do initial guess with no iteration
      high_constraint.initialize();
      guess_low_scale_parameters();
      guess_high_scale_parameters();
   } catch (const std::string& str) {
#ifdef ENABLE_VERBOSE
      WARNING("CNE6SSM_initial_guesser<Semianalytic>::guess_susy_parameters(): "
              "error encountered in tree level iteration!");
#endif
      // do initial guess with no iteration
      high_constraint.initialize();
      guess_low_scale_parameters();
      guess_high_scale_parameters();
   } catch (const char* str) {
#ifdef ENABLE_VERBOSE
      WARNING("CNE6SSM_initial_guesser<Semianalytic>::guess_susy_parameters(): "
              "error encountered in tree level iteration!");
#endif
      // do initial guess with no iteration
      high_constraint.initialize();
      guess_low_scale_parameters();
      guess_high_scale_parameters();
   } catch (const std::exception& error) {
#ifdef ENABLE_VERBOSE
      WARNING("CNE6SSM_initial_guesser<Semianalytic>::guess_susy_parameters(): "
              "error encountered in tree level iteration!");
#endif
      // do initial guess with no iteration
      high_constraint.initialize();
      guess_low_scale_parameters();
      guess_high_scale_parameters();
   }

   model->set_thresholds(threshold_loop_order);
}

/**
 * Guesses the soft parameters at the low-scale. First the coefficients
 * in the semianalytic expansions of the soft parameters are computed
 * at the low-scale. Then, EWSB is solved at tree level at the low-scale
 * to obtain the values of the input parameters.
 */
void CNE6SSM_initial_guesser<Semianalytic>::guess_soft_parameters()
{
   const double high_scale_guess = high_constraint.get_scale();

   // calculate semianalytic coefficients
   model->calculate_coefficients(high_scale_guess);

   // apply EWSB constraint
   model->solve_ewsb_tree_level();

   // calculate tree-level spectrum
   model->calculate_DRbar_masses();
}

} // namespace flexiblesusy
