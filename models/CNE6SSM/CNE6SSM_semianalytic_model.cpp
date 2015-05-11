// ====================================================================
// Test implementation for a class to solve EWSB and do spectrum
// calculation using semianalytic algorithm
// ====================================================================

/**
 * @file CNE6SSM_semianalytic_model.cpp
 * @brief implementation of the CNE6SSM semianalytic model class
 *
 * Contains the definition of the CNE6SSM semianalytic model class
 * methods which solve EWSB and calculate pole masses and mixings
 * from DRbar parameters.
 */

#include "CNE6SSM_semianalytic_model.hpp"
#include "numerics.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "root_finder.hpp"
#include "fixed_point_iterator.hpp"
#include "gsl_utils.hpp"
#include "config.h"
#include "functors.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

#include <gsl/gsl_multiroots.h>

namespace flexiblesusy {

using namespace CNE6SSM_info;

#define CLASSNAME CNE6SSM<Semianalytic>

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model->get_input().parameter
#define LOCALINPUT(parameter) input.parameter

CLASSNAME::CNE6SSM(const CNE6SSM_input_parameters<Semianalytic>& input_)
   : Semianalytic_model()
   , CNE6SSM_mass_eigenstates()
   , input(input_)
   , number_of_ewsb_iterations(200)
   , ewsb_loop_order(2)
   , precision(1.0e-3)
   , ewsb_iteration_precision(1.0e-5)
   , ewsb_solution(Eigen::Array<double,number_of_ewsb_equations,1>::Zero())
   , TYd_Azero_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TYd_m12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , ThE_Azero_coeff(Eigen::Matrix<double,3,2>::Zero())
   , ThE_m12_coeff(Eigen::Matrix<double,3,2>::Zero())
   , TYe_Azero_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TYe_m12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TSigmaL_Azero_coeff(0), TSigmaL_m12_coeff(0)
   , TKappaPr_Azero_coeff(0), TKappaPr_m12_coeff(0)
   , TSigmax_Azero_coeff(0), TSigmax_m12_coeff(0)
   , TgD_Azero_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TgD_m12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TKappa_Azero_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TKappa_m12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TLambda12_Azero_coeff(Eigen::Matrix<double,2,2>::Zero())
   , TLambda12_m12_coeff(Eigen::Matrix<double,2,2>::Zero())
   , TLambdax_Azero_coeff(0), TLambdax_m12_coeff(0)
   , Tfu_Azero_coeff(Eigen::Matrix<double,3,2>::Zero())
   , Tfu_m12_coeff(Eigen::Matrix<double,3,2>::Zero())
   , Tfd_Azero_coeff(Eigen::Matrix<double,3,2>::Zero())
   , Tfd_m12_coeff(Eigen::Matrix<double,3,2>::Zero())
   , TYu_Azero_coeff(Eigen::Matrix<double,3,3>::Zero())
   , TYu_m12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mq2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mq2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mq2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mq2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , ml2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , ml2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , ml2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , ml2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mHd2_m02_coeff(0), mHd2_m122_coeff(0), mHd2_Azerom12_coeff(0)
   , mHd2_Azero2_coeff(0), mHu2_m02_coeff(0), mHu2_m122_coeff(0)
   , mHu2_Azerom12_coeff(0), mHu2_Azero2_coeff(0)
   , md2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , md2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , md2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , md2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mu2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mu2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mu2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mu2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , me2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , me2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , me2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , me2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , ms2_m02_coeff(0), ms2_m122_coeff(0), ms2_Azerom12_coeff(0)
   , ms2_Azero2_coeff(0), msbar2_m02_coeff(0), msbar2_m122_coeff(0)
   , msbar2_Azerom12_coeff(0), msbar2_Azero2_coeff(0)
   , mH1I2_m02_coeff(Eigen::Matrix<double,2,2>::Zero())
   , mH1I2_m122_coeff(Eigen::Matrix<double,2,2>::Zero())
   , mH1I2_Azerom12_coeff(Eigen::Matrix<double,2,2>::Zero())
   , mH1I2_Azero2_coeff(Eigen::Matrix<double,2,2>::Zero())
   , mH2I2_m02_coeff(Eigen::Matrix<double,2,2>::Zero())
   , mH2I2_m122_coeff(Eigen::Matrix<double,2,2>::Zero())
   , mH2I2_Azerom12_coeff(Eigen::Matrix<double,2,2>::Zero())
   , mH2I2_Azero2_coeff(Eigen::Matrix<double,2,2>::Zero())
   , mSI2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mSI2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mSI2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mSI2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mDx2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mDx2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mDx2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mDx2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mDxbar2_m02_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mDxbar2_m122_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mDxbar2_Azerom12_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mDxbar2_Azero2_coeff(Eigen::Matrix<double,3,3>::Zero())
   , mHp2_m02_coeff(0), mHp2_m122_coeff(0), mHp2_Azerom12_coeff(0)
   , mHp2_Azero2_coeff(0), mHpbar2_m02_coeff(0), mHpbar2_m122_coeff(0)
   , mHpbar2_Azerom12_coeff(0), mHpbar2_Azero2_coeff(0), mphi2_m02_coeff(0)
   , mphi2_m122_coeff(0), mphi2_Azerom12_coeff(0), mphi2_Azero2_coeff(0)
   , MassB_Azero_coeff(0), MassB_m12_coeff(0), MassWB_Azero_coeff(0)
   , MassWB_m12_coeff(0), MassG_Azero_coeff(0), MassG_m12_coeff(0)
   , MassBp_Azero_coeff(0), MassBp_m12_coeff(0)
{
}

CLASSNAME::~CNE6SSM()
{
}

void CLASSNAME::set_ewsb_loop_order(unsigned loop_order)
{
   ewsb_loop_order = loop_order;
}

void CLASSNAME::set_number_of_ewsb_iterations(std::size_t iterations)
{
   number_of_ewsb_iterations = iterations;
}

void CLASSNAME::set_precision(double precision_)
{
   precision = precision_;
   ewsb_iteration_precision = precision_;
   diagonalization_precision = precision_;
}

void CLASSNAME::set_ewsb_iteration_precision(double precision)
{
   ewsb_iteration_precision = precision;
}

double CLASSNAME::get_ewsb_iteration_precision() const
{
   return ewsb_iteration_precision;
}

double CLASSNAME::get_ewsb_loop_order() const
{
   return ewsb_loop_order;
}

const CNE6SSM_input_parameters<Semianalytic>& CLASSNAME::get_input() const
{
   return input;
}

void CLASSNAME::set_input_parameters(const CNE6SSM_input_parameters<Semianalytic>& input_)
{
   input = input_;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @param tadpole array of tadpole
 */
void CLASSNAME::tadpole_equations(double tadpole[number_of_ewsb_equations]) const
{
   tadpole[0] = get_ewsb_eq_hh_1();
   tadpole[1] = get_ewsb_eq_hh_2();
   tadpole[2] = get_ewsb_eq_hh_3();
   tadpole[3] = get_ewsb_eq_hh_4();
   tadpole[4] = get_ewsb_eq_hh_5();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh(0));
      tadpole[1] -= Re(tadpole_hh(1));
      tadpole[2] -= Re(tadpole_hh(2));
      tadpole[3] -= Re(tadpole_hh(3));
      tadpole[4] -= Re(tadpole_hh(4));

      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[3];
         tadpole_hh_2loop(two_loop_tadpole);
         tadpole[0] -= two_loop_tadpole[0];
         tadpole[1] -= two_loop_tadpole[1];
         tadpole[2] -= two_loop_tadpole[2];

      }
   }
}

/**
 * Method which calculates the tadpoles at loop order specified in the
 * pointer to the CLASSNAME::EWSB_args struct.
 *
 * @param x GSL vector of EWSB output parameters
 * @param params pointer to CLASSNAME::EWSB_args struct
 * @param f GSL vector with tadpoles
 *
 * @return GSL_EDOM if x contains Nans, GSL_SUCCESS otherwise.
 */
int CLASSNAME::tadpole_equations(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (contains_nan(x, number_of_ewsb_equations)) {
      for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
         gsl_vector_set(f, i, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::EWSB_args* ewsb_args
      = static_cast<CLASSNAME::EWSB_args*>(params);
   CNE6SSM* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   // EWSB output parameters are TanTheta, m0, vphi, XiF, LXiF
   const double s = model->get_input().sInput;
   const double m12 = model->get_input().m12;
   const double Azero = model->get_input().Azero;

   const double m0 = gsl_vector_get(x, 0);
   const double TanTheta = gsl_vector_get(x, 1);
   const double vphi = gsl_vector_get(x, 2);
   const double XiF = gsl_vector_get(x, 3);
   const double LXiF = gsl_vector_get(x, 4);

   model->set_soft_parameters_at_current_scale(m0, m12, Azero);
   model->set_vs(s * Cos(ArcTan(TanTheta)));
   model->set_vsb(s * Sin(ArcTan(TanTheta)));
   model->set_vphi(vphi);
   model->set_XiF(XiF);
   model->set_LXiF(LXiF);

   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   double tadpole[number_of_ewsb_equations] = { 0. };

   model->tadpole_equations(tadpole);

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, tadpole[i]);

   bool is_finite = true;

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      is_finite = is_finite && std::isfinite(tadpole[i]);

   return (is_finite ? GSL_SUCCESS : GSL_EDOM);
}

/**
 * This method solves the EWSB conditions iteratively, trying several
 * root finding methods until a solution is found.
 */
int CLASSNAME::solve_ewsb_iteratively()
{
   EWSB_args params = {this, ewsb_loop_order};

   EWSB_solver* solvers[] = {
      new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative>(CLASSNAME::ewsb_step, &params, number_of_ewsb_iterations, ewsb_iteration_precision),
      new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_hybrid),
      new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_hybrids),
      new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_broyden)
   };

   const std::size_t number_of_solvers = sizeof(solvers)/sizeof(*solvers);
   double x_init[number_of_ewsb_equations];

#ifdef ENABLE_VERBOSE
   std::cout << "Getting EWSB initial guess ...\n";
#endif

   ewsb_initial_guess(x_init);

#ifdef ENABLE_VERBOSE
   std::cout << "Solving EWSB equations ...\n"
      "\tInitial guess: x_init =";
   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      std::cout << ' ' << x_init[i];
   std::cout << '\n';
#endif

   int status;
   for (std::size_t i = 0; i < number_of_solvers; ++i) {
      VERBOSE_MSG("\tStarting EWSB iteration using solver " << i);
      status = solve_ewsb_iteratively_with(solvers[i], x_init);
      if (status == EWSB_solver::SUCCESS) {
         VERBOSE_MSG("\tSolver " << i << " finished successfully!");
         break;
      }
#ifdef ENABLE_VERBOSE
      else {
         WARNING("\tSolver " << i << " could not find a solution!"
                 " (requested precision: " << ewsb_iteration_precision << ")");
      }
#endif
   }

   if (status == EWSB_solver::SUCCESS) {
      problems.unflag_no_ewsb();
   } else {
      problems.flag_no_ewsb();
#ifdef ENABLE_VERBOSE
      WARNING("\tCould not find a solution to the EWSB equations!"
              " (requested precision: " << ewsb_iteration_precision << ")");
#endif
   }

   for_each(solvers, solvers + number_of_solvers, Delete_object());

   return status;
}

/**
 * Solves EWSB equations with given EWSB solver
 *
 * @param solver EWSB solver
 * @param x_init initial values
 *
 * @return status of the EWSB solver
 */
int CLASSNAME::solve_ewsb_iteratively_with(
   EWSB_solver* solver,
   const double x_init[number_of_ewsb_equations]
)
{
   const int status = solver->solve(x_init);

   const double s = LOCALINPUT(sInput);
   const double m12 = LOCALINPUT(m12);
   const double Azero = LOCALINPUT(Azero);

   ewsb_solution(0) = solver->get_solution(0);
   ewsb_solution(1) = solver->get_solution(1);
   ewsb_solution(2) = solver->get_solution(2);
   ewsb_solution(3) = solver->get_solution(3);
   ewsb_solution(4) = solver->get_solution(4);

   set_soft_parameters_at_current_scale(ewsb_solution(0), m12, Azero);

   vs = s * Cos(ArcTan(ewsb_solution(1)));
   vsb = s * Sin(ArcTan(ewsb_solution(1)));
   vphi = ewsb_solution(2);
   XiF = ewsb_solution(3);
   LXiF = ewsb_solution(4);


   return status;
}

int CLASSNAME::check_ewsb_solution(double precision)
{
   double tadpole[number_of_ewsb_equations];

   if (ewsb_loop_order > 0) {
      calculate_DRbar_masses();
   }

   tadpole_equations(tadpole);

   double residual = Abs(tadpole[0]);

   for (std::size_t i = 1; i < number_of_ewsb_equations; ++i) {
      residual += Abs(tadpole[i]);
   } 
   
   return (residual < precision ? EWSB_solver::SUCCESS : EWSB_solver::FAIL);
}

int CLASSNAME::solve_ewsb_iteratively(unsigned loop_order)
{
   // temporarily set `ewsb_loop_order' to `loop_order' and do
   // iteration
   const unsigned old_loop_order = ewsb_loop_order;
   ewsb_loop_order = loop_order;
   const int status = solve_ewsb_iteratively();
   ewsb_loop_order = old_loop_order;
   return status;
}


int CLASSNAME::solve_ewsb_tree_level()
{
   int error = 0;

   error = solve_ewsb_iteratively(0);


   return error;
}

int CLASSNAME::solve_ewsb_one_loop()
{
   return solve_ewsb_iteratively(1);
}

int CLASSNAME::solve_ewsb()
{
   VERBOSE_MSG("\tSolving EWSB at " << ewsb_loop_order << "-loop order");

   if (ewsb_loop_order == 0)
      return solve_ewsb_tree_level();

   return solve_ewsb_iteratively(ewsb_loop_order);
}

void CLASSNAME::ewsb_initial_guess(double x_init[number_of_ewsb_equations])
{
   const double s = LOCALINPUT(sInput);
   const double m12 = LOCALINPUT(m12);
   const double Azero = LOCALINPUT(Azero);

   // initially guess m0, assuming TanTheta = 1 (or approximately so)
   x_init[0] = (mHd2_m122_coeff * Sqr(vd) - mHu2_m122_coeff * Sqr(vu)) * Sqr(m12)
      + (mHd2_Azerom12_coeff * Sqr(vd) - mHu2_Azerom12_coeff * Sqr(vu)) * m12 * Azero
      + (mHd2_Azero2_coeff * Sqr(vd) - mHu2_Azero2_coeff * Sqr(vu)) * Sqr(Azero)
      + 0.25 * Sqr(Lambdax) * Sqr(s) * (Sqr(vd) - Sqr(vu)) +
      0.125 * Sqr(g2) * (Power(vd, 4) - Power(vu, 4)) + 0.075 * Sqr(g1) *
      (Power(vd, 4) - Power(vu, 4)) + 0.0125 * Sqr(g1p) *
      (9.0 * Power(vd, 4) - 4.0 * Power(vu, 4));
   
   x_init[0] /= (mHu2_m02_coeff * Sqr(vu) - mHd2_m02_coeff * Sqr(vd));

#ifdef ENABLE_VERBOSE
   if (x_init[0] < 0.) {
      WARNING("\tSqr(m0) < 0 in initial guess!"); 
   }
#endif
   x_init[0] = AbsSqrt(x_init[0]);

   // set soft masses using initial guess for m0
   const double m0 = x_init[0];
   double mHd2 = mHd2_m02_coeff * Sqr(m0) + mHd2_m122_coeff * Sqr(m12)
      + mHd2_Azerom12_coeff * Azero * m12 + mHd2_Azero2_coeff * Sqr(Azero);
   double ms2 = ms2_m02_coeff * Sqr(m0) + ms2_m122_coeff * Sqr(m12)
      + ms2_Azerom12_coeff * Azero * m12 + ms2_Azero2_coeff * Sqr(Azero);
   double msbar2 = msbar2_m02_coeff * Sqr(m0) + msbar2_m122_coeff * Sqr(m12)
      + msbar2_Azerom12_coeff * Azero * m12 + msbar2_Azero2_coeff * Sqr(Azero);
   double mphi2 = mphi2_m02_coeff * Sqr(m0) + mphi2_m122_coeff * Sqr(m12)
      + mphi2_Azerom12_coeff * Azero * m12 + mphi2_Azero2_coeff * Sqr(Azero);

   // initial guess for TanTheta
   x_init[1] = (ms2 + 0.0125 * Sqr(g1p) * Sqr(QS) * Sqr(s)) 
      / (msbar2 + 0.0125 * Sqr(g1p) * Sqr(QS) * Sqr(s));

#ifdef ENABLE_VERBOSE
   if (x_init[1] < 0.) {
      WARNING("\tSqr(TanTheta) < 0 in initial guess!");
   }
#endif
   x_init[1] = AbsSqrt(x_init[1]);

   double vs = s * Cos(ArcTan(x_init[1]));
   double vsb = s * Sin(ArcTan(x_init[1]));

   // initial guess for vphi
   x_init[2] = (-4. / (vsb * vu * Lambdax * Conj(Sigmax) +
                       vsb * vu * Sigmax * Conj(Lambdax))) *
      (mHd2 * vd - 0.35355339059327373 * vs * vu * TLambdax
       - 0.35355339059327373 * vs * vu * Conj(TLambdax) + 0.5 * AbsSqr(Lambdax) * vd *
       (Sqr(vu) + Sqr(vs)) + 0.125 * Sqr(g2) * Power(vd,3) + 0.075 * Sqr(g1) * Power(vd,3)
       - 0.125 * Sqr(g2) * vd * Sqr(vu) - 0.075 * Sqr(g1) * vd * Sqr(vu) + 0.1125 * Sqr(g1p) *
       Power(vd,3) + 0.075 * Sqr(g1p) * vd * Sqr(vu) - 0.0375 * Sqr(g1p) * QS * vd *
       (Sqr(vs) - Sqr(vsb)));

   double vphi = x_init[2];

   // initial guess for XiF
   x_init[3] = (2. / (vs * Sigmax + vs * Conj(Sigmax))) *
      (msbar2 * vsb - 0.35355339059327373 * vphi * vs * MuPhi * Conj(Sigmax)
       - 0.35355339059327373 * vphi * vs * Sigmax * Conj(MuPhi) - 0.35355339059327373 * vphi *
       vs * TSigmax - 0.35355339059327373 * vphi * vs * Conj(TSigmax) + 0.25 * vphi * vd * vu *
       Lambdax * Conj(Sigmax) + 0.25 * vphi * vd * vu * Sigmax * Conj(Lambdax) 
       + 0.5 * AbsSqr(Sigmax) * vsb * Sqr(vphi) + 0.5 * AbsSqr(Sigmax) * vsb * Sqr(vs) - 0.25 * vs
       * Sqr(vphi) * KappaPr * Conj(Sigmax) - 0.25 * vs * Sqr(vphi) * Sigmax * Conj(KappaPr)
       + 0.0375 * Sqr(g1p) * QS * vsb * Sqr(vd) + 0.025 * Sqr(g1p) * QS * vsb * Sqr(vu) - 0.0125 *
       Sqr(g1p) * Sqr(QS) * vsb * Sqr(vs) + 0.0125 * Sqr(g1p) * Sqr(QS) * Power(vsb,3));

   double XiF = x_init[3];

   // initial guess for LXiF
   x_init[4] = -0.7071067811865475 *
      (mphi2 * vphi + vphi * AbsSqr(MuPhi) + Power(vphi,3) * AbsSqr(
         KappaPr) + 0.5 * vphi * BMuPhi + 0.5 * vphi * Conj(BMuPhi) + 0.7071067811865475 *
       MuPhi * Conj(XiF) - 0.35355339059327373 * MuPhi * vs * vsb * Conj(Sigmax) 
       - 0.35355339059327373 * vs * vsb * Conj(TSigmax) +
       vphi * Conj(XiF) * KappaPr - 0.5 * vphi * vs * vsb * Conj(Sigmax) *KappaPr + 0.25 * vd * vsb *
       vu * Conj(Sigmax) * Lambdax + 0.7071067811865475 * Conj(MuPhi) * XiF + vphi * Conj(
          KappaPr) * XiF - 0.35355339059327373 * vs * vsb * Conj(MuPhi) * Sigmax - 0.5 * vphi *
       vs * vsb * Conj(KappaPr) * Sigmax + 0.25 * vd * vsb * vu * Conj(Lambdax) * Sigmax
       + 1.0606601717798212 * MuPhi * Conj(KappaPr) * Sqr(vphi) +
       0.35355339059327373 * Conj(TKappaPr) * Sqr(vphi) + 1.0606601717798212 * Conj(
          MuPhi) * KappaPr * Sqr(vphi) + 0.5 * vphi * AbsSqr(Sigmax) * Sqr(vs) + 0.5 * vphi * AbsSqr
       (Sigmax) * Sqr(vsb) + 0.35355339059327373 * Sqr(vphi) * TKappaPr -
       0.35355339059327373 * vs * vsb *TSigmax);
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * @param ewsb_parameters new EWSB output parameters.  \a
 * ewsb_parameters is only modified if all new parameters are finite.
 *
 * @return GSL_SUCCESS if new EWSB output parameters are finite,
 * GSL_EDOM otherwise.
 */
int CLASSNAME::ewsb_step(double ewsb_parameters[number_of_ewsb_equations])
{

   // save initial parameters
   CNE6SSM_soft_parameters saved_pars;
   saved_pars.set(get());

   int error;

   const double s = LOCALINPUT(sInput);
   const double m12 = LOCALINPUT(m12);
   const double Azero = LOCALINPUT(Azero);

   // update m0
   double m0 = (mHd2_m122_coeff * Sqr(vd) - mHu2_m122_coeff * Sqr(vu)) * Sqr(m12)
      + (mHd2_Azerom12_coeff * Sqr(vd) - mHu2_Azerom12_coeff * Sqr(vu)) * m12 * Azero
      + (mHd2_Azero2_coeff * Sqr(vd) - mHu2_Azero2_coeff * Sqr(vu)) * Sqr(Azero)
      + 0.25 * Sqr(Lambdax) * Sqr(vs) * (Sqr(vd) - Sqr(vu)) +
      0.125 * Sqr(g2) * (Power(vd, 4) - Power(vu, 4)) + 0.075 * Sqr(g1) *
      (Power(vd, 4) - Power(vu, 4)) + 0.0125 * Sqr(g1p) *
      (9.0 * Power(vd, 4) - 4.0 * Power(vu, 4)) + 0.0125 * Sqr(g1p) * QS *
      (-3.0 * Sqr(vd) - 2.0 * Sqr(vu)) * (Sqr(vs) - Sqr(vsb));

   if (ewsb_loop_order > 0) {
      m0 -= (vd * Re(tadpole_hh(0)) - vu * Re(tadpole_hh(1)));
      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[3];
         tadpole_hh_2loop(two_loop_tadpole);
         m0 -= (vd * two_loop_tadpole[0] - vu * two_loop_tadpole[1]);
      }
   }

   m0 = AbsSqrt(m0 / (mHu2_m02_coeff * Sqr(vu) - mHd2_m02_coeff * Sqr(vd)));

   set_soft_parameters_at_current_scale(m0, m12, Azero);

   // update TanTheta
   double delta = 0.5 * AbsSqr(Lambdax) * Sqr(vs) * (Sqr(vd) + Sqr(vu))
      + 0.5 * AbsSqr(Sigmax) * Sqr(vs) * Sqr(vphi) - 0.5 * AbsSqr(Sigmax) * Sqr(vsb) * Sqr(vphi)
      - 0.35355339059327373 * vd * vu * vs * TLambdax - 0.35355339059327373 * vd * vu * vs * Conj(TLambdax) 
      - 0.25 * vphi * vsb * vd * vu * Lambdax * Conj(Sigmax) - 0.25 * vphi * vsb * vd * vu * Sigmax * Conj(Lambdax)
      - 0.0375 * QS * Sqr(g1p) * Sqr(s) * Sqr(vd) - 0.025 * QS * Sqr(g1p) * Sqr(s) * Sqr(vu);

   if (ewsb_loop_order > 0) {
      delta -= (vs * Re(tadpole_hh(2)) - vsb * Re(tadpole_hh(3)));
      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[3];
         tadpole_hh_2loop(two_loop_tadpole);
         delta -= vs * two_loop_tadpole[2];
      }
   }

   const double TanTheta = AbsSqrt((ms2 + 0.0125 * Sqr(g1p) * Sqr(QS) * Sqr(s) + delta / (Sqr(vs)))
                                   / (msbar2 + 0.0125 * Sqr(g1p) * Sqr(QS) * Sqr(s)));

   vs = s * Cos(ArcTan(TanTheta));
   vsb = s * Sin(ArcTan(TanTheta));

   // update vphi
   double rhs_vphi = mHd2 * vd - 0.35355339059327373 * vs * vu * TLambdax
      - 0.35355339059327373 * vs * vu * Conj(TLambdax) + 0.5 * AbsSqr(Lambdax) * vd *
      (Sqr(vu) + Sqr(vs)) + 0.125 * Sqr(g2) * Power(vd,3) + 0.075 * Sqr(g1) * Power(vd,3)
      - 0.125 * Sqr(g2) * vd * Sqr(vu) - 0.075 * Sqr(g1) * vd * Sqr(vu) + 0.1125 * Sqr(g1p) *
      Power(vd,3) + 0.075 * Sqr(g1p) * vd * Sqr(vu) - 0.0375 * Sqr(g1p) * QS * vd *
      (Sqr(vs) - Sqr(vsb));

   if (ewsb_loop_order > 0) {
      rhs_vphi -= Re(tadpole_hh(0));
      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[3];
         tadpole_hh_2loop(two_loop_tadpole);
         rhs_vphi -= two_loop_tadpole[0];
      }
   }

   rhs_vphi *= (-4. / (vsb * vu * Lambdax * Conj(Sigmax) + vsb * vu * Sigmax * Conj(Lambdax)));

   vphi = rhs_vphi;

   // update XiF
   double rhs_XiF = msbar2 * vsb - 0.35355339059327373 * vphi * vs * MuPhi * Conj(Sigmax)
      - 0.35355339059327373 * vphi * vs * Sigmax * Conj(MuPhi) - 0.35355339059327373 * vphi *
      vs * TSigmax - 0.35355339059327373 * vphi * vs * Conj(TSigmax) + 0.25 * vphi * vd * vu *
      Lambdax * Conj(Sigmax) + 0.25 * vphi * vd * vu * Sigmax * Conj(Lambdax) 
      + 0.5 * AbsSqr(Sigmax) * vsb * Sqr(vphi) + 0.5 * AbsSqr(Sigmax) * vsb * Sqr(vs) - 0.25 * vs
      * Sqr(vphi) * KappaPr * Conj(Sigmax) - 0.25 * vs * Sqr(vphi) * Sigmax * Conj(KappaPr)
      + 0.0375 * Sqr(g1p) * QS * vsb * Sqr(vd) + 0.025 * Sqr(g1p) * QS * vsb * Sqr(vu) - 0.0125 *
      Sqr(g1p) * Sqr(QS) * vsb * Sqr(vs) + 0.0125 * Sqr(g1p) * Sqr(QS) * Power(vsb,3);

   if (ewsb_loop_order > 0) {
      rhs_XiF -= Re(tadpole_hh(3));
      if (ewsb_loop_order > 1) {

      }
   }

   rhs_XiF *= ( 2. / (vs * Sigmax + vs * Conj(Sigmax)));

   XiF = rhs_XiF;

   // update LXiF
   double rhs_LXiF = mphi2 * vphi + vphi * AbsSqr(MuPhi) + Power(vphi,3) * AbsSqr(
      KappaPr) + 0.5 * vphi * BMuPhi + 0.5 * vphi * Conj(BMuPhi) + 0.7071067811865475 *
      MuPhi * Conj(XiF) - 0.35355339059327373 * MuPhi * vs * vsb * Conj(Sigmax) 
      - 0.35355339059327373 * vs * vsb * Conj(TSigmax) +
      vphi * Conj(XiF) * KappaPr - 0.5 * vphi * vs * vsb * Conj(Sigmax) * KappaPr + 0.25 * vd * vsb *
      vu * Conj(Sigmax) * Lambdax + 0.7071067811865475 * Conj(MuPhi) * XiF + vphi * Conj(
      KappaPr) * XiF - 0.35355339059327373 * vs * vsb * Conj(MuPhi) * Sigmax - 0.5 * vphi * vs *
      vsb * Conj(KappaPr) * Sigmax + 0.25 * vd * vsb * vu * Conj(Lambdax) * Sigmax
      + 1.0606601717798212 * MuPhi * Conj(KappaPr) * Sqr(vphi) +
      0.35355339059327373 * Conj(TKappaPr) * Sqr(vphi) + 1.0606601717798212 * Conj(
      MuPhi) * KappaPr * Sqr(vphi) + 0.5 * vphi * AbsSqr(Sigmax) * Sqr(vs) + 0.5 * vphi * AbsSqr
      (Sigmax) * Sqr(vsb) + 0.35355339059327373 * Sqr(vphi) * TKappaPr -
      0.35355339059327373 * vs * vsb * TSigmax;

   if (ewsb_loop_order > 0) {
      rhs_LXiF -= Re(tadpole_hh(4));
      if (ewsb_loop_order > 1) {

      }
   }

   rhs_LXiF *= -0.7071067811865475;

   LXiF = rhs_LXiF;

   const bool is_finite = std::isfinite(m0) && std::isfinite(TanTheta)
      && std::isfinite(vphi) && std::isfinite(XiF) && std::isfinite(LXiF);

   if (is_finite) {
      error = GSL_SUCCESS;

      ewsb_parameters[0] = m0;
      ewsb_parameters[1] = TanTheta;
      ewsb_parameters[2] = vphi;
      ewsb_parameters[3] = XiF;
      ewsb_parameters[4] = LXiF;

   } else {
      error = GSL_EDOM;
   }

   // reset old parameters
   set(saved_pars.get());

   return error;
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * @param x old EWSB output parameters
 * @param params further function parameters
 * @param f new EWSB output parameters
 *
 * @return Returns status of CLASSNAME::ewsb_step
 */
int CLASSNAME::ewsb_step(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (contains_nan(x, number_of_ewsb_equations)) {
      for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
         gsl_vector_set(f, i, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::EWSB_args* ewsb_args
      = static_cast<CLASSNAME::EWSB_args*>(params);
   CNE6SSM* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double s = INPUT(sInput);
   const double m12 = INPUT(m12);
   const double Azero = INPUT(Azero);

   const double m0 = gsl_vector_get(x, 0);
   const double TanTheta = gsl_vector_get(x, 1);
   const double vphi = gsl_vector_get(x, 2);
   const double XiF = gsl_vector_get(x, 3);
   const double LXiF = gsl_vector_get(x, 4);

   model->set_soft_parameters_at_current_scale(m0, m12, Azero);
   model->set_vs(s * Cos(ArcTan(TanTheta)));
   model->set_vsb(s * Sin(ArcTan(TanTheta)));
   model->set_vphi(vphi);
   model->set_XiF(XiF);
   model->set_LXiF(LXiF);


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   double ewsb_parameters[number_of_ewsb_equations] =
      { m0, TanTheta, vphi, XiF, LXiF};

   const int status = model->ewsb_step(ewsb_parameters);

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, ewsb_parameters[i]);

   return status;
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "CNE6SSM (solver type: semianalytic)\n"
           "========================================\n";
   CNE6SSM_mass_eigenstates::print(ostr);
}

/**
 * calculates spectrum for model once the DRbar parameters at
 * at low energies are known
 */
void CLASSNAME::calculate_spectrum()
{
   calculate_DRbar_masses();
   if (pole_mass_loop_order > 0)
      calculate_pole_masses();

   // move goldstone bosons to the front
   reorder_DRbar_masses();
   if (pole_mass_loop_order == 0)
      copy_DRbar_masses_to_pole_masses();
   else
      reorder_pole_masses();

   if (problems.have_problem() && !force_output) {
      clear_DRbar_parameters();
      physical.clear();
   }
}

void CLASSNAME::clear_problems()
{
   problems.unflag_all_tachyons();
}

void CLASSNAME::clear()
{
   CNE6SSM_mass_eigenstates::clear();
}

std::string CLASSNAME::name() const
{
   return "CNE6SSM";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   CNE6SSM_mass_eigenstates::run_to(scale, eps);
}

double CLASSNAME::get_parameter(unsigned parameter) const
{
   if (parameter >= CNE6SSM_info::NUMBER_OF_PARAMETERS)
      throw UnknownModelParameterError(parameter);

   switch (parameter) {

   case CNE6SSM_info::Yd00:
      return Yd(0,0);
   case CNE6SSM_info::Yd01:
      return Yd(0,1);
   case CNE6SSM_info::Yd02:
      return Yd(0,2);
   case CNE6SSM_info::Yd10:
      return Yd(1,0); 
   case CNE6SSM_info::Yd11:
      return Yd(1,1); 
   case CNE6SSM_info::Yd12:
      return Yd(1,2);
   case CNE6SSM_info::Yd20:
      return Yd(2,0);
   case CNE6SSM_info::Yd21:
      return Yd(2,1);
   case CNE6SSM_info::Yd22:
      return Yd(2,2);
   case CNE6SSM_info::hE00:
      return hE(0,0);
   case CNE6SSM_info::hE01:
      return hE(0,1);
   case CNE6SSM_info::hE10:
      return hE(1,0); 
   case CNE6SSM_info::hE11:
      return hE(1,1); 
   case CNE6SSM_info::hE20:
      return hE(2,0); 
   case CNE6SSM_info::hE21:
      return hE(2,1);
   case CNE6SSM_info::Ye00:
      return Ye(0,0); 
   case CNE6SSM_info::Ye01:
      return Ye(0,1);
   case CNE6SSM_info::Ye02:
      return Ye(0,2); 
   case CNE6SSM_info::Ye10:
      return Ye(1,0); 
   case CNE6SSM_info::Ye11:
      return Ye(1,1); 
   case CNE6SSM_info::Ye12:
      return Ye(1,2);
   case CNE6SSM_info::Ye20:
      return Ye(2,0); 
   case CNE6SSM_info::Ye21:
      return Ye(2,1); 
   case CNE6SSM_info::Ye22:
      return Ye(2,2); 
   case CNE6SSM_info::SigmaL:
      return SigmaL; 
   case CNE6SSM_info::KappaPr:
      return KappaPr; 
   case CNE6SSM_info::Sigmax:
      return Sigmax; 
   case CNE6SSM_info::gD00:
      return gD(0,0); 
   case CNE6SSM_info::gD01:
      return gD(0,1); 
   case CNE6SSM_info::gD02:
      return gD(0,2); 
   case CNE6SSM_info::gD10:
      return gD(1,0); 
   case CNE6SSM_info::gD11:
      return gD(1,1);
   case CNE6SSM_info::gD12:
      return gD(1,2); 
   case CNE6SSM_info::gD20:
      return gD(2,0);
   case CNE6SSM_info::gD21:
      return gD(2,1); 
   case CNE6SSM_info::gD22:
      return gD(2,2); 
   case CNE6SSM_info::Kappa00:
      return Kappa(0,0); 
   case CNE6SSM_info::Kappa01:
      return Kappa(0,1); 
   case CNE6SSM_info::Kappa02:
      return Kappa(0,2); 
   case CNE6SSM_info::Kappa10:
      return Kappa(1,0); 
   case CNE6SSM_info::Kappa11:
      return Kappa(1,1); 
   case CNE6SSM_info::Kappa12:
      return Kappa(1,2);
   case CNE6SSM_info::Kappa20:
      return Kappa(2,0); 
   case CNE6SSM_info::Kappa21:
      return Kappa(2,1); 
   case CNE6SSM_info::Kappa22:
      return Kappa(2,2); 
   case CNE6SSM_info::Lambda1200:
      return Lambda12(0,0); 
   case CNE6SSM_info::Lambda1201:
      return Lambda12(0,1); 
   case CNE6SSM_info::Lambda1210:
      return Lambda12(1,0); 
   case CNE6SSM_info::Lambda1211:
      return Lambda12(1,1);
   case CNE6SSM_info::Lambdax:
      return Lambdax;
   case CNE6SSM_info::fu00:
      return fu(0,0);
   case CNE6SSM_info::fu01:
      return fu(0,1); 
   case CNE6SSM_info::fu10:
      return fu(1,0);
   case CNE6SSM_info::fu11:
      return fu(1,1); 
   case CNE6SSM_info::fu20:
      return fu(2,0); 
   case CNE6SSM_info::fu21:
      return fu(2,1); 
   case CNE6SSM_info::fd00:
      return fd(0,0); 
   case CNE6SSM_info::fd01:
      return fd(0,1);
   case CNE6SSM_info::fd10:
      return fd(1,0); 
   case CNE6SSM_info::fd11:
      return fd(1,1); 
   case CNE6SSM_info::fd20:
      return fd(2,0);
   case CNE6SSM_info::fd21:
      return fd(2,1); 
   case CNE6SSM_info::Yu00:
      return Yu(0,0); 
   case CNE6SSM_info::Yu01:
      return Yu(0,1); 
   case CNE6SSM_info::Yu02:
      return Yu(0,2); 
   case CNE6SSM_info::Yu10:
      return Yu(1,0); 
   case CNE6SSM_info::Yu11:
      return Yu(1,1); 
   case CNE6SSM_info::Yu12:
      return Yu(1,2); 
   case CNE6SSM_info::Yu20:
      return Yu(2,0); 
   case CNE6SSM_info::Yu21:
      return Yu(2,1); 
   case CNE6SSM_info::Yu22:
      return Yu(2,2); 
   case CNE6SSM_info::MuPr:
      return MuPr;
   case CNE6SSM_info::MuPhi:
      return MuPhi; 
   case CNE6SSM_info::XiF:
      return XiF;
   case CNE6SSM_info::g1:
      return g1; 
   case CNE6SSM_info::g2:
      return g2;
   case CNE6SSM_info::g3:
      return g3; 
   case CNE6SSM_info::g1p:
      return g1p; 
   case CNE6SSM_info::vd:
      return vd; 
   case CNE6SSM_info::vu:
      return vu; 
   case CNE6SSM_info::vs:
      return vs; 
   case CNE6SSM_info::vsb:
      return vsb; 
   case CNE6SSM_info::vphi:
      return vphi; 
   case CNE6SSM_info::TYd00:
      return TYd(0,0); 
   case CNE6SSM_info::TYd01:
      return TYd(0,1); 
   case CNE6SSM_info::TYd02:
      return TYd(0,2); 
   case CNE6SSM_info::TYd10:
      return TYd(1,0); 
   case CNE6SSM_info::TYd11:
      return TYd(1,1);
   case CNE6SSM_info::TYd12:
      return TYd(1,2); 
   case CNE6SSM_info::TYd20:
      return TYd(2,0); 
   case CNE6SSM_info::TYd21:
      return TYd(2,1);
   case CNE6SSM_info::TYd22:
      return TYd(2,2); 
   case CNE6SSM_info::ThE00:
      return ThE(0,0); 
   case CNE6SSM_info::ThE01:
      return ThE(0,1); 
   case CNE6SSM_info::ThE10:
      return ThE(1,0); 
   case CNE6SSM_info::ThE11:
      return ThE(1,1); 
   case CNE6SSM_info::ThE20:
      return ThE(2,0); 
   case CNE6SSM_info::ThE21:
      return ThE(2,1); 
   case CNE6SSM_info::TYe00:
      return TYe(0,0);
   case CNE6SSM_info::TYe01:
      return TYe(0,1); 
   case CNE6SSM_info::TYe02:
      return TYe(0,2); 
   case CNE6SSM_info::TYe10:
      return TYe(1,0); 
   case CNE6SSM_info::TYe11:
      return TYe(1,1); 
   case CNE6SSM_info::TYe12:
      return TYe(1,2); 
   case CNE6SSM_info::TYe20:
      return TYe(2,0); 
   case CNE6SSM_info::TYe21:
      return TYe(2,1); 
   case CNE6SSM_info::TYe22:
      return TYe(2,2); 
   case CNE6SSM_info::TSigmaL:
      return TSigmaL; 
   case CNE6SSM_info::TKappaPr:
      return TKappaPr;
   case CNE6SSM_info::TSigmax:
      return TSigmax;
   case CNE6SSM_info::TgD00:
      return TgD(0,0); 
   case CNE6SSM_info::TgD01:
      return TgD(0,1); 
   case CNE6SSM_info::TgD02:
      return TgD(0,2); 
   case CNE6SSM_info::TgD10:
      return TgD(1,0); 
   case CNE6SSM_info::TgD11:
      return TgD(1,1); 
   case CNE6SSM_info::TgD12:
      return TgD(1,2); 
   case CNE6SSM_info::TgD20:
      return TgD(2,0); 
   case CNE6SSM_info::TgD21:
      return TgD(2,1); 
   case CNE6SSM_info::TgD22:
      return TgD(2,2);
   case CNE6SSM_info::TKappa00:
      return TKappa(0,0); 
   case CNE6SSM_info::TKappa01:
      return TKappa(0,1); 
   case CNE6SSM_info::TKappa02:
      return TKappa(0,2); 
   case CNE6SSM_info::TKappa10:
      return TKappa(1,0); 
   case CNE6SSM_info::TKappa11:
      return TKappa(1,1); 
   case CNE6SSM_info::TKappa12:
      return TKappa(1,2); 
   case CNE6SSM_info::TKappa20:
      return TKappa(2,0);
   case CNE6SSM_info::TKappa21:
      return TKappa(2,1); 
   case CNE6SSM_info::TKappa22:
      return TKappa(2,2);
   case CNE6SSM_info::TLambda1200:
      return TLambda12(0,0); 
   case CNE6SSM_info::TLambda1201:
      return TLambda12(0,1); 
   case CNE6SSM_info::TLambda1210:
      return TLambda12(1,0); 
   case CNE6SSM_info::TLambda1211:
      return TLambda12(1,1);
   case CNE6SSM_info::TLambdax:
      return TLambdax;
   case CNE6SSM_info::Tfu00:
      return Tfu(0,0); 
   case CNE6SSM_info::Tfu01:
      return Tfu(0,1); 
   case CNE6SSM_info::Tfu10:
      return Tfu(1,0);
   case CNE6SSM_info::Tfu11:
      return Tfu(1,1); 
   case CNE6SSM_info::Tfu20:
      return Tfu(2,0); 
   case CNE6SSM_info::Tfu21:
      return Tfu(2,1); 
   case CNE6SSM_info::Tfd00:
      return Tfd(0,0); 
   case CNE6SSM_info::Tfd01:
      return Tfd(0,1); 
   case CNE6SSM_info::Tfd10:
      return Tfd(1,0);
   case CNE6SSM_info::Tfd11:
      return Tfd(1,1); 
   case CNE6SSM_info::Tfd20:
      return Tfd(2,0); 
   case CNE6SSM_info::Tfd21:
      return Tfd(2,1);
   case CNE6SSM_info::TYu00:
      return TYu(0,0); 
   case CNE6SSM_info::TYu01:
      return TYu(0,1); 
   case CNE6SSM_info::TYu02:
      return TYu(0,2); 
   case CNE6SSM_info::TYu10:
      return TYu(1,0); 
   case CNE6SSM_info::TYu11:
      return TYu(1,1); 
   case CNE6SSM_info::TYu12:
      return TYu(1,2); 
   case CNE6SSM_info::TYu20:
      return TYu(2,0); 
   case CNE6SSM_info::TYu21:
      return TYu(2,1);
   case CNE6SSM_info::TYu22:
      return TYu(2,2); 
   case CNE6SSM_info::BMuPr:
      return BMuPr;
   case CNE6SSM_info::BMuPhi:
      return BMuPhi; 
   case CNE6SSM_info::LXiF:
      return LXiF; 
   case CNE6SSM_info::mq200:
      return mq2(0,0); 
   case CNE6SSM_info::mq201:
      return mq2(0,1); 
   case CNE6SSM_info::mq202:
      return mq2(0,2); 
   case CNE6SSM_info::mq210:
      return mq2(1,0); 
   case CNE6SSM_info::mq211:
      return mq2(1,1); 
   case CNE6SSM_info::mq212:
      return mq2(1,2); 
   case CNE6SSM_info::mq220:
      return mq2(2,0);
   case CNE6SSM_info::mq221:
      return mq2(2,1); 
   case CNE6SSM_info::mq222:
      return mq2(2,2); 
   case CNE6SSM_info::ml200:
      return ml2(0,0);
   case CNE6SSM_info::ml201:
      return ml2(0,1); 
   case CNE6SSM_info::ml202:
      return ml2(0,2); 
   case CNE6SSM_info::ml210:
      return ml2(1,0); 
   case CNE6SSM_info::ml211:
      return ml2(1,1); 
   case CNE6SSM_info::ml212:
      return ml2(1,2); 
   case CNE6SSM_info::ml220:
      return ml2(2,0); 
   case CNE6SSM_info::ml221:
      return ml2(2,1);
   case CNE6SSM_info::ml222:
      return ml2(2,2);
   case CNE6SSM_info::mHd2:
      return mHd2; 
   case CNE6SSM_info::mHu2:
      return mHu2; 
   case CNE6SSM_info::md200:
      return md2(0,0); 
   case CNE6SSM_info::md201:
      return md2(0,1); 
   case CNE6SSM_info::md202:
      return md2(0,2); 
   case CNE6SSM_info::md210:
      return md2(1,0); 
   case CNE6SSM_info::md211:
      return md2(1,1); 
   case CNE6SSM_info::md212:
      return md2(1,2);
   case CNE6SSM_info::md220:
      return md2(2,0); 
   case CNE6SSM_info::md221:
      return md2(2,1);
   case CNE6SSM_info::md222:
      return md2(2,2); 
   case CNE6SSM_info::mu200:
      return mu2(0,0); 
   case CNE6SSM_info::mu201:
      return mu2(0,1); 
   case CNE6SSM_info::mu202:
      return mu2(0,2); 
   case CNE6SSM_info::mu210:
      return mu2(1,0); 
   case CNE6SSM_info::mu211:
      return mu2(1,1); 
   case CNE6SSM_info::mu212:
      return mu2(1,2); 
   case CNE6SSM_info::mu220:
      return mu2(2,0); 
   case CNE6SSM_info::mu221:
      return mu2(2,1); 
   case CNE6SSM_info::mu222:
      return mu2(2,2); 
   case CNE6SSM_info::me200:
      return me2(0,0);
   case CNE6SSM_info::me201:
      return me2(0,1); 
   case CNE6SSM_info::me202:
      return me2(0,2); 
   case CNE6SSM_info::me210:
      return me2(1,0); 
   case CNE6SSM_info::me211:
      return me2(1,1); 
   case CNE6SSM_info::me212:
      return me2(1,2); 
   case CNE6SSM_info::me220:
      return me2(2,0);
   case CNE6SSM_info::me221:
      return me2(2,1); 
   case CNE6SSM_info::me222:
      return me2(2,2); 
   case CNE6SSM_info::ms2:
      return ms2; 
   case CNE6SSM_info::msbar2:
      return msbar2;
   case CNE6SSM_info::mH1I200:
      return mH1I2(0,0); 
   case CNE6SSM_info::mH1I201:
      return mH1I2(0,1); 
   case CNE6SSM_info::mH1I210:
      return mH1I2(1,0); 
   case CNE6SSM_info::mH1I211:
      return mH1I2(1,1); 
   case CNE6SSM_info::mH2I200:
      return mH2I2(0,0); 
   case CNE6SSM_info::mH2I201:
      return mH2I2(0,1); 
   case CNE6SSM_info::mH2I210:
      return mH2I2(1,0); 
   case CNE6SSM_info::mH2I211:
      return mH2I2(1,1);
   case CNE6SSM_info::mSI200:
      return mSI2(0,0); 
   case CNE6SSM_info::mSI201:
      return mSI2(0,1); 
   case CNE6SSM_info::mSI202:
      return mSI2(0,2); 
   case CNE6SSM_info::mSI210:
      return mSI2(1,0); 
   case CNE6SSM_info::mSI211:
      return mSI2(1,1); 
   case CNE6SSM_info::mSI212:
      return mSI2(1,2); 
   case CNE6SSM_info::mSI220:
      return mSI2(2,0);
   case CNE6SSM_info::mSI221:
      return mSI2(2,1); 
   case CNE6SSM_info::mSI222:
      return mSI2(2,2);
   case CNE6SSM_info::mDx200:
      return mDx2(0,0); 
   case CNE6SSM_info::mDx201:
      return mDx2(0,1); 
   case CNE6SSM_info::mDx202:
      return mDx2(0,2); 
   case CNE6SSM_info::mDx210:
      return mDx2(1,0); 
   case CNE6SSM_info::mDx211:
      return mDx2(1,1); 
   case CNE6SSM_info::mDx212:
      return mDx2(1,2); 
   case CNE6SSM_info::mDx220:
      return mDx2(2,0); 
   case CNE6SSM_info::mDx221:
      return mDx2(2,1); 
   case CNE6SSM_info::mDx222:
      return mDx2(2,2);
   case CNE6SSM_info::mDxbar200:
      return mDxbar2(0,0); 
   case CNE6SSM_info::mDxbar201:
      return mDxbar2(0,1); 
   case CNE6SSM_info::mDxbar202:
      return mDxbar2(0,2); 
   case CNE6SSM_info::mDxbar210:
      return mDxbar2(1,0); 
   case CNE6SSM_info::mDxbar211:
      return mDxbar2(1,1); 
   case CNE6SSM_info::mDxbar212:
      return mDxbar2(1,2); 
   case CNE6SSM_info::mDxbar220:
      return mDxbar2(2,0);
   case CNE6SSM_info::mDxbar221:
      return mDxbar2(2,1);
   case CNE6SSM_info::mDxbar222:
      return mDxbar2(2,2); 
   case CNE6SSM_info::mHp2:
      return mHp2; 
   case CNE6SSM_info::mHpbar2:
      return mHpbar2; 
   case CNE6SSM_info::mphi2:
      return mphi2; 
   case CNE6SSM_info::MassB:
      return MassB; 
   case CNE6SSM_info::MassWB:
      return MassWB; 
   case CNE6SSM_info::MassG:
      return MassG; 
   case CNE6SSM_info::MassBp:
      return MassBp;

   default:
      throw UnknownModelParameterError(parameter);
   }
}

void CLASSNAME::set_parameter(unsigned parameter, double x)
{
   if (parameter >= CNE6SSM_info::NUMBER_OF_PARAMETERS)
      throw UnknownModelParameterError(parameter);

   switch (parameter) {

   case CNE6SSM_info::Yd00:
      Yd(0,0) = x;
      break;
   case CNE6SSM_info::Yd01:
      Yd(0,1) = x;
      break;
   case CNE6SSM_info::Yd02:
      Yd(0,2) = x;
      break;
   case CNE6SSM_info::Yd10:
      Yd(1,0) = x;
      break;
   case CNE6SSM_info::Yd11:
      Yd(1,1) = x;
      break;
   case CNE6SSM_info::Yd12:
      Yd(1,2) = x;
      break;
   case CNE6SSM_info::Yd20:
      Yd(2,0) = x;
      break;
   case CNE6SSM_info::Yd21:
      Yd(2,1) = x;
      break;
   case CNE6SSM_info::Yd22:
      Yd(2,2) = x;
      break;
   case CNE6SSM_info::hE00:
      hE(0,0) = x;
      break;
   case CNE6SSM_info::hE01:
      hE(0,1) = x;
      break;
   case CNE6SSM_info::hE10:
      hE(1,0) = x;
      break;
   case CNE6SSM_info::hE11:
      hE(1,1) = x;
      break;
   case CNE6SSM_info::hE20:
      hE(2,0) = x;
      break;
   case CNE6SSM_info::hE21:
      hE(2,1) = x;
      break;
   case CNE6SSM_info::Ye00:
      Ye(0,0) = x; 
      break;
   case CNE6SSM_info::Ye01:
      Ye(0,1) = x;
      break;
   case CNE6SSM_info::Ye02:
      Ye(0,2) = x;
      break;
   case CNE6SSM_info::Ye10:
      Ye(1,0) = x;
      break;
   case CNE6SSM_info::Ye11:
      Ye(1,1) = x;
      break;
   case CNE6SSM_info::Ye12:
      Ye(1,2) = x;
      break;
   case CNE6SSM_info::Ye20:
      Ye(2,0) = x;
      break;
   case CNE6SSM_info::Ye21:
      Ye(2,1) = x;
      break;
   case CNE6SSM_info::Ye22:
      Ye(2,2) = x;
      break;
   case CNE6SSM_info::SigmaL:
      SigmaL = x;
      break;
   case CNE6SSM_info::KappaPr:
      KappaPr = x;
      break;
   case CNE6SSM_info::Sigmax:
      Sigmax = x;
      break;
   case CNE6SSM_info::gD00:
      gD(0,0) = x;
      break;
   case CNE6SSM_info::gD01:
      gD(0,1) = x;
      break;
   case CNE6SSM_info::gD02:
      gD(0,2) = x;
      break;
   case CNE6SSM_info::gD10:
      gD(1,0) = x;
      break;
   case CNE6SSM_info::gD11:
      gD(1,1) = x;
      break;
   case CNE6SSM_info::gD12:
      gD(1,2) = x;
      break;
   case CNE6SSM_info::gD20:
      gD(2,0) = x;
      break;
   case CNE6SSM_info::gD21:
      gD(2,1) = x;
      break;
   case CNE6SSM_info::gD22:
      gD(2,2) = x;
      break;
   case CNE6SSM_info::Kappa00:
      Kappa(0,0) = x;
      break;
   case CNE6SSM_info::Kappa01:
      Kappa(0,1) = x;
      break;
   case CNE6SSM_info::Kappa02:
      Kappa(0,2) = x;
      break;
   case CNE6SSM_info::Kappa10:
      Kappa(1,0) = x;
      break;
   case CNE6SSM_info::Kappa11:
      Kappa(1,1) = x;
      break;
   case CNE6SSM_info::Kappa12:
      Kappa(1,2) = x;
      break;
   case CNE6SSM_info::Kappa20:
      Kappa(2,0) = x;
      break;
   case CNE6SSM_info::Kappa21:
      Kappa(2,1) = x;
      break;
   case CNE6SSM_info::Kappa22:
      Kappa(2,2) = x;
      break;
   case CNE6SSM_info::Lambda1200:
      Lambda12(0,0) = x;
      break;
   case CNE6SSM_info::Lambda1201:
      Lambda12(0,1) = x;
      break;
   case CNE6SSM_info::Lambda1210:
      Lambda12(1,0) = x;
      break;
   case CNE6SSM_info::Lambda1211:
      Lambda12(1,1) = x;
      break;
   case CNE6SSM_info::Lambdax:
      Lambdax = x;
      break;
   case CNE6SSM_info::fu00:
      fu(0,0) = x;
      break;
   case CNE6SSM_info::fu01:
      fu(0,1) = x;
      break;
   case CNE6SSM_info::fu10:
      fu(1,0) = x;
      break;
   case CNE6SSM_info::fu11:
      fu(1,1) = x;
      break;
   case CNE6SSM_info::fu20:
      fu(2,0) = x;
      break;
   case CNE6SSM_info::fu21:
      fu(2,1) = x;
      break;
   case CNE6SSM_info::fd00:
      fd(0,0) = x;
      break;
   case CNE6SSM_info::fd01:
      fd(0,1) = x;
      break;
   case CNE6SSM_info::fd10:
      fd(1,0) = x;
      break;
   case CNE6SSM_info::fd11:
      fd(1,1) = x;
      break;
   case CNE6SSM_info::fd20:
      fd(2,0) = x;
      break;
   case CNE6SSM_info::fd21:
      fd(2,1) = x;
      break;
   case CNE6SSM_info::Yu00:
      Yu(0,0) = x;
      break;
   case CNE6SSM_info::Yu01:
      Yu(0,1) = x;
      break;
   case CNE6SSM_info::Yu02:
      Yu(0,2) = x;
      break;
   case CNE6SSM_info::Yu10:
      Yu(1,0) = x;
      break;
   case CNE6SSM_info::Yu11:
      Yu(1,1) = x;
      break;
   case CNE6SSM_info::Yu12:
      Yu(1,2) = x;
      break;
   case CNE6SSM_info::Yu20:
      Yu(2,0) = x;
      break;
   case CNE6SSM_info::Yu21:
      Yu(2,1) = x;
      break;
   case CNE6SSM_info::Yu22:
      Yu(2,2) = x;
      break;
   case CNE6SSM_info::MuPr:
      MuPr = x;
      break;
   case CNE6SSM_info::MuPhi:
      MuPhi = x;
      break;
   case CNE6SSM_info::XiF:
      XiF = x;
      break;
   case CNE6SSM_info::g1:
      g1 = x;
      break;
   case CNE6SSM_info::g2:
      g2 = x;
      break;
   case CNE6SSM_info::g3:
      g3 = x;
      break;
   case CNE6SSM_info::g1p:
      g1p = x;
      break;
   case CNE6SSM_info::vd:
      vd = x;
      break;
   case CNE6SSM_info::vu:
      vu = x;
      break;
   case CNE6SSM_info::vs:
      vs = x;
      break;
   case CNE6SSM_info::vsb:
      vsb = x;
      break;
   case CNE6SSM_info::vphi:
      vphi = x;
      break;
   case CNE6SSM_info::TYd00:
      TYd(0,0) = x;
      break;
   case CNE6SSM_info::TYd01:
      TYd(0,1) = x;
      break;
   case CNE6SSM_info::TYd02:
      TYd(0,2) = x;
      break;
   case CNE6SSM_info::TYd10:
      TYd(1,0) = x;
      break;
   case CNE6SSM_info::TYd11:
      TYd(1,1) = x;
      break;
   case CNE6SSM_info::TYd12:
      TYd(1,2) = x;
      break; 
   case CNE6SSM_info::TYd20:
      TYd(2,0) = x;
      break;
   case CNE6SSM_info::TYd21:
      TYd(2,1) = x;
      break;
   case CNE6SSM_info::TYd22:
      TYd(2,2) = x;
      break;
   case CNE6SSM_info::ThE00:
      ThE(0,0) = x;
      break;
   case CNE6SSM_info::ThE01:
      ThE(0,1) = x;
      break;
   case CNE6SSM_info::ThE10:
      ThE(1,0) = x;
      break;
   case CNE6SSM_info::ThE11:
      ThE(1,1) = x;
      break;
   case CNE6SSM_info::ThE20:
      ThE(2,0) = x;
      break;
   case CNE6SSM_info::ThE21:
      ThE(2,1) = x;
      break;
   case CNE6SSM_info::TYe00:
      TYe(0,0) = x;
      break;
   case CNE6SSM_info::TYe01:
      TYe(0,1) = x;
      break;
   case CNE6SSM_info::TYe02:
      TYe(0,2) = x;
      break;
   case CNE6SSM_info::TYe10:
      TYe(1,0) = x;
      break;
   case CNE6SSM_info::TYe11:
      TYe(1,1) = x;
      break;
   case CNE6SSM_info::TYe12:
      TYe(1,2) = x;
      break;
   case CNE6SSM_info::TYe20:
      TYe(2,0) = x;
      break;
   case CNE6SSM_info::TYe21:
      TYe(2,1) = x;
      break;
   case CNE6SSM_info::TYe22:
      TYe(2,2) = x;
      break;
   case CNE6SSM_info::TSigmaL:
      TSigmaL = x;
      break;
   case CNE6SSM_info::TKappaPr:
      TKappaPr = x;
      break;
   case CNE6SSM_info::TSigmax:
      TSigmax = x;
      break;
   case CNE6SSM_info::TgD00:
      TgD(0,0) = x;
      break;
   case CNE6SSM_info::TgD01:
      TgD(0,1) = x;
      break;
   case CNE6SSM_info::TgD02:
      TgD(0,2) = x;
      break;
   case CNE6SSM_info::TgD10:
      TgD(1,0) = x;
      break;
   case CNE6SSM_info::TgD11:
      TgD(1,1) = x;
      break;
   case CNE6SSM_info::TgD12:
      TgD(1,2) = x;
      break;
   case CNE6SSM_info::TgD20:
      TgD(2,0) = x;
      break;
   case CNE6SSM_info::TgD21:
      TgD(2,1) = x;
      break;
   case CNE6SSM_info::TgD22:
      TgD(2,2) = x;
      break;
   case CNE6SSM_info::TKappa00:
      TKappa(0,0) = x;
      break;
   case CNE6SSM_info::TKappa01:
      TKappa(0,1) = x;
      break;
   case CNE6SSM_info::TKappa02:
      TKappa(0,2) = x;
      break;
   case CNE6SSM_info::TKappa10:
      TKappa(1,0) = x;
      break;
   case CNE6SSM_info::TKappa11:
      TKappa(1,1) = x;
      break;
   case CNE6SSM_info::TKappa12:
      TKappa(1,2) = x;
      break;
   case CNE6SSM_info::TKappa20:
      TKappa(2,0) = x;
      break;
   case CNE6SSM_info::TKappa21:
      TKappa(2,1) = x;
      break;
   case CNE6SSM_info::TKappa22:
      TKappa(2,2) = x;
      break;
   case CNE6SSM_info::TLambda1200:
      TLambda12(0,0) = x;
      break;
   case CNE6SSM_info::TLambda1201:
      TLambda12(0,1) = x;
      break;
   case CNE6SSM_info::TLambda1210:
      TLambda12(1,0) = x;
      break;
   case CNE6SSM_info::TLambda1211:
      TLambda12(1,1) = x;
      break;
   case CNE6SSM_info::TLambdax:
      TLambdax = x;
      break;
   case CNE6SSM_info::Tfu00:
      Tfu(0,0) = x;
      break;
   case CNE6SSM_info::Tfu01:
      Tfu(0,1) = x;
      break;
   case CNE6SSM_info::Tfu10:
      Tfu(1,0) = x;
      break;
   case CNE6SSM_info::Tfu11:
      Tfu(1,1) = x;
      break;
   case CNE6SSM_info::Tfu20:
      Tfu(2,0) = x;
      break;
   case CNE6SSM_info::Tfu21:
      Tfu(2,1) = x;
      break;
   case CNE6SSM_info::Tfd00:
      Tfd(0,0) = x;
      break;
   case CNE6SSM_info::Tfd01:
      Tfd(0,1) = x;
      break;
   case CNE6SSM_info::Tfd10:
      Tfd(1,0) = x;
      break;
   case CNE6SSM_info::Tfd11:
      Tfd(1,1) = x;
      break;
   case CNE6SSM_info::Tfd20:
      Tfd(2,0) = x;
      break;
   case CNE6SSM_info::Tfd21:
      Tfd(2,1) = x;
      break;
   case CNE6SSM_info::TYu00:
      TYu(0,0) = x;
      break;
   case CNE6SSM_info::TYu01:
      TYu(0,1) = x;
      break;
   case CNE6SSM_info::TYu02:
      TYu(0,2) = x;
      break;
   case CNE6SSM_info::TYu10:
      TYu(1,0) = x;
      break;
   case CNE6SSM_info::TYu11:
      TYu(1,1) = x;
      break;
   case CNE6SSM_info::TYu12:
      TYu(1,2) = x;
      break;
   case CNE6SSM_info::TYu20:
      TYu(2,0) = x;
      break;
   case CNE6SSM_info::TYu21:
      TYu(2,1) = x;
      break;
   case CNE6SSM_info::TYu22:
      TYu(2,2) = x;
      break;
   case CNE6SSM_info::BMuPr:
      BMuPr = x;
      break;
   case CNE6SSM_info::BMuPhi:
      BMuPhi = x;
      break;
   case CNE6SSM_info::LXiF:
      LXiF = x;
      break;
   case CNE6SSM_info::mq200:
      mq2(0,0) = x;
      break;
   case CNE6SSM_info::mq201:
      mq2(0,1) = x;
      break;
   case CNE6SSM_info::mq202:
      mq2(0,2) = x;
      break;
   case CNE6SSM_info::mq210:
      mq2(1,0) = x;
      break;
   case CNE6SSM_info::mq211:
      mq2(1,1) = x;
      break;
   case CNE6SSM_info::mq212:
      mq2(1,2) = x;
      break;
   case CNE6SSM_info::mq220:
      mq2(2,0) = x;
      break;
   case CNE6SSM_info::mq221:
      mq2(2,1) = x;
      break;
   case CNE6SSM_info::mq222:
      mq2(2,2) = x;
      break;
   case CNE6SSM_info::ml200:
      ml2(0,0) = x;
      break;
   case CNE6SSM_info::ml201:
      ml2(0,1) = x;
      break;
   case CNE6SSM_info::ml202:
      ml2(0,2) = x;
      break;
   case CNE6SSM_info::ml210:
      ml2(1,0) = x;
      break;
   case CNE6SSM_info::ml211:
      ml2(1,1) = x;
      break;
   case CNE6SSM_info::ml212:
      ml2(1,2) = x;
      break;
   case CNE6SSM_info::ml220:
      ml2(2,0) = x;
      break;
   case CNE6SSM_info::ml221:
      ml2(2,1) = x;
      break;
   case CNE6SSM_info::ml222:
      ml2(2,2) = x;
      break;
   case CNE6SSM_info::mHd2:
      mHd2 = x;
      break;
   case CNE6SSM_info::mHu2:
      mHu2 = x;
      break;
   case CNE6SSM_info::md200:
      md2(0,0) = x;
      break;
   case CNE6SSM_info::md201:
      md2(0,1) = x;
      break;
   case CNE6SSM_info::md202:
      md2(0,2) = x;
      break;
   case CNE6SSM_info::md210:
      md2(1,0) = x;
      break;
   case CNE6SSM_info::md211:
      md2(1,1) = x;
      break;
   case CNE6SSM_info::md212:
      md2(1,2) = x;
      break;
   case CNE6SSM_info::md220:
      md2(2,0) = x;
      break;
   case CNE6SSM_info::md221:
      md2(2,1) = x;
      break;
   case CNE6SSM_info::md222:
      md2(2,2) = x;
      break;
   case CNE6SSM_info::mu200:
      mu2(0,0) = x;
      break;
   case CNE6SSM_info::mu201:
      mu2(0,1) = x;
      break;
   case CNE6SSM_info::mu202:
      mu2(0,2) = x;
      break;
   case CNE6SSM_info::mu210:
      mu2(1,0) = x;
      break;
   case CNE6SSM_info::mu211:
      mu2(1,1) = x;
      break;
   case CNE6SSM_info::mu212:
      mu2(1,2) = x;
      break;
   case CNE6SSM_info::mu220:
      mu2(2,0) = x;
      break;
   case CNE6SSM_info::mu221:
      mu2(2,1) = x;
      break;
   case CNE6SSM_info::mu222:
      mu2(2,2) = x;
      break;
   case CNE6SSM_info::me200:
      me2(0,0) = x;
      break;
   case CNE6SSM_info::me201:
      me2(0,1) = x;
      break;
   case CNE6SSM_info::me202:
      me2(0,2) = x;
      break;
   case CNE6SSM_info::me210:
      me2(1,0) = x;
      break;
   case CNE6SSM_info::me211:
      me2(1,1) = x;
      break;
   case CNE6SSM_info::me212:
      me2(1,2) = x;
      break;
   case CNE6SSM_info::me220:
      me2(2,0) = x;
      break;
   case CNE6SSM_info::me221:
      me2(2,1) = x;
      break;
   case CNE6SSM_info::me222:
      me2(2,2) = x;
      break;
   case CNE6SSM_info::ms2:
      ms2 = x;
      break;
   case CNE6SSM_info::msbar2:
      msbar2 = x;
      break;
   case CNE6SSM_info::mH1I200:
      mH1I2(0,0) = x;
      break;
   case CNE6SSM_info::mH1I201:
      mH1I2(0,1) = x;
      break;
   case CNE6SSM_info::mH1I210:
      mH1I2(1,0) = x;
      break;
   case CNE6SSM_info::mH1I211:
      mH1I2(1,1) = x;
      break;
   case CNE6SSM_info::mH2I200:
      mH2I2(0,0) = x;
      break;
   case CNE6SSM_info::mH2I201:
      mH2I2(0,1) = x;
      break;
   case CNE6SSM_info::mH2I210:
      mH2I2(1,0) = x;
      break;
   case CNE6SSM_info::mH2I211:
      mH2I2(1,1) = x;
      break;
   case CNE6SSM_info::mSI200:
      mSI2(0,0) = x;
      break;
   case CNE6SSM_info::mSI201:
      mSI2(0,1) = x;
      break;
   case CNE6SSM_info::mSI202:
      mSI2(0,2) = x;
      break;
   case CNE6SSM_info::mSI210:
      mSI2(1,0) = x;
      break;
   case CNE6SSM_info::mSI211:
      mSI2(1,1) = x;
      break;
   case CNE6SSM_info::mSI212:
      mSI2(1,2) = x;
      break;
   case CNE6SSM_info::mSI220:
      mSI2(2,0) = x;
      break;
   case CNE6SSM_info::mSI221:
      mSI2(2,1) = x;
      break;
   case CNE6SSM_info::mSI222:
      mSI2(2,2) = x;
      break;
   case CNE6SSM_info::mDx200:
      mDx2(0,0) = x;
      break;
   case CNE6SSM_info::mDx201:
      mDx2(0,1) = x;
      break;
   case CNE6SSM_info::mDx202:
      mDx2(0,2) = x;
      break;
   case CNE6SSM_info::mDx210:
      mDx2(1,0) = x;
      break;
   case CNE6SSM_info::mDx211:
      mDx2(1,1) = x;
      break;
   case CNE6SSM_info::mDx212:
      mDx2(1,2) = x;
      break;
   case CNE6SSM_info::mDx220:
      mDx2(2,0) = x;
      break;
   case CNE6SSM_info::mDx221:
      mDx2(2,1) = x;
      break;
   case CNE6SSM_info::mDx222:
      mDx2(2,2) = x;
      break;
   case CNE6SSM_info::mDxbar200:
      mDxbar2(0,0) = x;
      break;
   case CNE6SSM_info::mDxbar201:
      mDxbar2(0,1) = x;
      break;
   case CNE6SSM_info::mDxbar202:
      mDxbar2(0,2) = x;
      break;
   case CNE6SSM_info::mDxbar210:
      mDxbar2(1,0) = x;
      break;
   case CNE6SSM_info::mDxbar211:
      mDxbar2(1,1) = x;
      break;
   case CNE6SSM_info::mDxbar212:
      mDxbar2(1,2) = x;
      break;
   case CNE6SSM_info::mDxbar220:
      mDxbar2(2,0) = x;
      break;
   case CNE6SSM_info::mDxbar221:
      mDxbar2(2,1) = x;
      break;
   case CNE6SSM_info::mDxbar222:
      mDxbar2(2,2) = x;
      break;
   case CNE6SSM_info::mHp2:
      mHp2 = x;
      break;
   case CNE6SSM_info::mHpbar2:
      mHpbar2 = x;
      break;
   case CNE6SSM_info::mphi2:
      mphi2 = x;
      break;
   case CNE6SSM_info::MassB:
      MassB = x;
      break;
   case CNE6SSM_info::MassWB:
      MassWB = x;
      break;
   case CNE6SSM_info::MassG:
      MassG = x;
      break;
   case CNE6SSM_info::MassBp:
      MassBp = x;
      break;

   default:
      throw UnknownModelParameterError(parameter);
   }
}

void CLASSNAME::calculate_coefficients(double input_scale)
{
   static const double fit_Azero_values[number_of_fit_points] = {0., 1., 0., 1.};
   static const double fit_m12_values[number_of_fit_points] = {1., 0., 0., 1.};
   static const double fit_m0_values[number_of_fit_points] = {0., 0., 1., 0.};

   // save current set of parameters
   CNE6SSM_soft_parameters saved_pars;
   saved_pars.set(get());

   const double current_scale = get_scale();

   run_to(input_scale, precision);

   CNE6SSM_soft_parameters input_scale_pars;
   input_scale_pars.set(get());

   std::vector<CNE6SSM_soft_parameters> parameter_values;

   Eigen::Matrix<double,number_of_fit_points,2> dimension_one_inputs;
   Eigen::Matrix<double,number_of_fit_points,4> dimension_two_inputs;

   for (std::size_t i = 0; i < number_of_fit_points; ++i) {
      dimension_one_inputs(i,0) = fit_Azero_values[i];
      dimension_one_inputs(i,1) = fit_m12_values[i];

      dimension_two_inputs(i,0) = Sqr(fit_m0_values[i]);
      dimension_two_inputs(i,1) = Sqr(fit_m12_values[i]);
      dimension_two_inputs(i,2) = fit_m12_values[i] * fit_Azero_values[i];
      dimension_two_inputs(i,3) = Sqr(fit_Azero_values[i]);

      set_soft_parameters_at_input_scale(fit_m0_values[i], fit_m12_values[i], fit_Azero_values[i]);

      run_to(current_scale, precision);

      CNE6SSM_soft_parameters params;
      params.set(get());

      parameter_values.push_back(params);

      set(input_scale_pars.get());
      set_scale(input_scale);
   }

   // solve for coefficients using least squares
   Eigen::JacobiSVD<Eigen::Matrix<double,number_of_fit_points,2> > dimension_one_svd(dimension_one_inputs);
   Eigen::JacobiSVD<Eigen::Matrix<double,number_of_fit_points,4> > dimension_two_svd(dimension_two_inputs);

   Eigen::Matrix<double,number_of_fit_points,1> rhs;
   Eigen::Matrix<double,2,1> dimension_one_solution;
   Eigen::Matrix<double,4,1> dimension_two_solution;

   // TODO temporary working (awful!) solution, needs to be replaced...
   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_TYd(i,j);
         }
         dimension_one_solution = dimension_one_svd.solve(rhs);
         TYd_Azero_coeff(i,j) = dimension_one_solution(0);
         TYd_m12_coeff(i,j) = dimension_one_solution(1);
      }
   }
   
   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 2; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_ThE(i,j);
         }
         dimension_one_solution = dimension_one_svd.solve(rhs);
         ThE_Azero_coeff(i,j) = dimension_one_solution(0);
         ThE_m12_coeff(i,j) = dimension_one_solution(1);
      }
   }
   
   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_TYe(i,j);
         }
         dimension_one_solution = dimension_one_svd.solve(rhs);
         TYe_Azero_coeff(i,j) = dimension_one_solution(0);
         TYe_m12_coeff(i,j) = dimension_one_solution(1);
      }
   }
   
   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_TSigmaL();
   }
   dimension_one_solution = dimension_one_svd.solve(rhs);
   TSigmaL_Azero_coeff = dimension_one_solution(0);
   TSigmaL_m12_coeff = dimension_one_solution(1);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_TKappaPr();
   }
   dimension_one_solution = dimension_one_svd.solve(rhs);
   TKappaPr_Azero_coeff = dimension_one_solution(0);
   TKappaPr_m12_coeff = dimension_one_solution(1);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_TSigmax();
   }
   dimension_one_solution = dimension_one_svd.solve(rhs);
   TSigmax_Azero_coeff = dimension_one_solution(0);
   TSigmax_m12_coeff = dimension_one_solution(1);

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_TgD(i,j);
         }
         dimension_one_solution = dimension_one_svd.solve(rhs);
         TgD_Azero_coeff(i,j) = dimension_one_solution(0);
         TgD_m12_coeff(i,j) = dimension_one_solution(1);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_TKappa(i,j);
         }
         dimension_one_solution = dimension_one_svd.solve(rhs);
         TKappa_Azero_coeff(i,j) = dimension_one_solution(0);
         TKappa_m12_coeff(i,j) = dimension_one_solution(1);
      }
   }

   for (std::size_t i = 0; i < 2; ++i) {
      for (std::size_t j = 0; j < 2; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_TLambda12(i,j);
         }
         dimension_one_solution = dimension_one_svd.solve(rhs);
         TLambda12_Azero_coeff(i,j) = dimension_one_solution(0);
         TLambda12_m12_coeff(i,j) = dimension_one_solution(1);
      }
   }

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_TLambdax();
   }
   dimension_one_solution = dimension_one_svd.solve(rhs);
   TLambdax_Azero_coeff = dimension_one_solution(0);
   TLambdax_m12_coeff = dimension_one_solution(1);

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 2; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_Tfu(i,j);
         }
         dimension_one_solution = dimension_one_svd.solve(rhs);
         Tfu_Azero_coeff(i,j) = dimension_one_solution(0);
         Tfu_m12_coeff(i,j) = dimension_one_solution(1);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 2; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_Tfd(i,j);
         }
         dimension_one_solution = dimension_one_svd.solve(rhs);
         Tfd_Azero_coeff(i,j) = dimension_one_solution(0);
         Tfd_m12_coeff(i,j) = dimension_one_solution(1);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_TYu(i,j);
         }
         dimension_one_solution = dimension_one_svd.solve(rhs);
         TYu_Azero_coeff(i,j) = dimension_one_solution(0);
         TYu_m12_coeff(i,j) = dimension_one_solution(1);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_mq2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         mq2_m02_coeff(i,j) = dimension_two_solution(0);
         mq2_m122_coeff(i,j) = dimension_two_solution(1);
         mq2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         mq2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_ml2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         ml2_m02_coeff(i,j) = dimension_two_solution(0);
         ml2_m122_coeff(i,j) = dimension_two_solution(1);
         ml2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         ml2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_mHd2();
   }
   dimension_two_solution = dimension_two_svd.solve(rhs);
   mHd2_m02_coeff = dimension_two_solution(0);
   mHd2_m122_coeff = dimension_two_solution(1);
   mHd2_Azerom12_coeff = dimension_two_solution(2);
   mHd2_Azero2_coeff = dimension_two_solution(3);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_mHu2();
   }
   dimension_two_solution = dimension_two_svd.solve(rhs);
   mHu2_m02_coeff = dimension_two_solution(0);
   mHu2_m122_coeff = dimension_two_solution(1);
   mHu2_Azerom12_coeff = dimension_two_solution(2);
   mHu2_Azero2_coeff = dimension_two_solution(3);

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_md2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         md2_m02_coeff(i,j) = dimension_two_solution(0);
         md2_m122_coeff(i,j) = dimension_two_solution(1);
         md2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         md2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_mu2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         mu2_m02_coeff(i,j) = dimension_two_solution(0);
         mu2_m122_coeff(i,j) = dimension_two_solution(1);
         mu2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         mu2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_me2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         me2_m02_coeff(i,j) = dimension_two_solution(0);
         me2_m122_coeff(i,j) = dimension_two_solution(1);
         me2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         me2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_ms2();
   }
   dimension_two_solution = dimension_two_svd.solve(rhs);
   ms2_m02_coeff = dimension_two_solution(0);
   ms2_m122_coeff = dimension_two_solution(1);
   ms2_Azerom12_coeff = dimension_two_solution(2);
   ms2_Azero2_coeff = dimension_two_solution(3);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_msbar2();
   }
   dimension_two_solution = dimension_two_svd.solve(rhs);
   msbar2_m02_coeff = dimension_two_solution(0);
   msbar2_m122_coeff = dimension_two_solution(1);
   msbar2_Azerom12_coeff = dimension_two_solution(2);
   msbar2_Azero2_coeff = dimension_two_solution(3);

   for (std::size_t i = 0; i < 2; ++i) {
      for (std::size_t j = 0; j < 2; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_mH1I2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         mH1I2_m02_coeff(i,j) = dimension_two_solution(0);
         mH1I2_m122_coeff(i,j) = dimension_two_solution(1);
         mH1I2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         mH1I2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t i = 0; i < 2; ++i) {
      for (std::size_t j = 0; j < 2; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_mH2I2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         mH2I2_m02_coeff(i,j) = dimension_two_solution(0);
         mH2I2_m122_coeff(i,j) = dimension_two_solution(1);
         mH2I2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         mH2I2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_mSI2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         mSI2_m02_coeff(i,j) = dimension_two_solution(0);
         mSI2_m122_coeff(i,j) = dimension_two_solution(1);
         mSI2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         mSI2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_mDx2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         mDx2_m02_coeff(i,j) = dimension_two_solution(0);
         mDx2_m122_coeff(i,j) = dimension_two_solution(1);
         mDx2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         mDx2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
         for (std::size_t k = 0; k < number_of_fit_points; ++k) {
            rhs(k) = parameter_values[k].get_mDxbar2(i,j);
         }
         dimension_two_solution = dimension_two_svd.solve(rhs);
         mDxbar2_m02_coeff(i,j) = dimension_two_solution(0);
         mDxbar2_m122_coeff(i,j) = dimension_two_solution(1);
         mDxbar2_Azerom12_coeff(i,j) = dimension_two_solution(2);
         mDxbar2_Azero2_coeff(i,j) = dimension_two_solution(3);
      }
   }

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_mHp2();
   }
   dimension_two_solution = dimension_two_svd.solve(rhs);
   mHp2_m02_coeff = dimension_two_solution(0);
   mHp2_m122_coeff = dimension_two_solution(1);
   mHp2_Azerom12_coeff = dimension_two_solution(2);
   mHp2_Azero2_coeff = dimension_two_solution(3);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_mHpbar2();
   }
   dimension_two_solution = dimension_two_svd.solve(rhs);
   mHpbar2_m02_coeff = dimension_two_solution(0);
   mHpbar2_m122_coeff = dimension_two_solution(1);
   mHpbar2_Azerom12_coeff = dimension_two_solution(2);
   mHpbar2_Azero2_coeff = dimension_two_solution(3);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_mphi2();
   }
   dimension_two_solution = dimension_two_svd.solve(rhs);
   mphi2_m02_coeff = dimension_two_solution(0);
   mphi2_m122_coeff = dimension_two_solution(1);
   mphi2_Azerom12_coeff = dimension_two_solution(2);
   mphi2_Azero2_coeff = dimension_two_solution(3);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_MassB();
   }
   dimension_one_solution = dimension_one_svd.solve(rhs);
   MassB_Azero_coeff = dimension_one_solution(0);
   MassB_m12_coeff = dimension_one_solution(1);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_MassWB();
   }
   dimension_one_solution = dimension_one_svd.solve(rhs);
   MassWB_Azero_coeff = dimension_one_solution(0);
   MassWB_m12_coeff = dimension_one_solution(1);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_MassG();
   }
   dimension_one_solution = dimension_one_svd.solve(rhs);
   MassG_Azero_coeff = dimension_one_solution(0);
   MassG_m12_coeff = dimension_one_solution(1);

   for (std::size_t k = 0; k < number_of_fit_points; ++k) {
      rhs(k) = parameter_values[k].get_MassBp();
   }
   dimension_one_solution = dimension_one_svd.solve(rhs);
   MassBp_Azero_coeff = dimension_one_solution(0);
   MassBp_m12_coeff = dimension_one_solution(1);

   // reset parameters at initial scale
   set(saved_pars.get());
   set_scale(current_scale);
}

void CLASSNAME::set_soft_parameters_at_input_scale(double m0, double m12, double Azero)
{
   TYd = Azero * Yd;
   ThE = Azero * hE;
   TYe = Azero * Ye;
   TSigmaL = Azero * SigmaL;
   TKappaPr = Azero * KappaPr;
   TSigmax = Azero * Sigmax;
   TgD = Azero * gD;
   TKappa = Azero * Kappa;
   TLambda12 = Azero * Lambda12;
   TLambdax = Azero * Lambdax;
   Tfu = Azero * fu;
   Tfd = Azero * fd;
   mq2 = Sqr(m0) * UNITMATRIX(3);
   ml2 = Sqr(m0) * UNITMATRIX(3);
   mHd2 = Sqr(m0);
   mHu2 = Sqr(m0);
   md2 = Sqr(m0) * UNITMATRIX(3);
   mu2 = Sqr(m0) * UNITMATRIX(3);
   me2 = Sqr(m0) * UNITMATRIX(3);
   ms2 = Sqr(m0);
   msbar2 = Sqr(m0);
   mH1I2 = Sqr(m0) * UNITMATRIX(2);
   mH2I2 = Sqr(m0) * UNITMATRIX(2);
   mSI2 = Sqr(m0) * UNITMATRIX(3);
   mDx2 = Sqr(m0) * UNITMATRIX(3);
   mDxbar2 = Sqr(m0) * UNITMATRIX(3);
   mHp2 = Sqr(m0);
   mHpbar2 = Sqr(m0);
   mphi2 = Sqr(m0);
   MassB = m12;
   MassWB = m12;
   MassG = m12;
   MassBp = m12;
}

void CLASSNAME::set_soft_parameters_at_current_scale(double m0, double m12, double Azero)
{
   TYd = TYd_Azero_coeff * Azero + TYd_m12_coeff * m12;
   ThE = ThE_Azero_coeff * Azero + ThE_m12_coeff * m12;
   TYe = TYe_Azero_coeff * Azero + TYe_m12_coeff * m12;
   TSigmaL = TSigmaL_Azero_coeff * Azero + TSigmaL_m12_coeff * m12;
   TKappaPr = TKappaPr_Azero_coeff * Azero + TKappaPr_m12_coeff * m12;
   TSigmax = TSigmax_Azero_coeff * Azero + TSigmax_m12_coeff * m12;
   TgD = TgD_Azero_coeff * Azero + TgD_m12_coeff * m12;
   TKappa = TKappa_Azero_coeff * Azero + TKappa_m12_coeff * m12;
   TLambda12 = TLambda12_Azero_coeff * Azero + TLambda12_m12_coeff * m12;
   TLambdax = TLambdax_Azero_coeff * Azero + TLambdax_m12_coeff * m12;
   Tfu = Tfu_Azero_coeff * Azero + Tfu_m12_coeff * m12;
   Tfd = Tfd_Azero_coeff * Azero + Tfd_m12_coeff * m12;
   mq2 = mq2_m02_coeff * Sqr(m0) + mq2_m122_coeff * Sqr(m12)
      + mq2_Azerom12_coeff * Azero * m12 + mq2_Azero2_coeff * Sqr(Azero);
   ml2 = ml2_m02_coeff * Sqr(m0) + ml2_m122_coeff * Sqr(m12)
      + ml2_Azerom12_coeff * Azero * m12 + ml2_Azero2_coeff * Sqr(Azero);
   mHd2 = mHd2_m02_coeff * Sqr(m0) + mHd2_m122_coeff * Sqr(m12)
      + mHd2_Azerom12_coeff * Azero * m12 + mHd2_Azero2_coeff * Sqr(Azero);
   mHu2 = mHu2_m02_coeff * Sqr(m0) + mHu2_m122_coeff * Sqr(m12)
      + mHu2_Azerom12_coeff * Azero * m12 + mHu2_Azero2_coeff * Sqr(Azero);
   md2 = md2_m02_coeff * Sqr(m0) + md2_m122_coeff * Sqr(m12)
      + md2_Azerom12_coeff * Azero * m12 + md2_Azero2_coeff * Sqr(Azero);
   mu2 = mu2_m02_coeff * Sqr(m0) + mu2_m122_coeff * Sqr(m12)
      + mu2_Azerom12_coeff * Azero * m12 + mu2_Azero2_coeff * Sqr(Azero);
   me2 = me2_m02_coeff * Sqr(m0) + me2_m122_coeff * Sqr(m12)
      + me2_Azerom12_coeff * Azero * m12 + me2_Azero2_coeff * Sqr(Azero);
   ms2 = ms2_m02_coeff * Sqr(m0) + ms2_m122_coeff * Sqr(m12)
      + ms2_Azerom12_coeff * Azero * m12 + ms2_Azero2_coeff * Sqr(Azero);
   msbar2 = msbar2_m02_coeff * Sqr(m0) + msbar2_m122_coeff * Sqr(m12)
      + msbar2_Azerom12_coeff * Azero * m12 + msbar2_Azero2_coeff * Sqr(Azero);
   mH1I2 = mH1I2_m02_coeff * Sqr(m0) + mH1I2_m122_coeff * Sqr(m12)
      + mH1I2_Azerom12_coeff * Azero * m12 + mH1I2_Azero2_coeff * Sqr(Azero);
   mH2I2 = mH2I2_m02_coeff * Sqr(m0) + mH2I2_m122_coeff * Sqr(m12)
      + mH2I2_Azerom12_coeff * Azero * m12 + mH2I2_Azero2_coeff * Sqr(Azero);;
   mSI2 = mSI2_m02_coeff * Sqr(m0) + mSI2_m122_coeff * Sqr(m12)
      + mSI2_Azerom12_coeff * Azero * m12 + mSI2_Azero2_coeff * Sqr(Azero);
   mDx2 = mDx2_m02_coeff * Sqr(m0) + mDx2_m122_coeff * Sqr(m12)
      + mDx2_Azerom12_coeff * Azero * m12 + mDx2_Azero2_coeff * Sqr(Azero);
   mDxbar2 = mDxbar2_m02_coeff * Sqr(m0) + mDxbar2_m122_coeff * Sqr(m12)
      + mDxbar2_Azerom12_coeff * Azero * m12 + mDxbar2_Azero2_coeff * Sqr(Azero);
   mHp2 = mHp2_m02_coeff * Sqr(m0) + mHp2_m122_coeff * Sqr(m12)
      + mHp2_Azerom12_coeff * Azero * m12 + mHp2_Azero2_coeff * Sqr(Azero);
   mHpbar2 = mHpbar2_m02_coeff * Sqr(m0) + mHpbar2_m122_coeff * Sqr(m12)
      + mHpbar2_Azerom12_coeff * Azero * m12 + mHpbar2_Azero2_coeff * Sqr(Azero);
   mphi2 = mphi2_m02_coeff * Sqr(m0) + mphi2_m122_coeff * Sqr(m12)
      + mphi2_Azerom12_coeff * Azero * m12 + mphi2_Azero2_coeff * Sqr(Azero);
   MassB = MassB_Azero_coeff * Azero + MassB_m12_coeff * m12;
   MassWB = MassWB_Azero_coeff * Azero + MassWB_m12_coeff * m12;
   MassG = MassG_Azero_coeff * Azero + MassG_m12_coeff * m12;
   MassBp = MassBp_Azero_coeff * Azero + MassBp_m12_coeff * m12;
}

std::ostream& operator<<(std::ostream& ostr, const CNE6SSM<Semianalytic>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
