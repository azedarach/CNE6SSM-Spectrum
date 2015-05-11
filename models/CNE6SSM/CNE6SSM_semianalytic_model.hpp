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

   const Eigen::Matrix<double,3,3>& get_Azero_coeff_TYd() const { return TYd_Azero_coeff; }
   double get_Azero_coeff_TYd(int i, int k) const { return TYd_Azero_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m12_coeff_TYd() const { return TYd_m12_coeff; }
   double get_m12_coeff_TYd(int i, int k) const { return TYd_m12_coeff(i,k); }
   const Eigen::Matrix<double,3,2>& get_Azero_coeff_ThE() const { return ThE_Azero_coeff; }
   double get_Azero_coeff_ThE(int i, int k) const { return ThE_Azero_coeff(i,k); }
   const Eigen::Matrix<double,3,2>& get_m12_coeff_ThE() const { return ThE_m12_coeff; }
   double get_m12_coeff_ThE(int i, int k) const { return ThE_m12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero_coeff_TYe() const { return TYe_Azero_coeff; }
   double get_Azero_coeff_TYe(int i, int k) const { return TYe_Azero_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m12_coeff_TYe() const { return TYe_m12_coeff; }
   double get_m12_coeff_TYe(int i, int k) const { return TYe_m12_coeff(i,k); }
   double get_Azero_coeff_TSigmaL() const { return TSigmaL_Azero_coeff; }
   double get_m12_coeff_TSigmaL() const { return TSigmaL_m12_coeff; }
   double get_Azero_coeff_TKappaPr() const { return TKappaPr_Azero_coeff; }
   double get_m12_coeff_TKappaPr() const { return TKappaPr_m12_coeff; }
   double get_Azero_coeff_TSigmax() const { return TSigmax_Azero_coeff; }
   double get_m12_coeff_TSigmax() const { return TSigmax_m12_coeff; }
   const Eigen::Matrix<double,3,3>& get_Azero_coeff_TgD() const { return TgD_Azero_coeff; }
   double get_Azero_coeff_TgD(int i, int k) const { return TgD_Azero_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m12_coeff_TgD() const { return TgD_m12_coeff; }
   double get_m12_coeff_TgD(int i, int k) const { return TgD_m12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero_coeff_TKappa() const { return TKappa_Azero_coeff; }
   double get_Azero_coeff_TKappa(int i, int k) const { return TKappa_Azero_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m12_coeff_TKappa() const { return TKappa_m12_coeff; }
   double get_m12_coeff_TKappa(int i, int k) const { return TKappa_m12_coeff(i,k); }
   const Eigen::Matrix<double,2,2>& get_Azero_coeff_TLambda12() const { return TLambda12_Azero_coeff; }
   double get_Azero_coeff_TLambda12(int i, int k) const { return TLambda12_Azero_coeff(i,k); }
   const Eigen::Matrix<double,2,2>& get_m12_coeff_TLambda12() const { return TLambda12_m12_coeff; }
   double get_m12_coeff_TLambda12(int i, int k) const { return TLambda12_m12_coeff(i,k); }
   double get_Azero_coeff_TLambdax() const { return TLambdax_Azero_coeff; }
   double get_m12_coeff_TLambdax() const { return TLambdax_m12_coeff; }
   const Eigen::Matrix<double,3,2>& get_Azero_coeff_Tfu() const { return Tfu_Azero_coeff; }
   double get_Azero_coeff_Tfu(int i, int k) const { return Tfu_Azero_coeff(i,k); }
   const Eigen::Matrix<double,3,2>& get_m12_coeff_Tfu() const { return Tfu_m12_coeff; }
   double get_m12_coeff_Tfu(int i, int k) const { return Tfu_m12_coeff(i,k); }
   const Eigen::Matrix<double,3,2>& get_Azero_coeff_Tfd() const { return Tfd_Azero_coeff; }
   double get_Azero_coeff_Tfd(int i, int k) const { return Tfd_Azero_coeff(i,k); }
   const Eigen::Matrix<double,3,2>& get_m12_coeff_Tfd() const { return Tfd_m12_coeff; }
   double get_m12_coeff_Tfd(int i, int k) const { return Tfd_m12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero_coeff_TYu() const { return TYu_Azero_coeff; }
   double get_Azero_coeff_TYu(int i, int k) const { return TYu_Azero_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m12_coeff_TYu() const { return TYu_m12_coeff; }
   double get_m12_coeff_TYu(int i, int k) const { return TYu_m12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_mq2() const { return mq2_m02_coeff; }
   double get_m02_coeff_mq2(int i, int k) const { return mq2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_mq2() const { return mq2_m122_coeff; }
   double get_m122_coeff_mq2(int i, int k) const { return mq2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_mq2() const { return mq2_Azerom12_coeff; }
   double get_Azerom12_coeff_mq2(int i, int k) const { return mq2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_mq2() const { return mq2_Azero2_coeff; }
   double get_Azero2_coeff_mq2(int i, int k) const { return mq2_Azero2_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_ml2() const { return ml2_m02_coeff; }
   double get_m02_coeff_ml2(int i, int k) const { return ml2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_ml2() const { return ml2_m122_coeff; }
   double get_m122_coeff_ml2(int i, int k) const { return ml2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_ml2() const { return ml2_Azerom12_coeff; }
   double get_Azerom12_coeff_ml2(int i, int k) const { return ml2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_ml2() const { return ml2_Azero2_coeff; }
   double get_Azero2_coeff_ml2(int i, int k) const { return ml2_Azero2_coeff(i,k); }
   double get_m02_coeff_mHd2() const { return mHd2_m02_coeff; }
   double get_m122_coeff_mHd2() const { return mHd2_m122_coeff; }
   double get_Azerom12_coeff_mHd2() const { return mHd2_Azerom12_coeff; }
   double get_Azero2_coeff_mHd2() const { return mHd2_Azero2_coeff; }
   double get_m02_coeff_mHu2() const { return mHu2_m02_coeff; }
   double get_m122_coeff_mHu2() const { return mHu2_m122_coeff; }
   double get_Azerom12_coeff_mHu2() const { return mHu2_Azerom12_coeff; }
   double get_Azero2_coeff_mHu2() const { return mHu2_Azero2_coeff; }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_md2() const { return md2_m02_coeff; }
   double get_m02_coeff_md2(int i, int k) const { return md2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_md2() const { return md2_m122_coeff; }
   double get_m122_coeff_md2(int i, int k) const { return md2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_md2() const { return md2_Azerom12_coeff; }
   double get_Azerom12_coeff_md2(int i, int k) const { return md2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_md2() const { return md2_Azero2_coeff; }
   double get_Azero2_coeff_md2(int i, int k) const { return md2_Azero2_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_mu2() const { return mu2_m02_coeff; }
   double get_m02_coeff_mu2(int i, int k) const { return mu2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_mu2() const { return mu2_m122_coeff; }
   double get_m122_coeff_mu2(int i, int k) const { return mu2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_mu2() const { return mu2_Azerom12_coeff; }
   double get_Azerom12_coeff_mu2(int i, int k) const { return mu2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_mu2() const { return mu2_Azero2_coeff; }
   double get_Azero2_coeff_mu2(int i, int k) const { return mu2_Azero2_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_me2() const { return me2_m02_coeff; }
   double get_m02_coeff_me2(int i, int k) const { return me2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_me2() const { return me2_m122_coeff; }
   double get_m122_coeff_me2(int i, int k) const { return me2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_me2() const { return me2_Azerom12_coeff; }
   double get_Azerom12_coeff_me2(int i, int k) const { return me2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_me2() const { return me2_Azero2_coeff; }
   double get_Azero2_coeff_me2(int i, int k) const { return me2_Azero2_coeff(i,k); }
   double get_m02_coeff_ms2() const { return ms2_m02_coeff; }
   double get_m122_coeff_ms2() const { return ms2_m122_coeff; }
   double get_Azerom12_coeff_ms2() const { return ms2_Azerom12_coeff; }
   double get_Azero2_coeff_ms2() const { return ms2_Azero2_coeff; }
   double get_m02_coeff_msbar2() const { return msbar2_m02_coeff; }
   double get_m122_coeff_msbar2() const { return msbar2_m122_coeff; }
   double get_Azerom12_coeff_msbar2() const { return msbar2_Azerom12_coeff; }
   double get_Azero2_coeff_msbar2() const { return msbar2_Azero2_coeff; }
   const Eigen::Matrix<double,2,2>& get_m02_coeff_mH1I2() const { return mH1I2_m02_coeff; }
   double get_m02_coeff_mH1I2(int i, int k) const { return mH1I2_m02_coeff(i,k); }
   const Eigen::Matrix<double,2,2>& get_m122_coeff_mH1I2() const { return mH1I2_m122_coeff; }
   double get_m122_coeff_mH1I2(int i, int k) const { return mH1I2_m122_coeff(i,k); }
   const Eigen::Matrix<double,2,2>& get_Azerom12_coeff_mH1I2() const { return mH1I2_Azerom12_coeff; }
   double get_Azerom12_coeff_mH1I2(int i, int k) const { return mH1I2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,2,2>& get_Azero2_coeff_mH1I2() const { return mH1I2_Azero2_coeff; }
   double get_Azero2_coeff_mH1I2(int i, int k) const { return mH1I2_Azero2_coeff(i,k); }
   const Eigen::Matrix<double,2,2>& get_m02_coeff_mH2I2() const { return mH2I2_m02_coeff; }
   double get_m02_coeff_mH2I2(int i, int k) const { return mH2I2_m02_coeff(i,k); }
   const Eigen::Matrix<double,2,2>& get_m122_coeff_mH2I2() const { return mH2I2_m122_coeff; }
   double get_m122_coeff_mH2I2(int i, int k) const { return mH2I2_m122_coeff(i,k); }
   const Eigen::Matrix<double,2,2>& get_Azerom12_coeff_mH2I2() const { return mH2I2_Azerom12_coeff; }
   double get_Azerom12_coeff_mH2I2(int i, int k) const { return mH2I2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,2,2>& get_Azero2_coeff_mH2I2() const { return mH2I2_Azero2_coeff; }
   double get_Azero2_coeff_mH2I2(int i, int k) const { return mH2I2_Azero2_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_mSI2() const { return mSI2_m02_coeff; }
   double get_m02_coeff_mSI2(int i, int k) const { return mSI2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_mSI2() const { return mSI2_m122_coeff; }
   double get_m122_coeff_mSI2(int i, int k) const { return mSI2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_mSI2() const { return mSI2_Azerom12_coeff; }
   double get_Azerom12_coeff_mSI2(int i, int k) const { return mSI2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_mSI2() const { return mSI2_Azero2_coeff; }
   double get_Azero2_coeff_mSI2(int i, int k) const { return mSI2_Azero2_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_mDx2() const { return mDx2_m02_coeff; }
   double get_m02_coeff_mDx2(int i, int k) const { return mDx2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_mDx2() const { return mDx2_m122_coeff; }
   double get_m122_coeff_mDx2(int i, int k) const { return mDx2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_mDx2() const { return mDx2_Azerom12_coeff; }
   double get_Azerom12_coeff_mDx2(int i, int k) const { return mDx2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_mDx2() const { return mDx2_Azero2_coeff; }
   double get_Azero2_coeff_mDx2(int i, int k) const { return mDx2_Azero2_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_mDxbar2() const { return mDxbar2_m02_coeff; }
   double get_m02_coeff_mDxbar2(int i, int k) const { return mDxbar2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_mDxbar2() const { return mDxbar2_m122_coeff; }
   double get_m122_coeff_mDxbar2(int i, int k) const { return mDxbar2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_mDxbar2() const { return mDxbar2_Azerom12_coeff; }
   double get_Azerom12_coeff_mDxbar2(int i, int k) const { return mDxbar2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_mDxbar2() const { return mDxbar2_Azero2_coeff; }
   double get_Azero2_coeff_mDxbar2(int i, int k) const { return mDxbar2_Azero2_coeff(i,k); }
   double get_m02_coeff_mHp2() const { return mHp2_m02_coeff; }
   double get_m122_coeff_mHp2() const { return mHp2_m122_coeff; }
   double get_Azerom12_coeff_mHp2() const { return mHp2_Azerom12_coeff; }
   double get_Azero2_coeff_mHp2() const { return mHp2_Azero2_coeff; }
   double get_m02_coeff_mHpbar2() const { return mHpbar2_m02_coeff; }
   double get_m122_coeff_mHpbar2() const { return mHpbar2_m122_coeff; }
   double get_Azerom12_coeff_mHpbar2() const { return mHpbar2_Azerom12_coeff; }
   double get_Azero2_coeff_mHpbar2() const { return mHpbar2_Azero2_coeff; }
   double get_m02_coeff_mphi2() const { return mphi2_m02_coeff; }
   double get_m122_coeff_mphi2() const { return mphi2_m122_coeff; }
   double get_Azerom12_coeff_mphi2() const { return mphi2_Azerom12_coeff; }
   double get_Azero2_coeff_mphi2() const { return mphi2_Azero2_coeff; }
   double get_Azero_coeff_MassB() const { return MassB_Azero_coeff; }
   double get_m12_coeff_MassB() const { return MassB_m12_coeff; }
   double get_Azero_coeff_MassWB() const { return MassWB_Azero_coeff; }
   double get_m12_coeff_MassWB() const { return MassWB_m12_coeff; }
   double get_Azero_coeff_MassG() const { return MassG_Azero_coeff; }
   double get_m12_coeff_MassG() const { return MassG_m12_coeff; }
   double get_Azero_coeff_MassBp() const { return MassBp_Azero_coeff; }
   double get_m12_coeff_MassBp() const { return MassBp_m12_coeff; }

   void calculate_coefficients(double);

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
   static const std::size_t number_of_fit_points = 4;

   Eigen::Array<double,number_of_ewsb_equations,1> ewsb_solution;

   int solve_ewsb_iteratively();
   int solve_ewsb_iteratively(unsigned);
   int solve_ewsb_iteratively_with(EWSB_solver*, const double[number_of_ewsb_equations]);
   int check_ewsb_solution(double);
   void ewsb_initial_guess(double[number_of_ewsb_equations]);
   int ewsb_step(double[number_of_ewsb_equations]);
   static int ewsb_step(const gsl_vector*, void*, gsl_vector*);
   void tadpole_equations(double[number_of_ewsb_equations]) const;
   static int tadpole_equations(const gsl_vector*, void*, gsl_vector*);

   void set_soft_parameters_at_input_scale(double,double,double);
   void set_soft_parameters_at_current_scale(double,double,double);

   Eigen::Matrix<double,3,3> TYd_Azero_coeff;
   Eigen::Matrix<double,3,3> TYd_m12_coeff;
   Eigen::Matrix<double,3,2> ThE_Azero_coeff;
   Eigen::Matrix<double,3,2> ThE_m12_coeff;
   Eigen::Matrix<double,3,3> TYe_Azero_coeff;
   Eigen::Matrix<double,3,3> TYe_m12_coeff;
   double TSigmaL_Azero_coeff;
   double TSigmaL_m12_coeff;
   double TKappaPr_Azero_coeff;
   double TKappaPr_m12_coeff;
   double TSigmax_Azero_coeff;
   double TSigmax_m12_coeff;
   Eigen::Matrix<double,3,3> TgD_Azero_coeff;
   Eigen::Matrix<double,3,3> TgD_m12_coeff;
   Eigen::Matrix<double,3,3> TKappa_Azero_coeff;
   Eigen::Matrix<double,3,3> TKappa_m12_coeff;
   Eigen::Matrix<double,2,2> TLambda12_Azero_coeff;
   Eigen::Matrix<double,2,2> TLambda12_m12_coeff;
   double TLambdax_Azero_coeff;
   double TLambdax_m12_coeff;
   Eigen::Matrix<double,3,2> Tfu_Azero_coeff;
   Eigen::Matrix<double,3,2> Tfu_m12_coeff;
   Eigen::Matrix<double,3,2> Tfd_Azero_coeff;
   Eigen::Matrix<double,3,2> Tfd_m12_coeff;
   Eigen::Matrix<double,3,3> TYu_Azero_coeff;
   Eigen::Matrix<double,3,3> TYu_m12_coeff;
   Eigen::Matrix<double,3,3> mq2_m02_coeff;
   Eigen::Matrix<double,3,3> mq2_m122_coeff;
   Eigen::Matrix<double,3,3> mq2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> mq2_Azero2_coeff;
   Eigen::Matrix<double,3,3> ml2_m02_coeff;
   Eigen::Matrix<double,3,3> ml2_m122_coeff;
   Eigen::Matrix<double,3,3> ml2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> ml2_Azero2_coeff;
   double mHd2_m02_coeff;
   double mHd2_m122_coeff;
   double mHd2_Azerom12_coeff;
   double mHd2_Azero2_coeff;
   double mHu2_m02_coeff;
   double mHu2_m122_coeff;
   double mHu2_Azerom12_coeff;
   double mHu2_Azero2_coeff;
   Eigen::Matrix<double,3,3> md2_m02_coeff;
   Eigen::Matrix<double,3,3> md2_m122_coeff;
   Eigen::Matrix<double,3,3> md2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> md2_Azero2_coeff;
   Eigen::Matrix<double,3,3> mu2_m02_coeff;
   Eigen::Matrix<double,3,3> mu2_m122_coeff;
   Eigen::Matrix<double,3,3> mu2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> mu2_Azero2_coeff;
   Eigen::Matrix<double,3,3> me2_m02_coeff;
   Eigen::Matrix<double,3,3> me2_m122_coeff;
   Eigen::Matrix<double,3,3> me2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> me2_Azero2_coeff;
   double ms2_m02_coeff;
   double ms2_m122_coeff;
   double ms2_Azerom12_coeff;
   double ms2_Azero2_coeff;
   double msbar2_m02_coeff;
   double msbar2_m122_coeff;
   double msbar2_Azerom12_coeff;
   double msbar2_Azero2_coeff;
   Eigen::Matrix<double,2,2> mH1I2_m02_coeff;
   Eigen::Matrix<double,2,2> mH1I2_m122_coeff;
   Eigen::Matrix<double,2,2> mH1I2_Azerom12_coeff;
   Eigen::Matrix<double,2,2> mH1I2_Azero2_coeff;
   Eigen::Matrix<double,2,2> mH2I2_m02_coeff;
   Eigen::Matrix<double,2,2> mH2I2_m122_coeff;
   Eigen::Matrix<double,2,2> mH2I2_Azerom12_coeff;
   Eigen::Matrix<double,2,2> mH2I2_Azero2_coeff;
   Eigen::Matrix<double,3,3> mSI2_m02_coeff;
   Eigen::Matrix<double,3,3> mSI2_m122_coeff;
   Eigen::Matrix<double,3,3> mSI2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> mSI2_Azero2_coeff;
   Eigen::Matrix<double,3,3> mDx2_m02_coeff;
   Eigen::Matrix<double,3,3> mDx2_m122_coeff;
   Eigen::Matrix<double,3,3> mDx2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> mDx2_Azero2_coeff;
   Eigen::Matrix<double,3,3> mDxbar2_m02_coeff;
   Eigen::Matrix<double,3,3> mDxbar2_m122_coeff;
   Eigen::Matrix<double,3,3> mDxbar2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> mDxbar2_Azero2_coeff;
   double mHp2_m02_coeff;
   double mHp2_m122_coeff;
   double mHp2_Azerom12_coeff;
   double mHp2_Azero2_coeff;
   double mHpbar2_m02_coeff;
   double mHpbar2_m122_coeff;
   double mHpbar2_Azerom12_coeff;
   double mHpbar2_Azero2_coeff;
   double mphi2_m02_coeff;
   double mphi2_m122_coeff;
   double mphi2_Azerom12_coeff;
   double mphi2_Azero2_coeff;
   double MassB_Azero_coeff;
   double MassB_m12_coeff;
   double MassWB_Azero_coeff;
   double MassWB_m12_coeff;
   double MassG_Azero_coeff;
   double MassG_m12_coeff;
   double MassBp_Azero_coeff;
   double MassBp_m12_coeff;
};

   std::ostream& operator<<(std::ostream&, const CNE6SSM<Semianalytic>&);

} // namespace flexiblesusy

#endif
