// ====================================================================
// Test implementation of low scale constraint for semianalytic
// algorithm
// ====================================================================

#ifndef CNE6SSM_SEMIANALYTIC_LOW_SCALE_CONSTRAINT_H
#define CNE6SSM_SEMIANALYTIC_LOW_SCALE_CONSTRAINT_H

#include "CNE6SSM_low_scale_constraint.hpp"
#include "semianalytic_constraint.hpp"
#include "lowe.h"
#include <Eigen/Core>

namespace flexiblesusy {

template <class T>
class CNE6SSM;

class Semianalytic;

template<>
class CNE6SSM_low_scale_constraint<Semianalytic> : public Constraint<Semianalytic> {
public:
   CNE6SSM_low_scale_constraint();
   CNE6SSM_low_scale_constraint(CNE6SSM<Semianalytic>*, const QedQcd&);
   virtual ~CNE6SSM_low_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Semianalytic_model*);

   void clear();
   const Eigen::Matrix<std::complex<double>,3,3>& get_ckm() const;
   const Eigen::Matrix<std::complex<double>,3,3>& get_pmns() const;
   double get_initial_scale_guess() const;
   void initialize();
   const QedQcd& get_sm_parameters() const;
   void set_sm_parameters(const QedQcd&);
   void set_threshold_corrections_loop_order(unsigned); ///< threshold corrections loop order
   void apply_susy_constraint_only(bool);

private:
   double scale;
   double initial_scale_guess;
   CNE6SSM<Semianalytic>* model;
   QedQcd oneset;
   Eigen::Matrix<std::complex<double>,3,3> ckm;
   Eigen::Matrix<std::complex<double>,3,3> pmns;
   Eigen::Matrix<double,3,3> neutrinoDRbar;
   double MWDRbar;
   double MZDRbar;
   double AlphaS;
   double EDRbar;
   double ThetaWDRbar;
   double new_g1, new_g2, new_g3;
   unsigned threshold_corrections_loop_order; ///< threshold corrections loop order
   bool only_apply_susy_constraint; ///< only apply constraint for SUSY parameters

   void calculate_threshold_corrections();
   void calculate_DRbar_gauge_couplings();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
   void calculate_MNeutrino_DRbar();
   double calculate_delta_alpha_em(double) const;
   double calculate_delta_alpha_s(double) const;
   void update_scale();
};

} // namespace flexiblesusy

#endif
