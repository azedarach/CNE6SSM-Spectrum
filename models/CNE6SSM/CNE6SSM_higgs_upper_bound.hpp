// ====================================================================
// Computes an approximate 1-loop upper bound on the lightest CP-even
// Higgs
// ====================================================================

#ifndef CNE6SSM_HIGGS_UPPER_BOUND_H
#define CNE6SSM_HIGGS_UPPER_BOUND_H

#include "CNE6SSM_soft_parameters.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

class CNE6SSM_higgs_upper_bound {
public:
   CNE6SSM_higgs_upper_bound(const CNE6SSM_soft_parameters&);
   ~CNE6SSM_higgs_upper_bound();

   // if false (default), include only the third
   // generation contributions
   void set_include_all_SM_generations(bool flag) { include_all_gens = flag; }
   void set_include_up_tadpoles(bool flag) { include_ups = flag; }
   void set_include_exotic_tadpoles(bool flag) { include_exotics = flag; }
   void set_include_inert_singlet_tadpoles(bool flag) { include_inert_singlets = flag; }
   void set_include_inert_neutral_higgs_tadpoles(bool flag) {
      include_inert_neutral_higgs = flag;
   }
   void set_include_inert_charged_higgs_tadpoles(bool flag) {
      include_inert_charged_higgs = flag;
   }

   double get_upper_bound() const { return upper_bound; }

   double calculate_upper_bound();

   double get_tadpole_vd() const;
   double get_tadpole_vu() const;
   double get_up_contribution(unsigned) const;
   double get_exotic_contribution(unsigned) const;
   double get_inert_singlet_contribution(unsigned) const;
   double get_inert_neutral_higgs_contribution(unsigned) const;
   double get_inert_charged_higgs_contribution(unsigned) const;

   double get_unrotated_up_contribution(unsigned,unsigned,unsigned) const;
   double get_unrotated_exotic_contribution(unsigned,unsigned,unsigned) const;
   double get_unrotated_inert_singlet_contribution(unsigned,unsigned,unsigned) const;
   double get_unrotated_inert_neutral_higgs_contribution(unsigned,unsigned,unsigned) const;
   double get_unrotated_inert_charged_higgs_contribution(unsigned,unsigned,unsigned) const;

private:
   CNE6SSM_soft_parameters model;
   double upper_bound;
   bool include_all_gens;
   bool include_ups;
   bool include_exotics;
   bool include_inert_singlets;
   bool include_inert_neutral_higgs;
   bool include_inert_charged_higgs;

   Eigen::Matrix<double,2,2> get_mass_matrix_Su(unsigned) const;
   Eigen::Array<double,2,1> calculate_MSu2(unsigned) const;
   double calculate_Sin2ThetaSu(unsigned) const;
   double calculate_Cos2ThetaSu(unsigned) const;
   double calculate_MFu2(unsigned) const;

   Eigen::Matrix<double,2,2> get_mass_matrix_SDX(unsigned) const;
   Eigen::Array<double,2,1> calculate_MSDX2(unsigned) const;
   double calculate_Sin2ThetaSDX(unsigned) const;
   double calculate_Cos2ThetaSDX(unsigned) const;
   double calculate_MFDX2(unsigned) const;

   Eigen::Matrix<double,2,2> get_mass_matrix_HI0(unsigned) const;
   Eigen::Array<double,2,1> calculate_MHI02(unsigned) const;
   double calculate_Sin2ThetaHI0(unsigned) const;
   double calculate_Cos2ThetaHI0(unsigned) const;
   double calculate_MFHI02(unsigned) const;

   double calculate_MSI02(unsigned gen) const;

   Eigen::Matrix<double,2,2> get_mass_matrix_HIPM(unsigned) const;
   Eigen::Array<double,2,1> calculate_MHIPM2(unsigned) const;
   double calculate_Sin2ThetaHIPM(unsigned) const;
   double calculate_Cos2ThetaHIPM(unsigned) const;
   double calculate_MFHIPM2(unsigned) const;

   // convenient helper functions for derivatives
   Eigen::Matrix<double,2,2> get_dmass_matrix_Su_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_dmass_matrix_Su_dvu(unsigned) const;

   Eigen::Matrix<double,2,2> get_dmass_matrix_SDX_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_dmass_matrix_SDX_dvu(unsigned) const;

   Eigen::Matrix<double,2,2> get_dmass_matrix_HI0_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_dmass_matrix_HI0_dvu(unsigned) const;

   Eigen::Matrix<double,2,2> get_dmass_matrix_HIPM_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_dmass_matrix_HIPM_dvu(unsigned) const;

   double get_dV1lp_up_dvd(unsigned) const;
   double get_dV1lp_up_dvu(unsigned) const;

   double get_dV1lp_exotic_dvd(unsigned) const;
   double get_dV1lp_exotic_dvu(unsigned) const;

   double get_dV1lp_inert_singlet_dvd(unsigned) const;
   double get_dV1lp_inert_singlet_dvu(unsigned) const;

   double get_dV1lp_inert_neutral_higgs_dvd(unsigned) const;
   double get_dV1lp_inert_neutral_higgs_dvu(unsigned) const;

   double get_dV1lp_inert_charged_higgs_dvd(unsigned) const;
   double get_dV1lp_inert_charged_higgs_dvu(unsigned) const;
};

} // namespace flexiblesusy

#endif
