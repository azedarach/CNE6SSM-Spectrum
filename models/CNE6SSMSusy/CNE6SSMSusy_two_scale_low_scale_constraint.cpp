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

// File generated at Sun 19 Apr 2015 20:31:40

#include "CNE6SSMSusy_two_scale_low_scale_constraint.hpp"
#include "CNE6SSMSusy_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "weinberg_angle.hpp"

#include <cassert>
#include <cmath>
#include <limits>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME CNE6SSMSusy<Two_scale>
#define CKM ckm
#define PMNS pmns
#define THETAW theta_w
#define ALPHA_EM_DRBAR alpha_em_drbar
#define CALCULATE_DRBAR_MASSES() model->calculate_DRbar_masses()

CNE6SSMSusy_low_scale_constraint<Two_scale>::CNE6SSMSusy_low_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , oneset()
   , ckm()
   , pmns()
   , MWDRbar(0.)
   , MZDRbar(0.)
   , AlphaS(0.)
   , EDRbar(0.)
   , ThetaWDRbar(0.)
   , new_g1(0.)
   , new_g2(0.)
   , new_g3(0.)
   , self_energy_w_at_mw(0.)
   , threshold_corrections_loop_order(1)
{
   ckm << 1., 0., 0.,
          0., 1., 0.,
          0., 0., 1.;

   pmns << 1., 0., 0.,
           0., 1., 0.,
           0., 0., 1.;
}

CNE6SSMSusy_low_scale_constraint<Two_scale>::CNE6SSMSusy_low_scale_constraint(
   CNE6SSMSusy<Two_scale>* model_, const QedQcd& oneset_)
   : Constraint<Two_scale>()
   , model(model_)
   , oneset(oneset_)
   , new_g1(0.)
   , new_g2(0.)
   , new_g3(0.)
   , self_energy_w_at_mw(0.)
{
   initialize();
}

CNE6SSMSusy_low_scale_constraint<Two_scale>::~CNE6SSMSusy_low_scale_constraint()
{
}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: CNE6SSMSusy_low_scale_constraint::apply():"
          " model pointer must not be zero");

   update_scale();
   calculate_DRbar_gauge_couplings();

   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto QSInput = INPUTPARAMETER(QSInput);
   const auto XiFInput = INPUTPARAMETER(XiFInput);

   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
   MODEL->set_vd((2*MZDRbar)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)
      )));
   MODEL->set_vu((2*MZDRbar*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(
      TanBeta))));
   MODEL->set_QS(QSInput);
   MODEL->set_XiF(XiFInput);

   model->set_g1(new_g1);
   model->set_g2(new_g2);
   model->set_g3(new_g3);

   recalculate_mw_pole();
}

const Eigen::Matrix<std::complex<double>,3,3>& CNE6SSMSusy_low_scale_constraint<Two_scale>::get_ckm()
{
   return ckm;
}

const Eigen::Matrix<std::complex<double>,3,3>& CNE6SSMSusy_low_scale_constraint<Two_scale>::get_pmns()
{
   return pmns;
}

double CNE6SSMSusy_low_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double CNE6SSMSusy_low_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<CNE6SSMSusy<Two_scale>*>(model_);
}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::set_sm_parameters(const QedQcd& oneset_)
{
   oneset = oneset_;
}

const QedQcd& CNE6SSMSusy_low_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return oneset;
}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
   oneset = QedQcd();
   MWDRbar = 0.;
   MZDRbar = 0.;
   AlphaS = 0.;
   EDRbar = 0.;
   ThetaWDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
   self_energy_w_at_mw = 0.;
}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::initialize()
{
   assert(model && "CNE6SSMSusy_low_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   initial_scale_guess = LowEnergyConstant(MZ);

   scale = initial_scale_guess;

   MWDRbar = 0.;
   MZDRbar = 0.;
   AlphaS = 0.;
   EDRbar = 0.;
   ThetaWDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
   ckm = oneset.get_complex_ckm();
   pmns = oneset.get_complex_pmns();
   self_energy_w_at_mw = 0.;
}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "CNE6SSMSusy_low_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   scale = LowEnergyConstant(MZ);


}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::calculate_threshold_corrections()
{
   assert(oneset.displayMu() == get_scale() && "Error: low-energy data"
          " set is not defined at the same scale as the low-energy"
          " constraint.  You need to run the low-energy data set to this"
          " scale!");
   assert(model && "CNE6SSMSusy_low_scale_constraint<Two_scale>::"
          "calculate_threshold_corrections(): model pointer is zero");

   const double alpha_em = oneset.displayAlpha(ALPHA);
   const double alpha_s  = oneset.displayAlpha(ALPHAS);
   const double mw_pole  = oneset.displayPoleMW();
   const double mz_pole  = oneset.displayPoleMZ();

   double delta_alpha_em = 0.;
   double delta_alpha_s  = 0.;

   if (model->get_thresholds()) {
      // no thresholds
   }

   const double alpha_em_drbar = alpha_em / (1.0 - delta_alpha_em);
   const double alpha_s_drbar  = alpha_s  / (1.0 - delta_alpha_s);
   const double e_drbar        = Sqrt(4.0 * Pi * alpha_em_drbar);

   // interface variables
   MZDRbar = mz_pole;
   MWDRbar = mw_pole;

   // no thresholds

   AlphaS = alpha_s_drbar;
   EDRbar = e_drbar;
   ThetaWDRbar = calculate_theta_w(alpha_em_drbar);
}

double CNE6SSMSusy_low_scale_constraint<Two_scale>::calculate_theta_w(double alpha_em_drbar)
{
   assert(model && "CNE6SSMSusy_low_scale_constraint<Two_scale>::"
          "calculate_theta_w(): model pointer is zero");

   double theta_w = 0.;

   THETAW = ArcSin(Sqrt(1 - Sqr(MWDRbar)/Sqr(MZDRbar)));


   return theta_w;
}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::calculate_DRbar_gauge_couplings()
{
   assert(model && "CNE6SSMSusy_low_scale_constraint<Two_scale>::"
          "calculate_DRbar_gauge_couplings(): model pointer is zero");

   calculate_threshold_corrections();

   new_g1 = 1.2909944487358056*EDRbar*Sec(ThetaWDRbar);
   new_g2 = EDRbar*Csc(ThetaWDRbar);
   new_g3 = 3.5449077018110318*Sqrt(AlphaS);

}

double CNE6SSMSusy_low_scale_constraint<Two_scale>::calculate_delta_alpha_em(double alphaEm) const
{
   assert(model && "CNE6SSMSusy_low_scale_constraint<Two_scale>::"
          "calculate_delta_alpha_em(): model pointer is zero");

   // no thresholds

   return 0.;

}

double CNE6SSMSusy_low_scale_constraint<Two_scale>::calculate_delta_alpha_s(double alphaS) const
{
   assert(model && "CNE6SSMSusy_low_scale_constraint<Two_scale>::"
          "calculate_delta_alpha_s(): model pointer is zero");

   // no thresholds

   return 0.;

}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::calculate_DRbar_yukawa_couplings()
{
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::calculate_Yu_DRbar()
{
   assert(model && "CNE6SSMSusy_low_scale_constraint<Two_scale>::"
          "calculate_Yu_DRbar(): model pointer is zero");

   Eigen::Matrix<double,3,3> topDRbar(Eigen::Matrix<double,3,3>::Zero());
   topDRbar(0,0)      = oneset.displayMass(mUp);
   topDRbar(1,1)      = oneset.displayMass(mCharm);
   topDRbar(2,2)      = oneset.displayMass(mTop);

   // no thresholds

   const auto vu = MODELPARAMETER(vu);
   MODEL->set_Yu(Diag((1.4142135623730951*topDRbar)/vu));

}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::calculate_Yd_DRbar()
{
   assert(model && "CNE6SSMSusy_low_scale_constraint<Two_scale>::"
          "calculate_Yd_DRbar(): model pointer is zero");

   Eigen::Matrix<double,3,3> bottomDRbar(Eigen::Matrix<double,3,3>::Zero());
   bottomDRbar(0,0)   = oneset.displayMass(mDown);
   bottomDRbar(1,1)   = oneset.displayMass(mStrange);
   bottomDRbar(2,2)   = oneset.displayMass(mBottom);

   // no thresholds

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Yd(Diag((1.4142135623730951*bottomDRbar)/vd));

}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::calculate_Ye_DRbar()
{
   assert(model && "CNE6SSMSusy_low_scale_constraint<Two_scale>::"
          "calculate_Ye_DRbar(): model pointer is zero");

   Eigen::Matrix<double,3,3> electronDRbar(Eigen::Matrix<double,3,3>::Zero());
   electronDRbar(0,0) = oneset.displayMass(mElectron);
   electronDRbar(1,1) = oneset.displayMass(mMuon);
   electronDRbar(2,2) = oneset.displayMass(mTau);

   // no thresholds

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Ye(Diag((1.4142135623730951*electronDRbar)/vd));

}

void CNE6SSMSusy_low_scale_constraint<Two_scale>::calculate_MNeutrino_DRbar()
{
   assert(model && "CNE6SSMSusy_low_scale_constraint<Two_scale>::"
          "calculate_MNeutrino_DRbar(): model pointer is zero");

   neutrinoDRbar.setZero();
   neutrinoDRbar(0,0) = oneset.displayNeutrinoPoleMass(1);
   neutrinoDRbar(1,1) = oneset.displayNeutrinoPoleMass(2);
   neutrinoDRbar(2,2) = oneset.displayNeutrinoPoleMass(3);
}

/**
 * Recalculates the W boson pole mass using the new gauge couplings.
 */
void CNE6SSMSusy_low_scale_constraint<Two_scale>::recalculate_mw_pole()
{
   assert(model && "CNE6SSMSusy_low_scale_constraint<Two_scale>::"
          "recalculate_mw_pole(): model pointer is zero");

   if (!model->get_thresholds())
      return;


}

} // namespace flexiblesusy
