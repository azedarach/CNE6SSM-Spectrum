
#include "CNE6SSM_higgs_upper_bound.hpp"

#include "pv.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

CNE6SSM_higgs_upper_bound::CNE6SSM_higgs_upper_bound(const CNE6SSM_soft_parameters& model_)
   : model(model_)
   , upper_bound(0.)
   , include_all_gens(false)
   , include_ups(true)
   , include_exotics(true)
   , include_inert_singlets(true)
   , include_inert_neutral_higgs(true)
   , include_inert_charged_higgs(true)
{
}

CNE6SSM_higgs_upper_bound::~CNE6SSM_higgs_upper_bound()
{
}

double CNE6SSM_higgs_upper_bound::calculate_upper_bound()
{
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   const double v = Sqrt(Sqr(vd) + Sqr(vu));
   const double TanBeta = vu / vd;
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;
   const double Cos2Beta = Sqr(CosBeta) - Sqr(SinBeta);

   upper_bound = 0.5 * Sqr(Lambdax) * Sqr(v) * Sqr(Sin2Beta)
      + 0.25 * Sqr(Sqr(g2) + 0.6 * Sqr(g1)) * Sqr(v)
      * Sqr(Cos2Beta) + 0.025 * Sqr(g1p) * Sqr(v)
      * Sqr(-3.0 * Sqr(CosBeta) - 2.0 * Sqr(SinBeta));

   // additional contributions at 1-loop order
   if (model.get_loops() > 0) {
      upper_bound += Sqr(CosBeta) * get_tadpole_vd() / vd;
      upper_bound += Sqr(SinBeta) * get_tadpole_vu() / vu;
      if (include_ups) {
         if (include_all_gens) {
            upper_bound += get_up_contribution(0);
            upper_bound += get_up_contribution(1);
         }
         upper_bound += get_up_contribution(2);
      }
      if (include_exotics) {
         for (unsigned gen = 0; gen < 3; ++gen)
            upper_bound += get_exotic_contribution(gen);
      }
      if (include_inert_singlets) {
         for (unsigned gen = 0; gen < 3; ++gen)
            upper_bound += get_inert_singlet_contribution(gen);
      }
      if (include_inert_neutral_higgs) {
         for (unsigned gen = 0; gen < 4; ++gen)
            upper_bound += get_inert_neutral_higgs_contribution(gen);
      }
      if (include_inert_charged_higgs) {
         for (unsigned gen = 0; gen < 2; ++gen)
            upper_bound += get_inert_charged_higgs_contribution(gen);
      }
   }

   return upper_bound;
}

double CNE6SSM_higgs_upper_bound::get_tadpole_vd() const
{
   double result = 0.;

   if (include_ups) {
      if (include_all_gens) {
         result -= get_dV1lp_up_dvd(0);
         result -= get_dV1lp_up_dvd(1);
      }
      result -= get_dV1lp_up_dvd(2);
   }

   if (include_exotics) {
      for (unsigned gen = 0; gen < 3; ++gen)
         result -= get_dV1lp_exotic_dvd(gen);
   }

   if (include_inert_singlets) {
      for (unsigned gen = 0; gen < 3; ++gen)
         result -= get_dV1lp_inert_singlet_dvd(gen);
   }
   if (include_inert_neutral_higgs) {
      for (unsigned gen = 0; gen < 2; ++gen)
         result -= get_dV1lp_inert_neutral_higgs_dvd(gen);
   }
   if (include_inert_charged_higgs) {
      for (unsigned gen = 0; gen < 2; ++ gen)
         result -= get_dV1lp_inert_charged_higgs_dvd(gen);
   }

   return result;
}

double CNE6SSM_higgs_upper_bound::get_tadpole_vu() const
{
   double result = 0.;

   if (include_ups) {
      if (include_all_gens) {
         result -= get_dV1lp_up_dvu(0);
         result -= get_dV1lp_up_dvu(1);
      }
      result -= get_dV1lp_up_dvu(2);
   }

   if (include_exotics) {
      for (unsigned gen = 0; gen < 3; ++gen)
         result -= get_dV1lp_exotic_dvu(gen);
   }

   if (include_inert_singlets) {
      for (unsigned gen = 0; gen < 3; ++gen)
         result -= get_dV1lp_inert_singlet_dvu(gen);
   }
   if (include_inert_neutral_higgs) {
      for (unsigned gen = 0; gen < 2; ++ gen)
         result -= get_dV1lp_inert_neutral_higgs_dvu(gen);
   }
   if (include_inert_charged_higgs) {
      for (unsigned gen = 0; gen < 2; ++ gen)
         result -= get_dV1lp_inert_charged_higgs_dvu(gen);
   }

   return result;
}

double CNE6SSM_higgs_upper_bound::get_up_contribution(unsigned gen) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;

   const double Delta00Pr = get_unrotated_up_contribution(gen, 0, 0);
   const double Delta01Pr = get_unrotated_up_contribution(gen, 0, 1);
   const double Delta11Pr = get_unrotated_up_contribution(gen, 1, 1);

   const double result = Sqr(CosBeta) * Delta00Pr + Sqr(SinBeta)
      * Delta11Pr + Sin2Beta * Delta01Pr;

   return result;
}

double CNE6SSM_higgs_upper_bound::get_exotic_contribution(unsigned gen) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;

   const double Delta00Pr = get_unrotated_exotic_contribution(gen, 0, 0);
   const double Delta01Pr = get_unrotated_exotic_contribution(gen, 0, 1);
   const double Delta11Pr = get_unrotated_exotic_contribution(gen, 1, 1);

   const double result = Sqr(CosBeta) * Delta00Pr + Sqr(SinBeta)
      * Delta11Pr + Sin2Beta * Delta01Pr;

   return result;
}

double CNE6SSM_higgs_upper_bound::get_inert_singlet_contribution(unsigned gen) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;

   const double Delta00Pr = get_unrotated_inert_singlet_contribution(gen, 0, 0);
   const double Delta01Pr = get_unrotated_inert_singlet_contribution(gen, 0, 1);
   const double Delta11Pr = get_unrotated_inert_singlet_contribution(gen, 1, 1);

   const double result = Sqr(CosBeta) * Delta00Pr + Sqr(SinBeta)
      * Delta11Pr + Sin2Beta * Delta01Pr;

   return result;
}

double CNE6SSM_higgs_upper_bound::get_inert_neutral_higgs_contribution(unsigned gen) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;

   const double Delta00Pr = get_unrotated_inert_neutral_higgs_contribution(gen, 0, 0);
   const double Delta01Pr = get_unrotated_inert_neutral_higgs_contribution(gen, 0, 1);
   const double Delta11Pr = get_unrotated_inert_neutral_higgs_contribution(gen, 1, 1);

   const double result = Sqr(CosBeta) * Delta00Pr + Sqr(SinBeta)
      * Delta11Pr + Sin2Beta * Delta01Pr;

   return result;
}

double CNE6SSM_higgs_upper_bound::get_inert_charged_higgs_contribution(unsigned gen) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;

   const double Delta00Pr = get_unrotated_inert_charged_higgs_contribution(gen, 0, 0);
   const double Delta01Pr = get_unrotated_inert_charged_higgs_contribution(gen, 0, 1);
   const double Delta11Pr = get_unrotated_inert_charged_higgs_contribution(gen, 1, 1);

   const double result = Sqr(CosBeta) * Delta00Pr + Sqr(SinBeta)
      * Delta11Pr + Sin2Beta * Delta01Pr;

   return result;
}


double CNE6SSM_higgs_upper_bound::get_unrotated_up_contribution(unsigned gen, unsigned i, unsigned j) const
{
   double result = 0.;

   if (i == 0 && j == 0) {
      result = get_d2V1lp_up_dvd_dvd(gen);
   } else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
      result = get_d2V1lp_up_dvd_dvu(gen);
   } else if (i == 1 && j == 1) {
      result = get_d2V1lp_up_dvu_dvu(gen);
   }

   return result;
}

double CNE6SSM_higgs_upper_bound::get_unrotated_exotic_contribution(unsigned gen, unsigned i, unsigned j) const
{
   double result = 0.;

   if (i == 0 && j == 0) {
      result = get_d2V1lp_exotic_dvd_dvd(gen);
   } else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
      result = get_d2V1lp_exotic_dvd_dvu(gen);
   } else if (i == 1 && j == 1) {
      result = get_d2V1lp_exotic_dvu_dvu(gen);
   }

   return result;
}

double CNE6SSM_higgs_upper_bound::get_unrotated_inert_singlet_contribution(unsigned gen, unsigned i, unsigned j) const
{
   double result = 0.;

   if (i == 0 && j == 0) {
      result = get_d2V1lp_inert_singlet_dvd_dvd(gen);
   } else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
      result = get_d2V1lp_inert_singlet_dvd_dvu(gen);
   } else if (i == 1 && j == 1) {
      result = get_d2V1lp_inert_singlet_dvu_dvu(gen);
   }

   return result;
}

double CNE6SSM_higgs_upper_bound::get_unrotated_inert_neutral_higgs_contribution(unsigned gen, unsigned i, unsigned j) const
{
   double result = 0.;

   if (i == 0 && j == 0) {
      result = get_d2V1lp_inert_neutral_higgs_dvd_dvd(gen);
   } else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
      result = get_d2V1lp_inert_neutral_higgs_dvd_dvu(gen);
   } else if (i == 1 && j == 1) {
      result = get_d2V1lp_inert_neutral_higgs_dvu_dvu(gen);
   }

   return result;
}

double CNE6SSM_higgs_upper_bound::get_unrotated_inert_charged_higgs_contribution(unsigned gen, unsigned i, unsigned j) const
{
   double result = 0.;

   if (i == 0 && j == 0) {
      result = get_d2V1lp_inert_charged_higgs_dvd_dvd(gen);
   } else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
      result = get_d2V1lp_inert_charged_higgs_dvd_dvu(gen);
   } else if (i == 1 && j == 1) {
      result = get_d2V1lp_inert_charged_higgs_dvu_dvu(gen);
   }

   return result;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_mass_matrix_Su(unsigned gen) const
{
   const double mq2 = model.get_mq2(gen, gen);
   const double mu2 = model.get_mu2(gen, gen);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double Lambdax = model.get_Lambdax();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double yf = model.get_Yu(gen, gen);
   const double Tyf = model.get_TYu(gen, gen);
   const double QS = model.get_QS();

   const double mf2 = calculate_MFu2(gen);

   const double DeltaQ = 0.0125 * Sqr(g1p) * (-3.0 * Sqr(vd)
      - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));
   const double DeltaU = 0.0125 * Sqr(g1p) * (-3.0 * Sqr(vd)
      - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));

   Eigen::Matrix<double,2,2> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2 + 0.125 * (Sqr(g2) - 0.2 *
      Sqr(g1)) * (Sqr(vd) - Sqr(vu)) + DeltaQ + mf2;
   mass_matrix_Su(0,1) = Tyf * vu / Sqrt(2.0) - 0.5 *
      Lambdax * yf * vd * vs;
   mass_matrix_Su(1,0) = mass_matrix_Su(0,1);
   mass_matrix_Su(1,1) = mu2 + 0.1 * Sqr(g1) * (Sqr(vd)
      - Sqr(vu)) + DeltaU + mf2;

   return mass_matrix_Su;
}

Eigen::Array<double,2,1> CNE6SSM_higgs_upper_bound::calculate_MSu2(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_Su(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   const double mass_sum = mass_matrix(0,0) + mass_matrix(1,1);

   Eigen::Array<double,2,1> MSu2;

   MSu2(0) = 0.5 * (mass_sum - mass_diff);
   MSu2(1) = 0.5 * (mass_sum + mass_diff);

   return MSu2;
}

double CNE6SSM_higgs_upper_bound::calculate_Sin2ThetaSu(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_Su(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return 2.0 * mass_matrix(0,1) / mass_diff;
}

double CNE6SSM_higgs_upper_bound::calculate_Cos2ThetaSu(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_Su(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return (mass_matrix(1,1) - mass_matrix(0,0)) / mass_diff;
}

double CNE6SSM_higgs_upper_bound::calculate_MFu2(unsigned gen) const
{
   const double yf = model.get_Yu(gen, gen);
   const double vev = model.get_vu();

   return 0.5 * Sqr(yf * vev);
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_mass_matrix_SDX(unsigned gen) const
{
   const double mDx2 = model.get_mDx2(gen, gen);
   const double mDxbar2 = model.get_mDxbar2(gen, gen);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double Lambdax = model.get_Lambdax();
   const double Sigmax = model.get_Sigmax();
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();
   const double Kappa = model.get_Kappa(gen, gen);
   const double TKappa = model.get_TKappa(gen, gen);
   const double QS = model.get_QS();

   const double mf2 = calculate_MFDX2(gen);

   const double DeltaDx = -0.025 * Sqr(g1p) * (-3.0 * Sqr(vd)
      - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));
   const double DeltaDxbar = -0.0375 * Sqr(g1p) * (-3.0 *
      Sqr(vd) - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));

   Eigen::Matrix<double,2,2> mass_matrix_SDX;

   mass_matrix_SDX(0,0) = mDx2 + 0.05 * Sqr(g1) * (Sqr(vd)
      - Sqr(vu)) + DeltaDx + mf2;
   mass_matrix_SDX(0,1) = TKappa * vs / Sqrt(2.0) - 0.5 *
      Kappa * (Lambdax * vd * vu + Sigmax * vphi * vsb);
   mass_matrix_SDX(1,0) = mass_matrix_SDX(0,1);
   mass_matrix_SDX(1,1) = mDxbar2 - 0.05 * Sqr(g1) * (Sqr(vd)
      - Sqr(vu)) + DeltaDxbar + mf2;

   return mass_matrix_SDX;
}

Eigen::Array<double,2,1> CNE6SSM_higgs_upper_bound::calculate_MSDX2(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_SDX(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   const double mass_sum = mass_matrix(0,0) + mass_matrix(1,1);

   Eigen::Array<double,2,1> MSDX2;

   MSDX2(0) = 0.5 * (mass_sum - mass_diff);
   MSDX2(1) = 0.5 * (mass_sum + mass_diff);

   return MSDX2;
}

double CNE6SSM_higgs_upper_bound::calculate_Sin2ThetaSDX(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_SDX(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return 2.0 * mass_matrix(0,1) / mass_diff;
}

double CNE6SSM_higgs_upper_bound::calculate_Cos2ThetaSDX(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_SDX(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return (mass_matrix(1,1) - mass_matrix(0,0)) / mass_diff;
}

double CNE6SSM_higgs_upper_bound::calculate_MFDX2(unsigned gen) const
{
   const double yf = model.get_Kappa(gen, gen);
   const double vev = model.get_vs();

   return 0.5 * Sqr(yf * vev);
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_mass_matrix_HI0(unsigned gen) const
{
   const double mH1I2 = model.get_mH1I2(gen, gen);
   const double mH2I2 = model.get_mH2I2(gen, gen);
   const double Lambda = model.get_Lambda12(gen, gen);
   const double TLambda = model.get_TLambda12(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double Sigmax = model.get_Sigmax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();
   const double mf2 = calculate_MFHI02(gen);

   const double DeltaHd = -0.0375 * Sqr(g1p) * (-3.0 * Sqr(vd)
      - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));
   const double DeltaHu = -0.025 * Sqr(g1p) * (-3.0 *
      Sqr(vd) - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));

   Eigen::Matrix<double,2,2> mass_matrix_HI0;

   mass_matrix_HI0(0,0) = mH1I2 + 0.125 * (Sqr(g2) + 0.6 * Sqr(g1))
      * (Sqr(vd) - Sqr(vu)) + DeltaHd + mf2;
   mass_matrix_HI0(0,1) = -TLambda * vs / Sqrt(2.0) + 0.5 * Lambda
      * (Lambdax * vd * vu + Sigmax * vphi * vsb);
   mass_matrix_HI0(1,0) = mass_matrix_HI0(0,1);
   mass_matrix_HI0(1,1) = mH2I2 - 0.125 * (Sqr(g2) + 0.6 * Sqr(g1))
      * (Sqr(vd) - Sqr(vu)) + DeltaHu + mf2;

   return mass_matrix_HI0;
}

Eigen::Array<double,2,1> CNE6SSM_higgs_upper_bound::calculate_MHI02(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_HI0(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   const double mass_sum = mass_matrix(0,0) + mass_matrix(1,1);

   Eigen::Array<double,2,1> MHI02;

   MHI02(0) = 0.5 * (mass_sum - mass_diff);
   MHI02(1) = 0.5 * (mass_sum + mass_diff);

   return MHI02;
}

double CNE6SSM_higgs_upper_bound::calculate_Sin2ThetaHI0(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_HI0(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return 2.0 * mass_matrix(0,1) / mass_diff;
}

double CNE6SSM_higgs_upper_bound::calculate_Cos2ThetaHI0(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_HI0(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return (mass_matrix(1,1) - mass_matrix(0,0)) / mass_diff;
}

double CNE6SSM_higgs_upper_bound::calculate_MFHI02(unsigned gen) const
{
   const double Lambda = model.get_Lambda12(gen, gen);
   const double vev = model.get_vs();

   return 0.5 * Sqr(Lambda * vev);
}

double CNE6SSM_higgs_upper_bound::calculate_MSI02(unsigned gen) const
{
   const double mSI2 = model.get_mSI2(gen, gen);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();

   const double DeltaS = 0.0625 * Sqr(g1p) * (-3.0 *
      Sqr(vd) - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));

   const double MSI02 = mSI2 + DeltaS;

   return MSI02;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_mass_matrix_HIPM(unsigned gen) const
{
   const double mH1I2 = model.get_mH1I2(gen, gen);
   const double mH2I2 = model.get_mH2I2(gen, gen);
   const double Lambda = model.get_Lambda12(gen, gen);
   const double TLambda = model.get_TLambda12(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double Sigmax = model.get_Sigmax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();
   const double mf2 = calculate_MFHI02(gen);

   const double DeltaHd = -0.0375 * Sqr(g1p) * (-3.0 * Sqr(vd)
      - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));
   const double DeltaHu = -0.025 * Sqr(g1p) * (-3.0 *
      Sqr(vd) - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));

   Eigen::Matrix<double,2,2> mass_matrix_HIPM;

   mass_matrix_HIPM(0,0) = mH1I2 - 0.125 * (Sqr(g2) - 0.6 * Sqr(g1))
      * (Sqr(vd) - Sqr(vu)) + DeltaHd + mf2;
   mass_matrix_HIPM(0,1) = TLambda * vs / Sqrt(2.0) - 0.5 * Lambda
      * (Lambdax * vd * vu + Sigmax * vphi * vsb);
   mass_matrix_HIPM(1,0) = mass_matrix_HIPM(0,1);
   mass_matrix_HIPM(1,1) = mH2I2 + 0.125 * (Sqr(g2) - 0.6 * Sqr(g1))
      * (Sqr(vd) - Sqr(vu)) + DeltaHu + mf2;

   return mass_matrix_HIPM;
}

Eigen::Array<double,2,1> CNE6SSM_higgs_upper_bound::calculate_MHIPM2(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_HIPM(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   const double mass_sum = mass_matrix(0,0) + mass_matrix(1,1);

   Eigen::Array<double,2,1> MHIPM2;

   MHIPM2(0) = 0.5 * (mass_sum - mass_diff);
   MHIPM2(1) = 0.5 * (mass_sum + mass_diff);

   return MHIPM2;
}

double CNE6SSM_higgs_upper_bound::calculate_Sin2ThetaHIPM(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_HIPM(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return 2.0 * mass_matrix(0,1) / mass_diff;
}

double CNE6SSM_higgs_upper_bound::calculate_Cos2ThetaHIPM(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_HIPM(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return (mass_matrix(1,1) - mass_matrix(0,0)) / mass_diff;
}

double CNE6SSM_higgs_upper_bound::calculate_MFHIPM2(unsigned gen) const
{
   const double Lambda = model.get_Lambda12(gen, gen);
   const double vev = model.get_vs();

   return 0.5 * Sqr(Lambda * vev);
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_dmass_matrix_Su_dvd(unsigned gen) const
{
   const double yf = model.get_Yu(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vs = model.get_vs();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) * vd - 0.075 * Sqr(g1p) * vd;
   derivatives(0,1) = -0.5 * Lambdax * yf * vs;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.2 * Sqr(g1) * vd - 0.075 * Sqr(g1p) * vd;

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_dmass_matrix_Su_dvu(unsigned gen) const
{
   const double yf = model.get_Yu(gen, gen);
   const double Tyf = model.get_TYu(gen, gen);
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) * vu - 0.05 * Sqr(g1p) * vu
      + Sqr(yf) * vu;
   derivatives(0,1) = Tyf / Sqrt(2.0);
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.2 * Sqr(g1) * vu - 0.05 * Sqr(g1p) * vu + Sqr(yf) * vu;

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_dmass_matrix_SDX_dvd(unsigned gen) const
{
   const double Lambdax = model.get_Lambdax();
   const double Kappa = model.get_Kappa(gen, gen);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.1 * Sqr(g1) * vd + 0.15 * Sqr(g1p) * vd;
   derivatives(0,1) = -0.5 * Kappa * Lambdax * vu;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.1 * Sqr(g1) * vd + 0.225 * Sqr(g1p) * vd;

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_dmass_matrix_SDX_dvu(unsigned gen) const
{
   const double Lambdax = model.get_Lambdax();
   const double Kappa = model.get_Kappa(gen, gen);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.1 * Sqr(g1) * vu + 0.1 * Sqr(g1p) * vu;
   derivatives(0,1) = -0.5 * Kappa * Lambdax * vd;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.1 * Sqr(g1) * vu + 0.15 * Sqr(g1p) * vu;

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_dmass_matrix_HI0_dvd(unsigned gen) const
{
   const double Lambda = model.get_Lambda12(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) * vd + 0.225 * Sqr(g1p) * vd;
   derivatives(0,1) = 0.5 * Lambda * Lambdax * vu;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) * vd + 0.15 * Sqr(g1p) * vd;

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_dmass_matrix_HI0_dvu(unsigned gen) const
{
   const double Lambda = model.get_Lambda12(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) * vu + 0.15 * Sqr(g1p) * vu;
   derivatives(0,1) = 0.5 * Lambda * Lambdax * vd;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) * vu + 0.1 * Sqr(g1p) * vu;

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_dmass_matrix_HIPM_dvd(unsigned gen) const
{
   const double Lambda = model.get_Lambda12(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) * vd + 0.225 * Sqr(g1p) * vd;
   derivatives(0,1) = -0.5 * Lambda * Lambdax * vu;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) * vd + 0.15 * Sqr(g1p) * vd;

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_dmass_matrix_HIPM_dvu(unsigned gen) const
{
   const double Lambda = model.get_Lambda12(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) * vu + 0.15 * Sqr(g1p) * vu;
   derivatives(0,1) = -0.5 * Lambda * Lambdax * vd;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) * vu + 0.1 * Sqr(g1p) * vu;

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_d2mass_matrix_Su_dvd_dvd(unsigned) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) - 0.075 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.2 * Sqr(g1) - 0.075 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_d2mass_matrix_Su_dvd_dvu(unsigned) const
{
   return Eigen::Matrix<double,2,2>::Zero();
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_d2mass_matrix_Su_dvu_dvu(unsigned gen) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double yf = model.get_Yu(gen, gen);

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) - 0.05 * Sqr(g1p) + Sqr(yf);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.2 * Sqr(g1) - 0.05 * Sqr(g1p) + Sqr(yf);

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_d2mass_matrix_SDX_dvd_dvd(unsigned) const
{
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.1 * Sqr(g1) + 0.15 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.1 * Sqr(g1) + 0.225 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_d2mass_matrix_SDX_dvd_dvu(unsigned gen) const
{
   const double Lambdax = model.get_Lambdax();
   const double Kappa = model.get_Kappa(gen, gen);

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.;
   derivatives(0,1) = -0.5 * Kappa * Lambdax;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.;

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_d2mass_matrix_SDX_dvu_dvu(unsigned) const
{
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.1 * Sqr(g1) + 0.1 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.1 * Sqr(g1) + 0.15 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_d2mass_matrix_HI0_dvd_dvd(unsigned) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) + 0.225 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) + 0.15 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_d2mass_matrix_HI0_dvd_dvu(unsigned gen) const
{
   const double Lambdax = model.get_Lambdax();
   const double Lambda = model.get_Lambda12(gen, gen);

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.;
   derivatives(0,1) = 0.5 * Lambda * Lambdax;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.;

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_d2mass_matrix_HI0_dvu_dvu(unsigned) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) + 0.15 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) + 0.1 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_d2mass_matrix_HIPM_dvd_dvd(unsigned) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) + 0.225 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) + 0.15 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_d2mass_matrix_HIPM_dvd_dvu(unsigned gen) const
{
   const double Lambdax = model.get_Lambdax();
   const double Lambda = model.get_Lambda12(gen, gen);

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.;
   derivatives(0,1) = -0.5 * Lambda * Lambdax;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.;

   return derivatives;
}

Eigen::Matrix<double,2,2> CNE6SSM_higgs_upper_bound::get_d2mass_matrix_HIPM_dvu_dvu(unsigned) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) + 0.15 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) + 0.1 * Sqr(g1p);

   return derivatives;
}

double CNE6SSM_higgs_upper_bound::get_dV1lp_up_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_Su_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaSu(gen);
   const double Cos2Theta = calculate_Cos2ThetaSu(gen);

   const double dMS20_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSu2(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * 3.0 * (dMS20_dvd * A0MS0
      + dMS21_dvd * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_dV1lp_up_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_Su_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSu(gen);
   const double Cos2Theta = calculate_Cos2ThetaSu(gen);

   const double dMS20_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));
   const double yf = model.get_Yu(gen, gen);
   const double vu = model.get_vu();
   const double dMF2_dvu = Sqr(yf) * vu;

   const Eigen::Array<double,2,1> MS2(calculate_MSu2(gen));
   const double MF2 = calculate_MFu2(gen);
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));
   const double A0MF = passarino_veltman::ReA0(MF2, Sqr(scale));

   const double result = -oneOver16PiSqr * 3.0 * (dMS20_dvu * A0MS0
      + dMS21_dvu * A0MS1 - 2.0 * dMF2_dvu * A0MF);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_dV1lp_exotic_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_SDX_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaSDX(gen);
   const double Cos2Theta = calculate_Cos2ThetaSDX(gen);

   const double dMS20_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSDX2(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * 3.0 * (dMS20_dvd * A0MS0
      + dMS21_dvd * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_dV1lp_exotic_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_SDX_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSDX(gen);
   const double Cos2Theta = calculate_Cos2ThetaSDX(gen);

   const double dMS20_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSDX2(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * 3.0 * (dMS20_dvu * A0MS0
      + dMS21_dvu * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_dV1lp_inert_singlet_dvd(unsigned gen) const
{
   const double vd = model.get_vd();
   const double g1p = model.get_g1p();
   const double dMS2_dvd = -0.375 * Sqr(g1p) * vd;
   
   const double MS2 = calculate_MSI02(gen);
   const double scale = model.get_scale();
   const double A0MS = passarino_veltman::ReA0(MS2, Sqr(scale));
   
   const double result = -oneOver16PiSqr * dMS2_dvd * A0MS;
   
   return result;
}

double CNE6SSM_higgs_upper_bound::get_dV1lp_inert_singlet_dvu(unsigned gen) const
{
   const double vu = model.get_vu();
   const double g1p = model.get_g1p();
   const double dMS2_dvu = -0.25 * Sqr(g1p) * vu;
   
   const double MS2 = calculate_MSI02(gen);
   const double scale = model.get_scale();
   const double A0MS = passarino_veltman::ReA0(MS2, Sqr(scale));
   
   const double result = -oneOver16PiSqr * dMS2_dvu * A0MS;
   
   return result;
}

double CNE6SSM_higgs_upper_bound::get_dV1lp_inert_neutral_higgs_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_HI0_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaHI0(gen);
   const double Cos2Theta = calculate_Cos2ThetaHI0(gen);

   const double dMS20_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHI02(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * (dMS20_dvd * A0MS0
      + dMS21_dvd * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_dV1lp_inert_neutral_higgs_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_HI0_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaHI0(gen);
   const double Cos2Theta = calculate_Cos2ThetaHI0(gen);

   const double dMS20_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHI02(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * (dMS20_dvu * A0MS0
      + dMS21_dvu * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_dV1lp_inert_charged_higgs_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_HIPM_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaHIPM(gen);
   const double Cos2Theta = calculate_Cos2ThetaHIPM(gen);

   const double dMS20_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHIPM2(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * (dMS20_dvd * A0MS0
      + dMS21_dvd * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_dV1lp_inert_charged_higgs_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_HIPM_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaHIPM(gen);
   const double Cos2Theta = calculate_Cos2ThetaHIPM(gen);

   const double dMS20_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHIPM2(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * (dMS20_dvu * A0MS0
      + dMS21_dvu * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_up_dvd_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_Su_dvd(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_Su_dvd_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaSu(gen);
   const double Cos2Theta = calculate_Cos2ThetaSu(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSu2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (Sqr(dMS20_dvd) * logMS20Q2
      - d2MS20_dvd_dvd * A0MS0 + Sqr(dMS21_dvd) * logMS21Q2 - d2MS21_dvd_dvd
      * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_up_dvd_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> vd_first_derivatives(get_dmass_matrix_Su_dvd(gen));
   const Eigen::Matrix<double,2,2> vu_first_derivatives(get_dmass_matrix_Su_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_Su_dvd_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSu(gen);
   const double Cos2Theta = calculate_Cos2ThetaSu(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) + Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) - 2.0 * Sin2Theta * vd_first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) - Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) + 2.0 * Sin2Theta * vd_first_derivatives(0,1));

   const double dMS20_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) + Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) - 2.0 * Sin2Theta * vu_first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) - Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) + 2.0 * Sin2Theta * vu_first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSu2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) + Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) - 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double d2MS21_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) - Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) + 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (dMS20_dvd * dMS20_dvu *
      logMS20Q2 - d2MS20_dvd_dvu * A0MS0 + dMS21_dvd * dMS21_dvu * logMS21Q2
      - d2MS21_dvd_dvu * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_up_dvu_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_Su_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_Su_dvu_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSu(gen);
   const double Cos2Theta = calculate_Cos2ThetaSu(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double yf = model.get_Yu(gen, gen);
   const double vu = model.get_vu();

   const double dMS20_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMF2_dvu = Sqr(yf) * vu;

   const Eigen::Array<double,2,1> MS2(calculate_MSu2(gen));
   const double MF2 = calculate_MFu2(gen);
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MF2_dvu_dvu = Sqr(yf);

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double logMF2Q2 = Log(MF2 / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));
   const double A0MF = passarino_veltman::ReA0(MF2, Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (Sqr(dMS20_dvu) * logMS20Q2
      - d2MS20_dvu_dvu * A0MS0 + Sqr(dMS21_dvu) * logMS21Q2 - d2MS21_dvu_dvu
      * A0MS1 - 2.0 * Sqr(dMF2_dvu) * logMF2Q2 + 2.0 * d2MF2_dvu_dvu * A0MF);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_exotic_dvd_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_SDX_dvd(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_SDX_dvd_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaSDX(gen);
   const double Cos2Theta = calculate_Cos2ThetaSDX(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSDX2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (Sqr(dMS20_dvd) * logMS20Q2
      - d2MS20_dvd_dvd * A0MS0 + Sqr(dMS21_dvd) * logMS21Q2 - d2MS21_dvd_dvd
      * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_exotic_dvd_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> vd_first_derivatives(get_dmass_matrix_SDX_dvd(gen));
   const Eigen::Matrix<double,2,2> vu_first_derivatives(get_dmass_matrix_SDX_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_SDX_dvd_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSDX(gen);
   const double Cos2Theta = calculate_Cos2ThetaSDX(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) + Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) - 2.0 * Sin2Theta * vd_first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) - Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) + 2.0 * Sin2Theta * vd_first_derivatives(0,1));

   const double dMS20_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) + Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) - 2.0 * Sin2Theta * vu_first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) - Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) + 2.0 * Sin2Theta * vu_first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSDX2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) + Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) - 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double d2MS21_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) - Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) + 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (dMS20_dvd * dMS20_dvu *
      logMS20Q2 - d2MS20_dvd_dvu * A0MS0 + dMS21_dvd * dMS21_dvu * logMS21Q2
      - d2MS21_dvd_dvu * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_exotic_dvu_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_SDX_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_SDX_dvu_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSDX(gen);
   const double Cos2Theta = calculate_Cos2ThetaSDX(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSDX2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (Sqr(dMS20_dvu) * logMS20Q2
      - d2MS20_dvu_dvu * A0MS0 + Sqr(dMS21_dvu) * logMS21Q2 - d2MS21_dvu_dvu
      * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_inert_singlet_dvd_dvd(unsigned gen) const
{
   const double g1p = model.get_g1p();
   const double vd = model.get_vd();

   const double dMS2_dvd = -0.375 * Sqr(g1p) * vd;
   const double d2MS2_dvd_dvd = -0.375 * Sqr(g1p);

   const double MS2 = calculate_MSI02(gen);

   const double scale = model.get_scale();
   const double logMS2Q2 = Log(MS2 / Sqr(scale));
   const double A0MS = passarino_veltman::ReA0(MS2, Sqr(scale));

   const double result = oneOver16PiSqr * (Sqr(dMS2_dvd) * logMS2Q2
      - d2MS2_dvd_dvd * A0MS);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_inert_singlet_dvd_dvu(unsigned gen) const
{
   const double g1p = model.get_g1p();
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double dMS2_dvd = - 0.375 * Sqr(g1p) * vd;
   const double dMS2_dvu = -0.25 * Sqr(g1p) * vu;

   const double MS2 = calculate_MSI02(gen);

   const double scale = model.get_scale();
   const double logMS2Q2 = Log(MS2 / Sqr(scale));

   const double result = oneOver16PiSqr * (dMS2_dvd * dMS2_dvu *
      logMS2Q2);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_inert_singlet_dvu_dvu(unsigned gen) const
{
   const double g1p = model.get_g1p();
   const double vu = model.get_vu();

   const double dMS2_dvu = -0.25 * Sqr(g1p) * vu;
   const double d2MS2_dvu_dvu = -0.25 * Sqr(g1p);

   const double MS2 = calculate_MSI02(gen);

   const double scale = model.get_scale();
   const double logMS2Q2 = Log(MS2 / Sqr(scale));
   const double A0MS = passarino_veltman::ReA0(MS2, Sqr(scale));

   const double result = oneOver16PiSqr * (Sqr(dMS2_dvu) * logMS2Q2
      - d2MS2_dvu_dvu * A0MS);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_inert_neutral_higgs_dvd_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_HI0_dvd(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_HI0_dvd_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaHI0(gen);
   const double Cos2Theta = calculate_Cos2ThetaHI0(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHI02(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * (Sqr(dMS20_dvd) * logMS20Q2
      - d2MS20_dvd_dvd * A0MS0 + Sqr(dMS21_dvd) * logMS21Q2 - d2MS21_dvd_dvd
      * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_inert_neutral_higgs_dvd_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> vd_first_derivatives(get_dmass_matrix_HI0_dvd(gen));
   const Eigen::Matrix<double,2,2> vu_first_derivatives(get_dmass_matrix_HI0_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_HI0_dvd_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaHI0(gen);
   const double Cos2Theta = calculate_Cos2ThetaHI0(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) + Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) - 2.0 * Sin2Theta * vd_first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) - Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) + 2.0 * Sin2Theta * vd_first_derivatives(0,1));

   const double dMS20_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) + Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) - 2.0 * Sin2Theta * vu_first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) - Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) + 2.0 * Sin2Theta * vu_first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHI02(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) + Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) - 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double d2MS21_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) - Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) + 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * (dMS20_dvd * dMS20_dvu *
      logMS20Q2 - d2MS20_dvd_dvu * A0MS0 + dMS21_dvd * dMS21_dvu * logMS21Q2
      - d2MS21_dvd_dvu * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_inert_neutral_higgs_dvu_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_HI0_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_HI0_dvu_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaHI0(gen);
   const double Cos2Theta = calculate_Cos2ThetaHI0(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHI02(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * (Sqr(dMS20_dvu) * logMS20Q2
      - d2MS20_dvu_dvu * A0MS0 + Sqr(dMS21_dvu) * logMS21Q2 - d2MS21_dvu_dvu
      * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_inert_charged_higgs_dvd_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_HIPM_dvd(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_HIPM_dvd_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaHIPM(gen);
   const double Cos2Theta = calculate_Cos2ThetaHIPM(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHIPM2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * (Sqr(dMS20_dvd) * logMS20Q2
      - d2MS20_dvd_dvd * A0MS0 + Sqr(dMS21_dvd) * logMS21Q2 - d2MS21_dvd_dvd
      * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_inert_charged_higgs_dvd_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> vd_first_derivatives(get_dmass_matrix_HIPM_dvd(gen));
   const Eigen::Matrix<double,2,2> vu_first_derivatives(get_dmass_matrix_HIPM_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_HIPM_dvd_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaHIPM(gen);
   const double Cos2Theta = calculate_Cos2ThetaHIPM(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) + Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) - 2.0 * Sin2Theta * vd_first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) - Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) + 2.0 * Sin2Theta * vd_first_derivatives(0,1));

   const double dMS20_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) + Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) - 2.0 * Sin2Theta * vu_first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) - Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) + 2.0 * Sin2Theta * vu_first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHIPM2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) + Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) - 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double d2MS21_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) - Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) + 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * (dMS20_dvd * dMS20_dvu *
      logMS20Q2 - d2MS20_dvd_dvu * A0MS0 + dMS21_dvd * dMS21_dvu * logMS21Q2
      - d2MS21_dvd_dvu * A0MS1);

   return result;
}

double CNE6SSM_higgs_upper_bound::get_d2V1lp_inert_charged_higgs_dvu_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_HIPM_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_HIPM_dvu_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaHIPM(gen);
   const double Cos2Theta = calculate_Cos2ThetaHIPM(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHIPM2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * (Sqr(dMS20_dvu) * logMS20Q2
      - d2MS20_dvu_dvu * A0MS0 + Sqr(dMS21_dvu) * logMS21Q2 - d2MS21_dvu_dvu
      * A0MS1);

   return result;
}

} // namespace flexiblesusy
