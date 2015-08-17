// ====================================================================
// Test suite for implementation of approximate 1-loop upper bound
// on the Higgs mass in the CNE6SSM
// ====================================================================

#include <complex>
#include <iomanip>
#include <limits>
#include <random>
#include <chrono>
#include <Eigen/Core>
#include <gsl/gsl_deriv.h>

#include "CNE6SSM_higgs_upper_bound.hpp"
#include "CNE6SSM_mass_eigenstates.hpp"

#include "pv.hpp"
#include "wrappers.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CNE6SSM_higgs_upper_bound

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

namespace flexiblesusy {

struct A0_fn {
   A0_fn(double scale_) : scale(scale_) {}
   double operator()(double m) const {
      return passarino_veltman::ReA0(m*m, Sqr(scale));
   }
   double scale;
};

std::complex<double> partial_tadpole_hh(const CNE6SSM_mass_eigenstates& model, 
                                        unsigned gO1)
{
   Eigen::Array<double,6,1> MSu(model.get_MSu());
   Eigen::Array<double,3,1> MFu(model.get_MFu());
   Eigen::Array<double,6,1> MSd(model.get_MSd());
   Eigen::Array<double,3,1> MFd(model.get_MFd());
   Eigen::Array<double,6,1> MSDX(model.get_MSDX());
   Eigen::Array<double,3,1> MFDX(model.get_MFDX());
   Eigen::Array<double,7,1> MSHI0(model.get_MSHI0());
   Eigen::Array<double,7,1> MChiI(model.get_MChiI());
   Eigen::Array<double,4,1> MSHIPM(model.get_MSHIPM());
   Eigen::Array<double,2,1> MChaI(model.get_MChaI());

   A0_fn A0(model.get_scale());

   std::complex<double> result;

   std::complex<double> tmp_9736;
   std::complex<double> tmp_9737;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_9737 += A0(MChaI(gI1))*(model.CpUhhbarChaIChaIPL(gO1,gI1,gI1) +
         model.CpUhhbarChaIChaIPR(gO1,gI1,gI1))*MChaI(gI1);
   }
   tmp_9736 += tmp_9737;
   result += (2) * tmp_9736;
   std::complex<double> tmp_9743;
   std::complex<double> tmp_9744;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_9744 += A0(MFDX(gI1))*(model.CpUhhbarFDXFDXPL(gO1,gI1,gI1) +
         model.CpUhhbarFDXFDXPR(gO1,gI1,gI1))*MFDX(gI1);
   }
   tmp_9743 += tmp_9744;
   result += (6) * tmp_9743;
   std::complex<double> tmp_9747;
   std::complex<double> tmp_9748;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_9748 += A0(MFu(gI1))*(model.CpUhhbarFuFuPL(gO1,gI1,gI1)
         + model.CpUhhbarFuFuPR(gO1,gI1,gI1))*MFu(gI1);
   }
   tmp_9747 += tmp_9748;
   result += (6) * tmp_9747;
   std::complex<double> tmp_9749;
   std::complex<double> tmp_9750;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_9750 += A0(MSHIPM(gI1))*model.CpUhhconjSHIPMSHIPM(gO1,gI1,gI1);
   }
   tmp_9749 += tmp_9750;
   result += (-1) * tmp_9749;
   std::complex<double> tmp_9757;
   std::complex<double> tmp_9758;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_9758 += A0(MSDX(gI1))*model.CpUhhconjSDXSDX(gO1,gI1,gI1);
   }
   tmp_9757 += tmp_9758;
   result += (-3) * tmp_9757;
   std::complex<double> tmp_9761;
   std::complex<double> tmp_9762;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_9762 += A0(MSu(gI1))*model.CpUhhconjSuSu(gO1,gI1,gI1);
   }
   tmp_9761 += tmp_9762;
   result += (-3) * tmp_9761;
   std::complex<double> tmp_9763;
   std::complex<double> tmp_9764;
   for (unsigned gI1 = 0; gI1 < 7; ++gI1) {
      tmp_9764 += A0(MSHI0(gI1))*model.CpUhhconjSHI0SHI0(gO1,gI1,gI1);
   }
   tmp_9763 += tmp_9764;
   result += (-1) * tmp_9763;
   std::complex<double> tmp_9765;
   for (unsigned gI1 = 0; gI1 < 7; ++gI1) {
      tmp_9765 += A0(MChiI(gI1))*(model.CpUhhChiIChiIPL(gO1,gI1,gI1) +
         model.CpUhhChiIChiIPR(gO1,gI1,gI1))*MChiI(gI1);
   }
   result += tmp_9765;
   std::complex<double> tmp_9741;
   std::complex<double> tmp_9742;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_9742 += A0(MFd(gI1))*(model.CpUhhbarFdFdPL(gO1,gI1,gI1)
         + model.CpUhhbarFdFdPR(gO1,gI1,gI1))*MFd(gI1);
   }
   tmp_9741 += tmp_9742;
   result += (6) * tmp_9741;
   std::complex<double> tmp_9755;
   std::complex<double> tmp_9756;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_9756 += A0(MSd(gI1))*model.CpUhhconjSdSd(gO1,gI1,gI1);
   }
   tmp_9755 += tmp_9756;
   result += (-3) * tmp_9755;

   return result * oneOver16PiSqr;
}

} // namespace flexiblesusy

using namespace flexiblesusy;

void set_test_model_parameters(CNE6SSM_mass_eigenstates& model)
{
   model.set_scale(1000.0);//3864.263);
   model.set_QS(5.0);

   model.set_g1(0.4741);
   model.set_g2(0.6358);
   model.set_g3(1.0649);
   model.set_g1p(0.4550);

   model.set_Yu(0,0,7.007e-6);
   model.set_Yu(0,1,0.);
   model.set_Yu(0,2,0.);
   model.set_Yu(1,0,0.);
   model.set_Yu(1,1,3.203e-3);
   model.set_Yu(1,2,0.);
   model.set_Yu(2,0,0.);
   model.set_Yu(2,1,0.);
   model.set_Yu(2,2,0.8287);

   model.set_Yd(0,0,1.323e-4);
   model.set_Yd(0,1,0.);
   model.set_Yd(0,2,0.);
   model.set_Yd(1,0,0.);
   model.set_Yd(1,1,2.896e-3);
   model.set_Yd(1,2,0.);
   model.set_Yd(2,0,0.);
   model.set_Yd(2,1,0.);
   model.set_Yd(2,2,0.1273);

   model.set_Ye(0,0,2.818e-5);
   model.set_Ye(0,1,0.);
   model.set_Ye(0,2,0.);
   model.set_Ye(1,0,0.);
   model.set_Ye(1,1,5.826e-3);
   model.set_Ye(1,2,0.);
   model.set_Ye(2,0,0.);
   model.set_Ye(2,1,0.);
   model.set_Ye(2,2,0.0982);

   model.set_TYe(0,0,-0.03741);
   model.set_TYe(0,1,0.);
   model.set_TYe(0,2,0.);
   model.set_TYe(1,0,0.);
   model.set_TYe(1,1,-7.735);
   model.set_TYe(1,2,0.);
   model.set_TYe(2,0,0.);
   model.set_TYe(2,1,0.);
   model.set_TYe(2,2,-129.651);

   model.set_TYd(0,0,-0.7829);
   model.set_TYd(0,1,0.);
   model.set_TYd(0,2,0.);
   model.set_TYd(1,0,0.);
   model.set_TYd(1,1,-17.1405);
   model.set_TYd(1,2,0.);
   model.set_TYd(2,0,0.);
   model.set_TYd(2,1,0.);
   model.set_TYd(2,2,-703.529);

   model.set_TYu(0,0,-0.03342);
   model.set_TYu(0,1,0.);
   model.set_TYu(0,2,0.);
   model.set_TYu(1,0,0.);
   model.set_TYu(1,1,-15.275);
   model.set_TYu(1,2,0.);
   model.set_TYu(2,0,0.);
   model.set_TYu(2,1,0.);
   model.set_TYu(2,2,-3026.593);

   model.set_XiF(2.08775906e8);
   model.set_MuPhi(0.);
   model.set_KappaPr(0.01712);
   model.set_Sigmax(0.1132);
   model.set_LXiF(3.97264570e11);
   model.set_vd(25.324);
   model.set_vu(240.071);

   model.set_hE(0,0,0.);
   model.set_hE(0,1,0.);
   model.set_hE(1,0,0.);
   model.set_hE(1,1,0.);
   model.set_hE(2,0,0.);
   model.set_hE(2,1,0.);

   model.set_SigmaL(0.4102);
   model.set_TKappaPr(5.401);
   model.set_TSigmax(-27.296);
   model.set_TSigmaL(-388.052);
   model.set_BMuPhi(-1.064960e5);
   model.set_vs(2.83138703e4);
   model.set_vsb(2.82546412e4);
   model.set_vphi(-2.00812928e4);
   model.set_Lambdax(0.06993);
   model.set_TLambdax(-8.2418);
   model.set_MuPr(1.440102e4);
   model.set_BMuPr(-1.4807298e7);

   model.set_gD(0,0,0.);
   model.set_gD(0,1,0.);
   model.set_gD(0,2,0.);
   model.set_gD(1,0,0.);
   model.set_gD(1,1,0.);
   model.set_gD(1,2,0.);
   model.set_gD(2,0,0.);
   model.set_gD(2,1,0.);
   model.set_gD(2,2,0.);

   model.set_fu(0,0,1.41478056e-7);
   model.set_fu(0,1,0.);
   model.set_fu(1,0,0.);
   model.set_fu(1,1,1.41478056e-7);
   model.set_fu(2,0,1.41478056e-7);
   model.set_fu(2,1,0.);

   model.set_fd(0,0,1.67921241e-7);
   model.set_fd(0,1,0.);
   model.set_fd(1,0,0.);
   model.set_fd(1,1,1.67921241e-7);
   model.set_fd(2,0,0.);
   model.set_fd(2,1,1.67921241e-7);

   model.set_ThE(0,0,0.);
   model.set_ThE(0,1,0.);
   model.set_ThE(1,0,0.);
   model.set_ThE(1,1,0.);
   model.set_ThE(2,0,0.);
   model.set_ThE(2,1,0.);

   model.set_TgD(0,0,0.);
   model.set_TgD(0,1,0.);
   model.set_TgD(0,2,0.);
   model.set_TgD(1,0,0.);
   model.set_TgD(1,1,0.);
   model.set_TgD(1,2,0.);
   model.set_TgD(2,0,0.);
   model.set_TgD(2,1,0.);
   model.set_TgD(2,2,0.);

   model.set_Tfu(0,0,-2.22261220e-5);
   model.set_Tfu(0,1,0.);
   model.set_Tfu(1,0,0.);
   model.set_Tfu(1,1,-2.22261220e-5);
   model.set_Tfu(2,0,-2.22261220e-5);
   model.set_Tfu(2,1,0.);

   model.set_Tfd(0,0,-2.1934159e-4);
   model.set_Tfd(0,1,0.);
   model.set_Tfd(1,0,0.);
   model.set_Tfd(1,1,-2.1934159e-4);
   model.set_Tfd(2,0,0.);
   model.set_Tfd(2,1,-2.1934159e-4);

   model.set_mq2(0,0,2.49039302e7);
   model.set_mq2(0,1,0.);
   model.set_mq2(0,2,0.);
   model.set_mq2(1,0,0.);
   model.set_mq2(1,1,2.49037451e7);
   model.set_mq2(1,2,0.);
   model.set_mq2(2,0,0.);
   model.set_mq2(2,1,0.);
   model.set_mq2(2,2,1.86697107e7);

   model.set_me2(0,0,1.73783932e7);
   model.set_me2(0,1,0.);
   model.set_me2(0,2,0.);
   model.set_me2(1,0,0.);
   model.set_me2(1,1,1.73774736e7);
   model.set_me2(1,2,0.);
   model.set_me2(2,0,0.);
   model.set_me2(2,1,0.);
   model.set_me2(2,2,1.71172453e7);

   model.set_ml2(0,0,1.81909620e7);
   model.set_ml2(0,1,0.);
   model.set_ml2(0,2,0.);
   model.set_ml2(1,0,0.);
   model.set_ml2(1,1,1.81905060e7);
   model.set_ml2(1,2,0.);
   model.set_ml2(2,0,0.);
   model.set_ml2(2,1,0.);
   model.set_ml2(2,2,1.80614891e7);

   model.set_mu2(0,0,2.41984097e7);
   model.set_mu2(0,1,0.);
   model.set_mu2(0,2,0.);
   model.set_mu2(1,0,0.);
   model.set_mu2(1,1,2.41982117e7);
   model.set_mu2(1,2,0.);
   model.set_mu2(2,0,0.);
   model.set_mu2(2,1,0.);
   model.set_mu2(2,2,1.18579271e7);

   model.set_md2(0,0,2.41711967e7);
   model.set_md2(0,1,0.);
   model.set_md2(0,2,0.);
   model.set_md2(1,0,0.);
   model.set_md2(1,1,2.41710210e7);
   model.set_md2(1,2,0.);
   model.set_md2(2,0,0.);
   model.set_md2(2,1,0.);
   model.set_md2(2,2,2.38399555e7);

   model.set_mHd2(1.75316767e7);
   model.set_mHu2(-1.34277344e6);
   model.set_ms2(1.65701318e7);
   model.set_msbar2(1.70525287e7);
   model.set_mphi2(1.19184877e7);
   model.set_mHp2(1.58327783e7);
   model.set_mHpbar2(1.58347882e7);
   model.set_MassB(483.724424);
   model.set_MassWB(841.541566);
   model.set_MassG(1800.95665);
   model.set_MassBp(443.432203);

   model.set_mH1I2(0,0,1.81521e7);
   model.set_mH1I2(0,1,0.);
   model.set_mH1I2(1,0,0.);
   model.set_mH1I2(1,1,1.81521e7);

   model.set_mH2I2(0,0,1.80988e7);
   model.set_mH2I2(0,1,0.);
   model.set_mH2I2(1,0,0.);
   model.set_mH2I2(1,1,1.80988e7);

   model.set_mSI2(0,0,1.72979e7);
   model.set_mSI2(0,1,0.);
   model.set_mSI2(0,2,0.);
   model.set_mSI2(1,0,0.);
   model.set_mSI2(1,1,1.72979e7);
   model.set_mSI2(1,2,0.);
   model.set_mSI2(2,0,0.);
   model.set_mSI2(2,1,0.);
   model.set_mSI2(2,2,1.72979e7);

   model.set_mDx2(0,0,2.41453332e7);
   model.set_mDx2(0,1,0.);
   model.set_mDx2(0,2,0.);
   model.set_mDx2(1,0,0.);
   model.set_mDx2(1,1,2.41453332e7);
   model.set_mDx2(1,2,0.);
   model.set_mDx2(2,0,0.);
   model.set_mDx2(2,1,0.);
   model.set_mDx2(2,2,2.41453332e7);

   model.set_mDxbar2(0,0,2.41773923e7);
   model.set_mDxbar2(0,1,0.);
   model.set_mDxbar2(0,2,0.);
   model.set_mDxbar2(1,0,0.);
   model.set_mDxbar2(1,1,2.41773923e7);
   model.set_mDxbar2(1,2,0.);
   model.set_mDxbar2(2,0,0.);
   model.set_mDxbar2(2,1,0.);
   model.set_mDxbar2(2,2,2.41773923e7);

   model.set_Kappa(0,0,0.0037300647);
   model.set_Kappa(0,1,0.);
   model.set_Kappa(0,2,0.);
   model.set_Kappa(1,0,0.);
   model.set_Kappa(1,1,0.0037300647);
   model.set_Kappa(1,2,0.);
   model.set_Kappa(2,0,0.);
   model.set_Kappa(2,1,0.);
   model.set_Kappa(2,2,0.0037300647);

   model.set_TKappa(0,0,-19.1766933);
   model.set_TKappa(0,1,0.);
   model.set_TKappa(0,2,0.);
   model.set_TKappa(1,0,0.);
   model.set_TKappa(1,1,-19.1766933);
   model.set_TKappa(1,2,0.);
   model.set_TKappa(2,0,0.);
   model.set_TKappa(2,1,0.);
   model.set_TKappa(2,2,-19.1766933);

   model.set_Lambda12(0,0,0.0837147154);
   model.set_Lambda12(0,1,0.);
   model.set_Lambda12(1,0,0.);
   model.set_Lambda12(1,1,0.0837147154);

   model.set_TLambda12(0,0,-111.017188);
   model.set_TLambda12(0,1,0.);
   model.set_TLambda12(1,0,0.);
   model.set_TLambda12(1,1,-111.017188);
}

BOOST_AUTO_TEST_CASE( test_tadpole_vd )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   model.calculate_DRbar_masses();

   CNE6SSM_higgs_upper_bound upper_bound(model);

   upper_bound.set_include_all_SM_generations(true);

   const double upper_bound_tadpole = upper_bound.get_tadpole_vd();
   const double model_tadpole = Re(partial_tadpole_hh(model, 0));

   BOOST_CHECK_CLOSE(upper_bound_tadpole, model_tadpole, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_tadpole_vu )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   model.calculate_DRbar_masses();

   CNE6SSM_higgs_upper_bound upper_bound(model);

   upper_bound.set_include_all_SM_generations(true);

   const double upper_bound_tadpole = upper_bound.get_tadpole_vu();
   const double model_tadpole = Re(partial_tadpole_hh(model, 1));

   BOOST_CHECK_CLOSE(upper_bound_tadpole, model_tadpole, 1.0e-10);
}

double dV1lp_up_dvd_at_vd(double vd, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vd(vd);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_all_SM_generations(true);
   upper_bound.set_include_up_tadpoles(true);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vd();
}

double dV1lp_up_dvd_at_vu(double vu, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vu(vu);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_all_SM_generations(true);
   upper_bound.set_include_up_tadpoles(true);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vd();
}

double dV1lp_up_dvu_at_vd(double vd, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vd(vd);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_all_SM_generations(true);
   upper_bound.set_include_up_tadpoles(true);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vu();
}

double dV1lp_up_dvu_at_vu(double vu, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vu(vu);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_all_SM_generations(true);
   upper_bound.set_include_up_tadpoles(true);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vu();
}

BOOST_AUTO_TEST_CASE( test_up_contributions )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   upper_bound.set_include_all_SM_generations(true);
   upper_bound.set_include_up_tadpoles(true);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   double Delta00Pr_exact = 0.;
   double Delta01Pr_exact = 0.;
   double Delta11Pr_exact = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      Delta00Pr_exact += upper_bound.get_unrotated_up_contribution(gen, 0, 0);
      Delta01Pr_exact += upper_bound.get_unrotated_up_contribution(gen, 0, 1);
      Delta11Pr_exact += upper_bound.get_unrotated_up_contribution(gen, 1, 1);
   }

   const double h = 1.0e-5;

   double Delta00Pr_approx;
   double Delta00Pr_approx_err;

   gsl_function F1 = {&dV1lp_up_dvd_at_vd, &model};
   gsl_deriv_central(&F1, model.get_vd(), h, &Delta00Pr_approx,
                     &Delta00Pr_approx_err);

   double Delta01Pr_approx;
   double Delta01Pr_approx_err;

   gsl_function F2 = {&dV1lp_up_dvd_at_vu, &model};
   gsl_deriv_central(&F2, model.get_vu(), h, &Delta01Pr_approx,
                     &Delta01Pr_approx_err);

   double Delta11Pr_approx;
   double Delta11Pr_approx_err;

   gsl_function F3 = {&dV1lp_up_dvu_at_vu, &model};
   gsl_deriv_central(&F3, model.get_vu(), h, &Delta11Pr_approx,
                     &Delta11Pr_approx_err);

   BOOST_CHECK_LT(Abs(Delta00Pr_exact - Delta00Pr_approx), Delta00Pr_approx_err);
   BOOST_CHECK_LT(Abs(Delta01Pr_exact - Delta01Pr_approx), Delta01Pr_approx_err);
   BOOST_CHECK_LT(Abs(Delta11Pr_exact - Delta11Pr_approx), Delta11Pr_approx_err);
}

double dV1lp_down_dvd_at_vd(double vd, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vd(vd);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_all_SM_generations(true);
   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(true);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vd();
}

double dV1lp_down_dvd_at_vu(double vu, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vu(vu);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_all_SM_generations(true);
   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(true);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vd();
}

double dV1lp_down_dvu_at_vd(double vd, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vd(vd);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_all_SM_generations(true);
   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(true);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vu();
}

double dV1lp_down_dvu_at_vu(double vu, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vu(vu);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_all_SM_generations(true);
   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(true);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vu();
}

BOOST_AUTO_TEST_CASE( test_down_contributions )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   upper_bound.set_include_all_SM_generations(true);
   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(true);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   double Delta00Pr_exact = 0.;
   double Delta01Pr_exact = 0.;
   double Delta11Pr_exact = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      Delta00Pr_exact += upper_bound.get_unrotated_down_contribution(gen, 0, 0);
      Delta01Pr_exact += upper_bound.get_unrotated_down_contribution(gen, 0, 1);
      Delta11Pr_exact += upper_bound.get_unrotated_down_contribution(gen, 1, 1);
   }

   const double h = 1.0e-5;

   double Delta00Pr_approx;
   double Delta00Pr_approx_err;

   gsl_function F1 = {&dV1lp_down_dvd_at_vd, &model};
   gsl_deriv_central(&F1, model.get_vd(), h, &Delta00Pr_approx,
                     &Delta00Pr_approx_err);

   double Delta01Pr_approx;
   double Delta01Pr_approx_err;

   gsl_function F2 = {&dV1lp_down_dvd_at_vu, &model};
   gsl_deriv_central(&F2, model.get_vu(), h, &Delta01Pr_approx,
                     &Delta01Pr_approx_err);

   double Delta11Pr_approx;
   double Delta11Pr_approx_err;

   gsl_function F3 = {&dV1lp_down_dvu_at_vu, &model};
   gsl_deriv_central(&F3, model.get_vu(), h, &Delta11Pr_approx,
                     &Delta11Pr_approx_err);

   BOOST_CHECK_LT(Abs(Delta00Pr_exact - Delta00Pr_approx), Delta00Pr_approx_err);
   BOOST_CHECK_LT(Abs(Delta01Pr_exact - Delta01Pr_approx), Delta01Pr_approx_err);
   BOOST_CHECK_LT(Abs(Delta11Pr_exact - Delta11Pr_approx), Delta11Pr_approx_err);
}

double dV1lp_exotic_dvd_at_vd(double vd, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vd(vd);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(true);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vd();
}

double dV1lp_exotic_dvd_at_vu(double vu, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vu(vu);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(true);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vd();
}

double dV1lp_exotic_dvu_at_vd(double vd, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vd(vd);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(true);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vu();
}

double dV1lp_exotic_dvu_at_vu(double vu, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vu(vu);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(true);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vu();
}

BOOST_AUTO_TEST_CASE( test_exotic_contributions )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(true);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   double Delta00Pr_exact = 0.;
   double Delta01Pr_exact = 0.;
   double Delta11Pr_exact = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      Delta00Pr_exact +=
         upper_bound.get_unrotated_exotic_contribution(gen, 0, 0);
      Delta01Pr_exact +=
         upper_bound.get_unrotated_exotic_contribution(gen, 0, 1);
      Delta11Pr_exact +=
         upper_bound.get_unrotated_exotic_contribution(gen, 1, 1);
   }

   const double h = 1.0e-5;

   double Delta00Pr_approx;
   double Delta00Pr_approx_err;

   gsl_function F1 = {&dV1lp_exotic_dvd_at_vd, &model};
   gsl_deriv_central(&F1, model.get_vd(), h, &Delta00Pr_approx,
                     &Delta00Pr_approx_err);

   double Delta01Pr_approx;
   double Delta01Pr_approx_err;

   gsl_function F2 = {&dV1lp_exotic_dvd_at_vu, &model};
   gsl_deriv_central(&F2, model.get_vu(), h, &Delta01Pr_approx,
                     &Delta01Pr_approx_err);

   double Delta11Pr_approx;
   double Delta11Pr_approx_err;

   gsl_function F3 = {&dV1lp_exotic_dvu_at_vu, &model};
   gsl_deriv_central(&F3, model.get_vu(), h, &Delta11Pr_approx,
                     &Delta11Pr_approx_err);

   BOOST_CHECK_LT(Abs(Delta00Pr_exact - Delta00Pr_approx), Delta00Pr_approx_err);
   BOOST_CHECK_LT(Abs(Delta01Pr_exact - Delta01Pr_approx), Delta01Pr_approx_err);
   BOOST_CHECK_LT(Abs(Delta11Pr_exact - Delta11Pr_approx), Delta11Pr_approx_err);
}

double dV1lp_inert_singlet_dvd_at_vd(double vd, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vd(vd);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(true);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vd();
}

double dV1lp_inert_singlet_dvd_at_vu(double vu, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vu(vu);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(true);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vd();
}

double dV1lp_inert_singlet_dvu_at_vd(double vd, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vd(vd);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(true);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vu();
}

double dV1lp_inert_singlet_dvu_at_vu(double vu, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vu(vu);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(true);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vu();
}

BOOST_AUTO_TEST_CASE( test_inert_singlet_contributions )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(true);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   double Delta00Pr_exact = 0.;
   double Delta01Pr_exact = 0.;
   double Delta11Pr_exact = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      Delta00Pr_exact +=
         upper_bound.get_unrotated_inert_singlet_contribution(gen, 0, 0);
      Delta01Pr_exact +=
         upper_bound.get_unrotated_inert_singlet_contribution(gen, 0, 1);
      Delta11Pr_exact +=
         upper_bound.get_unrotated_inert_singlet_contribution(gen, 1, 1);
   }

   const double h = 1.0e-5;

   double Delta00Pr_approx;
   double Delta00Pr_approx_err;

   gsl_function F1 = {&dV1lp_inert_singlet_dvd_at_vd, &model};
   gsl_deriv_central(&F1, model.get_vd(), h, &Delta00Pr_approx,
                     &Delta00Pr_approx_err);

   double Delta01Pr_approx;
   double Delta01Pr_approx_err;

   gsl_function F2 = {&dV1lp_inert_singlet_dvd_at_vu, &model};
   gsl_deriv_central(&F2, model.get_vu(), h, &Delta01Pr_approx,
                     &Delta01Pr_approx_err);

   double Delta11Pr_approx;
   double Delta11Pr_approx_err;

   gsl_function F3 = {&dV1lp_inert_singlet_dvu_at_vu, &model};
   gsl_deriv_central(&F3, model.get_vu(), h, &Delta11Pr_approx,
                     &Delta11Pr_approx_err);

   BOOST_CHECK_LT(Abs(Delta00Pr_exact - Delta00Pr_approx), Delta00Pr_approx_err);
   BOOST_CHECK_LT(Abs(Delta01Pr_exact - Delta01Pr_approx), Delta01Pr_approx_err);
   BOOST_CHECK_LT(Abs(Delta11Pr_exact - Delta11Pr_approx), Delta11Pr_approx_err);
}

double dV1lp_inert_neutral_higgs_dvd_at_vd(double vd, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vd(vd);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(true);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vd();
}

double dV1lp_inert_neutral_higgs_dvd_at_vu(double vu, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vu(vu);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(true);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vd();
}

double dV1lp_inert_neutral_higgs_dvu_at_vd(double vd, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vd(vd);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(true);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vu();
}

double dV1lp_inert_neutral_higgs_dvu_at_vu(double vu, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vu(vu);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(true);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   return -upper_bound.get_tadpole_vu();
}

BOOST_AUTO_TEST_CASE( test_inert_neutral_higgs_contributions )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(true);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   double Delta00Pr_exact = 0.;
   double Delta01Pr_exact = 0.;
   double Delta11Pr_exact = 0.;

   for (unsigned gen = 0; gen < 2; ++gen) {
      Delta00Pr_exact +=
         upper_bound.get_unrotated_inert_neutral_higgs_contribution(gen, 0, 0);
      Delta01Pr_exact +=
         upper_bound.get_unrotated_inert_neutral_higgs_contribution(gen, 0, 1);
      Delta11Pr_exact +=
         upper_bound.get_unrotated_inert_neutral_higgs_contribution(gen, 1, 1);
   }

   const double h = 1.0e-5;

   double Delta00Pr_approx;
   double Delta00Pr_approx_err;

   gsl_function F1 = {&dV1lp_inert_neutral_higgs_dvd_at_vd, &model};
   gsl_deriv_central(&F1, model.get_vd(), h, &Delta00Pr_approx,
                     &Delta00Pr_approx_err);

   double Delta01Pr_approx;
   double Delta01Pr_approx_err;

   gsl_function F2 = {&dV1lp_inert_neutral_higgs_dvd_at_vu, &model};
   gsl_deriv_central(&F2, model.get_vu(), h, &Delta01Pr_approx,
                     &Delta01Pr_approx_err);

   double Delta11Pr_approx;
   double Delta11Pr_approx_err;

   gsl_function F3 = {&dV1lp_inert_neutral_higgs_dvu_at_vu, &model};
   gsl_deriv_central(&F3, model.get_vu(), h, &Delta11Pr_approx,
                     &Delta11Pr_approx_err);

   BOOST_CHECK_LT(Abs(Delta00Pr_exact - Delta00Pr_approx), Delta00Pr_approx_err);
   BOOST_CHECK_LT(Abs(Delta01Pr_exact - Delta01Pr_approx), Delta01Pr_approx_err);
   BOOST_CHECK_LT(Abs(Delta11Pr_exact - Delta11Pr_approx), Delta11Pr_approx_err);
}

double dV1lp_inert_charged_higgs_dvd_at_vd(double vd, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vd(vd);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(true);

   return -upper_bound.get_tadpole_vd();
}

double dV1lp_inert_charged_higgs_dvd_at_vu(double vu, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vu(vu);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(true);

   return -upper_bound.get_tadpole_vd();
}

double dV1lp_inert_charged_higgs_dvu_at_vd(double vd, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vd(vd);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(true);

   return -upper_bound.get_tadpole_vu();
}

double dV1lp_inert_charged_higgs_dvu_at_vu(double vu, void* params)
{
   CNE6SSM_mass_eigenstates* model
      = static_cast<CNE6SSM_mass_eigenstates*>(params);

   model->set_vu(vu);

   CNE6SSM_higgs_upper_bound upper_bound(*model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(true);

   return -upper_bound.get_tadpole_vu();
}

BOOST_AUTO_TEST_CASE( test_inert_charged_higgs_contributions )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(true);

   double Delta00Pr_exact = 0.;
   double Delta01Pr_exact = 0.;
   double Delta11Pr_exact = 0.;

   for (unsigned gen = 0; gen < 2; ++gen) {
      Delta00Pr_exact +=
         upper_bound.get_unrotated_inert_charged_higgs_contribution(gen, 0, 0);
      Delta01Pr_exact +=
         upper_bound.get_unrotated_inert_charged_higgs_contribution(gen, 0, 1);
      Delta11Pr_exact +=
         upper_bound.get_unrotated_inert_charged_higgs_contribution(gen, 1, 1);
   }

   const double h = 1.0e-5;

   double Delta00Pr_approx;
   double Delta00Pr_approx_err;

   gsl_function F1 = {&dV1lp_inert_charged_higgs_dvd_at_vd, &model};
   gsl_deriv_central(&F1, model.get_vd(), h, &Delta00Pr_approx,
                     &Delta00Pr_approx_err);

   double Delta01Pr_approx;
   double Delta01Pr_approx_err;

   gsl_function F2 = {&dV1lp_inert_charged_higgs_dvd_at_vu, &model};
   gsl_deriv_central(&F2, model.get_vu(), h, &Delta01Pr_approx,
                     &Delta01Pr_approx_err);

   double Delta11Pr_approx;
   double Delta11Pr_approx_err;

   gsl_function F3 = {&dV1lp_inert_charged_higgs_dvu_at_vu, &model};
   gsl_deriv_central(&F3, model.get_vu(), h, &Delta11Pr_approx,
                     &Delta11Pr_approx_err);

   BOOST_CHECK_LT(Abs(Delta00Pr_exact - Delta00Pr_approx), Delta00Pr_approx_err);
   BOOST_CHECK_LT(Abs(Delta01Pr_exact - Delta01Pr_approx), Delta01Pr_approx_err);
   BOOST_CHECK_LT(Abs(Delta11Pr_exact - Delta11Pr_approx), Delta11Pr_approx_err);
}

double get_tadpole_hh_0(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double gbar = Sqrt(Sqr(g2) + 0.6 * Sqr(g1));
   const double g1p = model.get_g1p();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();

   CNE6SSM_higgs_upper_bound upper_bound(model);

   A0_fn A0(model.get_scale());

   const double QHd = -3.0;
   const double QHu = -2.0;
   const double QQ = 1.0;
   const double Qu = 1.0;
   const double Qd = 2.0;
   const double QDX = -2.0;
   const double QDXbar = -3.0;
   const double QS = model.get_QS();

   double t1_up = 0.;
   double t1_down = 0.;
   double t1_exotic = 0.;
   double t1_inert_singlet = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      const double yu = model.get_Yu(gen, gen);
      const Eigen::Array<double,2,1> MSu2(upper_bound.calculate_MSu2(gen));
      const double Sin2ThetaSu = upper_bound.calculate_Sin2ThetaSu(gen);
      const double Cos2ThetaSu = upper_bound.calculate_Cos2ThetaSu(gen);

      t1_up += 1.5 * oneOver16PiSqr * vd * ((A0(Sqrt(MSu2(0))) + A0(Sqrt(MSu2(1))))
         * (0.25 * Sqr(gbar) + 0.025 * Sqr(g1p) * QHd * (QQ + Qu)) +
         (A0(Sqrt(MSu2(0))) - A0(Sqrt(MSu2(1)))) * (Lambdax * yu * vs * Sin2ThetaSu
         / vd + Cos2ThetaSu * (0.25 * (Sqr(g2) - Sqr(g1)) + 0.025 * Sqr(g1p) * QHd
         * (QQ - Qu))));

      const double yd = model.get_Yd(gen, gen);
      const double Tyd = model.get_TYd(gen, gen);
      const Eigen::Array<double,2,1> MSd2(upper_bound.calculate_MSd2(gen));
      const double MFd2 = upper_bound.calculate_MFd2(gen);
      const double Sin2ThetaSd = upper_bound.calculate_Sin2ThetaSd(gen);
      const double Cos2ThetaSd = upper_bound.calculate_Cos2ThetaSd(gen);

      t1_down += -6.0 * oneOver16PiSqr * vd * Sqr(yd) * A0(Sqrt(MFd2)) + 1.5 *
         oneOver16PiSqr * vd * ((A0(Sqrt(MSd2(0))) + A0(Sqrt(MSd2(1)))) * (-0.25
         * Sqr(gbar) + 0.025 * Sqr(g1p) * QHd * (QQ + Qd) + 2.0 * Sqr(yd)) +
         (A0(Sqrt(MSd2(1))) - A0(Sqrt(MSd2(0)))) * (Sqrt(2.0) * Tyd * Sin2ThetaSd
         / vd + Cos2ThetaSd * (0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) - 0.025 * Sqr(g1p)
         * QHd * (QQ - Qd))));

      const double Kappa = model.get_Kappa(gen, gen);
      const Eigen::Array<double,2,1> MSDX2(upper_bound.calculate_MSDX2(gen));
      const double Sin2ThetaSDX = upper_bound.calculate_Sin2ThetaSDX(gen);
      const double Cos2ThetaSDX = upper_bound.calculate_Cos2ThetaSDX(gen);

      t1_exotic += 1.5 * oneOver16PiSqr * vd * (0.025 * Sqr(g1p) * QHd * (QDX +
         QDXbar) * (A0(Sqrt(MSDX2(0))) + A0(Sqrt(MSDX2(1)))) + (A0(Sqrt(MSDX2(0)))
         - A0(Sqrt(MSDX2(1)))) * (Kappa * Lambdax * vu * Sin2ThetaSDX / vd +
         Cos2ThetaSDX * (0.2 * Sqr(g1) + 0.025 * Sqr(g1p) * QHd * (QDX - QDXbar))));

      const double MSI02 = upper_bound.calculate_MSI02(gen);

      t1_inert_singlet += 0.025 * oneOver16PiSqr * vd * Sqr(g1p) * QHd * QS *
         A0(Sqrt(MSI02));
   }

   double t1_inert_neutral = 0.;
   double t1_inert_charged = 0.;

   for (unsigned gen = 0; gen < 2; ++gen) {
      const double Lambda = model.get_Lambda12(gen, gen);
      const Eigen::Array<double,2,1> MHI02(upper_bound.calculate_MHI02(gen));
      const double Sin2ThetaHI0 = upper_bound.calculate_Sin2ThetaHI0(gen);
      const double Cos2ThetaHI0 = upper_bound.calculate_Cos2ThetaHI0(gen);

      t1_inert_neutral += 0.5 * oneOver16PiSqr * vd * (0.025 * Sqr(g1p) * QHd *
         (QHd + QHu) * (A0(Sqrt(MHI02(0))) + A0(Sqrt(MHI02(1)))) + (A0(Sqrt(MHI02(1)))
         - A0(Sqrt(MHI02(0)))) * (Lambda * Lambdax * vu * Sin2ThetaHI0 / vd -
         Cos2ThetaHI0 * (0.5 * Sqr(gbar) + 0.025 * Sqr(g1p) * QHd * (QHd - QHu))));

      const Eigen::Array<double,2,1> MHIPM2(upper_bound.calculate_MHIPM2(gen));
      const double Sin2ThetaHIPM = upper_bound.calculate_Sin2ThetaHIPM(gen);
      const double Cos2ThetaHIPM = upper_bound.calculate_Cos2ThetaHIPM(gen);

      t1_inert_charged += 0.5 * oneOver16PiSqr * vd * (0.025 * Sqr(g1p) * QHd *
         (QHd + QHu) * (A0(Sqrt(MHIPM2(0))) + A0(Sqrt(MHIPM2(1)))) +
         (A0(Sqrt(MHIPM2(0))) - A0(Sqrt(MHIPM2(1)))) * (Lambda * Lambdax * vu *
         Sin2ThetaHIPM / vd - Cos2ThetaHIPM * (0.5 * (Sqr(g2) - 0.6 * Sqr(g1)) -
         0.025 * Sqr(g1p) * QHd * (QHd - QHu))));
   }

   return t1_up + t1_down + t1_exotic + t1_inert_neutral + t1_inert_singlet
      + t1_inert_charged;
}

double get_tadpole_hh_1(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double gbar = Sqrt(Sqr(g2) + 0.6 * Sqr(g1));
   const double g1p = model.get_g1p();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();

   CNE6SSM_higgs_upper_bound upper_bound(model);

   A0_fn A0(model.get_scale());

   const double QHd = -3.0;
   const double QHu = -2.0;
   const double QQ = 1.0;
   const double Qu = 1.0;
   const double Qd = 2.0;
   const double QDX = -2.0;
   const double QDXbar = -3.0;
   const double QS = model.get_QS();

   double t2_up = 0.;
   double t2_down = 0.;
   double t2_exotic = 0.;
   double t2_inert_singlet = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      const double yu = model.get_Yu(gen, gen);
      const double Tyu = model.get_TYu(gen, gen);
      const Eigen::Array<double,2,1> MSu2(upper_bound.calculate_MSu2(gen));
      const double MFu2 = upper_bound.calculate_MFu2(gen);
      const double Sin2ThetaSu = upper_bound.calculate_Sin2ThetaSu(gen);
      const double Cos2ThetaSu = upper_bound.calculate_Cos2ThetaSu(gen);

      t2_up += -6.0 * oneOver16PiSqr * vu * Sqr(yu) * A0(Sqrt(MFu2)) + 1.5 *
         oneOver16PiSqr * vu * ((A0(Sqrt(MSu2(0))) + A0(Sqrt(MSu2(1)))) * (-0.25 *
         Sqr(gbar) + 0.025 * Sqr(g1p) * QHu * (QQ + Qu) + 2.0 * Sqr(yu)) +
         (A0(Sqrt(MSu2(1))) - A0(Sqrt(MSu2(0)))) * (Sqrt(2.0) * Tyu * Sin2ThetaSu
         / vu + Cos2ThetaSu * (0.25 * (Sqr(g2) - Sqr(g1)) - 0.025 * Sqr(g1p) * QHu
         * (QQ - Qu))));

      const double yd = model.get_Yd(gen, gen);
      const Eigen::Array<double,2,1> MSd2(upper_bound.calculate_MSd2(gen));
      const double Sin2ThetaSd = upper_bound.calculate_Sin2ThetaSd(gen);
      const double Cos2ThetaSd = upper_bound.calculate_Cos2ThetaSd(gen);

      t2_down += 1.5 * oneOver16PiSqr * vu * ((A0(Sqrt(MSd2(0))) + A0(Sqrt(MSd2(1))))
         * (0.25 * Sqr(gbar) + 0.025 * Sqr(g1p) * QHu * (QQ + Qd)) + (A0(Sqrt(MSd2(0)))
         - A0(Sqrt(MSd2(1)))) * (Lambdax * yd * vs * Sin2ThetaSd / vu + Cos2ThetaSd *
         (0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) + 0.025 * Sqr(g1p) * QHu * (QQ - Qd))));

      const double Kappa = model.get_Kappa(gen, gen);
      const Eigen::Array<double,2,1> MSDX2(upper_bound.calculate_MSDX2(gen));
      const double Sin2ThetaSDX = upper_bound.calculate_Sin2ThetaSDX(gen);
      const double Cos2ThetaSDX = upper_bound.calculate_Cos2ThetaSDX(gen);

      t2_exotic += 1.5 * oneOver16PiSqr * vu * (0.025 * Sqr(g1p) * QHu * (QDX + QDXbar)
         * (A0(Sqrt(MSDX2(0))) + A0(Sqrt(MSDX2(1)))) + (A0(Sqrt(MSDX2(0))) -
         A0(Sqrt(MSDX2(1)))) * (Kappa * Lambdax * vd * Sin2ThetaSDX / vu + Cos2ThetaSDX
         * (-0.2 * Sqr(g1) + 0.025 * Sqr(g1p) * QHu * (QDX - QDXbar))));

      const double MSI02 = upper_bound.calculate_MSI02(gen);

      t2_inert_singlet += 0.025 * oneOver16PiSqr * vu * Sqr(g1p) * QHu * QS *
         A0(Sqrt(MSI02));
   }

   double t2_inert_neutral = 0.;
   double t2_inert_charged = 0.;

   for (unsigned gen = 0; gen < 2; ++gen) {
      const double Lambda = model.get_Lambda12(gen, gen);
      const Eigen::Array<double,2,1> MHI02(upper_bound.calculate_MHI02(gen));
      const double Sin2ThetaHI0 = upper_bound.calculate_Sin2ThetaHI0(gen);
      const double Cos2ThetaHI0 = upper_bound.calculate_Cos2ThetaHI0(gen);

      t2_inert_neutral += 0.5 * oneOver16PiSqr * vu * (0.025 * Sqr(g1p) * QHu *
         (QHd + QHu) * (A0(Sqrt(MHI02(0))) + A0(Sqrt(MHI02(1)))) + (A0(Sqrt(MHI02(1)))
         - A0(Sqrt(MHI02(0)))) * (Lambda * Lambdax * vd * Sin2ThetaHI0 / vu +
         Cos2ThetaHI0 * (0.5 * Sqr(gbar) - 0.025 * Sqr(g1p) * QHu * (QHd - QHu))));

      const Eigen::Array<double,2,1> MHIPM2(upper_bound.calculate_MHIPM2(gen));
      const double Sin2ThetaHIPM = upper_bound.calculate_Sin2ThetaHIPM(gen);
      const double Cos2ThetaHIPM = upper_bound.calculate_Cos2ThetaHIPM(gen);

      t2_inert_charged += 0.5 * oneOver16PiSqr * vu * (0.025 * Sqr(g1p) * QHu *
         (QHd + QHu) * (A0(Sqrt(MHIPM2(0))) + A0(Sqrt(MHIPM2(1)))) +
         (A0(Sqrt(MHIPM2(0))) - A0(Sqrt(MHIPM2(1)))) * (Lambda * Lambdax * vd *
         Sin2ThetaHIPM / vu + Cos2ThetaHIPM * (0.5 * (Sqr(g2) - 0.6 * Sqr(g1)) +
         0.025 * Sqr(g1p) * QHu * (QHd - QHu))));
   }

   return t2_up + t2_down + t2_exotic + t2_inert_neutral + t2_inert_singlet
      + t2_inert_charged;
}

BOOST_AUTO_TEST_CASE( test_analytic_tadpoles )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   upper_bound.set_include_all_SM_generations(true);
   upper_bound.set_include_up_tadpoles(true);
   upper_bound.set_include_down_tadpoles(true);
   upper_bound.set_include_exotic_tadpoles(true);
   upper_bound.set_include_inert_singlet_tadpoles(true);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(true);
   upper_bound.set_include_inert_charged_higgs_tadpoles(true);

   const double t1 = upper_bound.get_tadpole_vd();
   const double t2 = upper_bound.get_tadpole_vu();

   const double expect_t1 = get_tadpole_hh_0(model);
   const double expect_t2 = get_tadpole_hh_1(model);

   BOOST_CHECK_CLOSE(t1, expect_t1, 1.0e-10);
   BOOST_CHECK_CLOSE(t2, expect_t2, 1.0e-10);
}

double get_d2V1lp_up_dvd_dvd(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vs = model.get_vs();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   const double gbar = Sqrt(Sqr(g2) + 0.6 * Sqr(g1));

   const double QHd = -3.0;
   const double QQ = 1.0;
   const double Qu = 1.0;

   const double scale = model.get_scale();
   A0_fn A0(scale);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double result = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      const Eigen::Array<double,2,1> MSu2(upper_bound.calculate_MSu2(gen));
      const double Sin2ThetaSu = upper_bound.calculate_Sin2ThetaSu(gen);
      const double Cos2ThetaSu = upper_bound.calculate_Cos2ThetaSu(gen);
      const double Sin4ThetaSu = 2.0 * Sin2ThetaSu * Cos2ThetaSu;

      const double inverse_mass_diff = 1.0 / (MSu2(1) - MSu2(0));

      const double yf = model.get_Yu(gen, gen);

      result += 1.5 * oneOver16PiSqr * (0.5 * Sqr(vd) * Log(MSu2(0)
         * MSu2(1) / (Power(scale, 4.0))) * (Sqr(0.25 * Sqr(gbar) +
         0.025 * Sqr(g1p) * QHd * (QQ + Qu)) + Sqr(Lambdax * yf * vs
         * Sin2ThetaSu / vd + Cos2ThetaSu * (0.25 * (Sqr(g2) - Sqr(g1))
         + 0.025 * Sqr(g1p) * QHd * (QQ - Qu)))) - Sqr(vd) * Log(MSu2(1)
         / MSu2(0)) * (0.25 * Sqr(gbar) + 0.025 * Sqr(g1p) * QHd * (QQ +
         Qu)) * (Lambdax * yf * vs * Sin2ThetaSu / vd + Cos2ThetaSu * (0.25
         * (Sqr(g2) - Sqr(g1)) + 0.025 * Sqr(g1p) * QHd * (QQ - Qu))) - (
         0.25 * Sqr(gbar) + 0.025 * Sqr(g1p) * QHd * (QQ + Qu)) * (A0(
         Sqrt(MSu2(0))) + A0(Sqrt(MSu2(1)))) + (A0(Sqrt(MSu2(1))) - A0(Sqrt(
         MSu2(0)))) * (Cos2ThetaSu * (0.25 * (Sqr(g2) - Sqr(g1)) + 0.025 *
         Sqr(g1p) * QHd * (QQ - Qu)) - Sqr(vd) * inverse_mass_diff * (Sqr(
         Sin2ThetaSu) * Sqr(0.25 * (Sqr(g2) - Sqr(g1)) + 0.025 * Sqr(g1p) *
         QHd * (QQ - Qu)) + Sqr(Lambdax) * Sqr(yf) * Sqr(vs) * Sqr(Cos2ThetaSu)
         / Sqr(vd) - Lambdax * yf * vs * Sin4ThetaSu * (0.25 * (Sqr(g2) -
         Sqr(g1)) + 0.025 * Sqr(g1p) * QHd * (QQ - Qu)) / vd)));
   }

   return result;
}

BOOST_AUTO_TEST_CASE( test_analytic_up_contributions )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double Delta00Pr = 0.;
   for (unsigned gen = 0; gen < 3; ++gen) {
      Delta00Pr += upper_bound.get_unrotated_up_contribution(gen, 0, 0);
   }

   const double Delta00Pr_expect = get_d2V1lp_up_dvd_dvd(model);

   BOOST_CHECK_CLOSE(Delta00Pr, Delta00Pr_expect, 1.0e-10);
}

double get_d2V1lp_down_dvu_dvu(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   const double gbar = Sqrt(Sqr(g2) + 0.6 * Sqr(g1));

   const double QHu = -2.0;
   const double QQ = 1.0;
   const double Qd = 2.0;

   const double scale = model.get_scale();
   A0_fn A0(scale);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double result = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      const Eigen::Array<double,2,1> MSd2(upper_bound.calculate_MSd2(gen));
      const double Sin2ThetaSd = upper_bound.calculate_Sin2ThetaSd(gen);
      const double Cos2ThetaSd = upper_bound.calculate_Cos2ThetaSd(gen);
      const double Sin4ThetaSd = 2.0 * Sin2ThetaSd * Cos2ThetaSd;

      const double inverse_mass_diff = 1.0 / (MSd2(1) - MSd2(0));

      const double yf = model.get_Yd(gen, gen);

      result += 1.5 * oneOver16PiSqr * (0.5 * Sqr(vu) * Log(MSd2(0) *
         MSd2(1) / Power(scale, 4)) * (Sqr(0.25 * Sqr(gbar) + 0.025 *
         Sqr(g1p) * QHu * (QQ + Qd)) + Sqr(Lambdax * yf * vs *
         Sin2ThetaSd / vu + Cos2ThetaSd * (0.25 * (Sqr(g2) - 0.2 *
         Sqr(g1)) + 0.025 * Sqr(g1p) * QHu * (QQ - Qd)))) - Sqr(vu) *
         Log(MSd2(1) / MSd2(0)) * (0.25 * Sqr(gbar) + 0.025 * Sqr(g1p)
         * QHu * (QQ + Qd)) * (Lambdax * yf * vs * Sin2ThetaSd / vu +
         Cos2ThetaSd * (0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) + 0.025 *
         Sqr(g1p) * QHu * (QQ - Qd))) - (A0(Sqrt(MSd2(0))) + A0(Sqrt(
         MSd2(1)))) * (0.25 * Sqr(gbar) + 0.025 * Sqr(g1p) * QHu * (QQ
         + Qd)) + (A0(Sqrt(MSd2(1))) - A0(Sqrt(MSd2(0)))) * (Cos2ThetaSd
         * (0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) + 0.025 * Sqr(g1p) * QHu *
         (QQ - Qd)) - Sqr(vu) * inverse_mass_diff * (Sqr(Sin2ThetaSd) *
         Sqr(0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) + 0.025 * Sqr(g1p) * QHu *
         (QQ - Qd)) + Sqr(Lambdax) * Sqr(yf) * Sqr(vs) * Sqr(Cos2ThetaSd)
         / Sqr(vu) - Lambdax * yf * vs * Sin4ThetaSd * (0.25 * (Sqr(g2) -
         0.2 * Sqr(g1)) + 0.025 * Sqr(g1p) * QHu * (QQ - Qd)) / vu)));
   }

   return result;
}

BOOST_AUTO_TEST_CASE( test_analytic_down_contributions )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double Delta11Pr = 0.;
   for (unsigned gen = 0; gen < 3; ++gen) {
      Delta11Pr += upper_bound.get_unrotated_down_contribution(gen, 1, 1);
   }

   const double Delta11Pr_expect = get_d2V1lp_down_dvu_dvu(model);

   BOOST_CHECK_CLOSE(Delta11Pr, Delta11Pr_expect, 1.0e-10);
}

double get_d2V1lp_exotic_dvd_dvd(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();

   const double QHd = -3.0;
   const double QDX = -2.0;
   const double QDXbar = -3.0;

   const double scale = model.get_scale();
   A0_fn A0(scale);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double result = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      const Eigen::Array<double,2,1> MSDX2(upper_bound.calculate_MSDX2(gen));
      const double Sin2ThetaSDX = upper_bound.calculate_Sin2ThetaSDX(gen);
      const double Cos2ThetaSDX = upper_bound.calculate_Cos2ThetaSDX(gen);
      const double Sin4ThetaSDX = 2.0 * Sin2ThetaSDX * Cos2ThetaSDX;

      const double inverse_mass_diff = 1.0 / (MSDX2(1) - MSDX2(0));

      const double Kappa = model.get_Kappa(gen, gen);

      result += 1.5 * oneOver16PiSqr * (0.5 * Sqr(vd) * Log(MSDX2(0) * MSDX2(1)
         / Power(scale, 4)) * (0.000625 * Power(g1p, 4) * Sqr(QHd) * Sqr(QDX +
         QDXbar) + Sqr(Kappa * Lambdax * vu * Sin2ThetaSDX / vd + Cos2ThetaSDX
         * (0.2 * Sqr(g1) + 0.025 * Sqr(g1p) * QHd * (QDX - QDXbar)))) - 0.025
         * Sqr(g1p) * QHd * (QDX + QDXbar) * Sqr(vd) * Log(MSDX2(1) / MSDX2(0))
         * (Kappa * Lambdax * vu * Sin2ThetaSDX / vd + Cos2ThetaSDX * (0.2 * Sqr(
         g1) + 0.025 * Sqr(g1p) * QHd * (QDX - QDXbar))) - 0.025 * Sqr(g1p) * QHd
         * (QDX + QDXbar) * (A0(Sqrt(MSDX2(0))) + A0(Sqrt(MSDX2(1)))) + (A0(Sqrt(
         MSDX2(1))) - A0(Sqrt(MSDX2(0)))) * (Cos2ThetaSDX * (0.2 * Sqr(g1) + 0.025
         * Sqr(g1p) * QHd * (QDX - QDXbar)) - Sqr(vd) * inverse_mass_diff * (Sqr(
         Sin2ThetaSDX) * Sqr(0.2 * Sqr(g1) + 0.025 * Sqr(g1p) * QHd * (QDX - QDXbar)
         ) + Sqr(Kappa) * Sqr(Lambdax) * Sqr(vu) * Sqr(Cos2ThetaSDX) / Sqr(vd) -
         Kappa * Lambdax * vu * Sin4ThetaSDX * (0.2 * Sqr(g1) + 0.025 * Sqr(g1p) *
         QHd * (QDX - QDXbar)) / vd )));
   }

   return result;
}

double get_d2V1lp_exotic_dvd_dvu(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();

   const double QHd = -3.0;
   const double QHu = -2.0;
   const double QDX = -2.0;
   const double QDXbar = -3.0;

   const double scale = model.get_scale();
   A0_fn A0(scale);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double result = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      const Eigen::Array<double,2,1> MSDX2(upper_bound.calculate_MSDX2(gen));
      const double Sin2ThetaSDX = upper_bound.calculate_Sin2ThetaSDX(gen);
      const double Cos2ThetaSDX = upper_bound.calculate_Cos2ThetaSDX(gen);
      const double Sin4ThetaSDX = 2.0 * Sin2ThetaSDX * Cos2ThetaSDX;

      const double inverse_mass_diff = 1.0 / (MSDX2(1) - MSDX2(0));

      const double Kappa = model.get_Kappa(gen, gen);

      result += 1.5 * oneOver16PiSqr * (0.5 * vd * vu * Log(MSDX2(0) * MSDX2(1)
         / Power(scale, 4)) * (0.000625 * Power(g1p, 4) * QHd * QHu * Sqr(QDX +
         QDXbar) + (Kappa * Lambdax * vu * Sin2ThetaSDX / vd + Cos2ThetaSDX *
         (0.2 * Sqr(g1) + 0.025 * Sqr(g1p) * QHd * (QDX - QDXbar))) * (Kappa *
         Lambdax * vd * Sin2ThetaSDX / vu - Cos2ThetaSDX * (0.2 * Sqr(g1) - 0.025
         * Sqr(g1p) * QHu * (QDX - QDXbar)))) - 0.0125 * Sqr(g1p) * (QDX + QDXbar)
         * vd * vu * Log(MSDX2(1) / MSDX2(0)) * (Kappa * Lambdax * (QHd * vd / vu
         + QHu * vu / vd) * Sin2ThetaSDX + 0.2 * Sqr(g1) * (QHu - QHd) *
         Cos2ThetaSDX + 0.05 * Sqr(g1p) * QHu * QHd * Cos2ThetaSDX * (QDX - QDXbar))
         + (A0(Sqrt(MSDX2(1))) - A0(Sqrt(MSDX2(0)))) * (Kappa * Lambdax *
         Sin2ThetaSDX + vd * vu * inverse_mass_diff * (Sqr(Sin2ThetaSDX) * (0.2 *
         Sqr(g1) + 0.025 * Sqr(g1p) * QHd * (QDX - QDXbar)) * (0.2 * Sqr(g1) - 0.025
         * Sqr(g1p) * QHu * (QDX - QDXbar)) - Sqr(Kappa) * Sqr(Lambdax) * Sqr(
         Cos2ThetaSDX) + 0.5 * Kappa * Lambdax * Sin4ThetaSDX * (vd * (0.2 * Sqr(g1)
         + 0.025 * Sqr(g1p) * QHd * (QDX- QDXbar)) / vu - vu * (0.2 * Sqr(g1) - 0.025
         * Sqr(g1p) * QHu * (QDX - QDXbar)) / vd ))));
   }

   return result;
}

double get_d2V1lp_exotic_dvu_dvu(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();

   const double QHu = -2.0;
   const double QDX = -2.0;
   const double QDXbar = -3.0;

   const double scale = model.get_scale();
   A0_fn A0(scale);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double result = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      const Eigen::Array<double,2,1> MSDX2(upper_bound.calculate_MSDX2(gen));
      const double Sin2ThetaSDX = upper_bound.calculate_Sin2ThetaSDX(gen);
      const double Cos2ThetaSDX = upper_bound.calculate_Cos2ThetaSDX(gen);
      const double Sin4ThetaSDX = 2.0 * Sin2ThetaSDX * Cos2ThetaSDX;

      const double inverse_mass_diff = 1.0 / (MSDX2(1) - MSDX2(0));

      const double Kappa = model.get_Kappa(gen, gen);

      result += 1.5 * oneOver16PiSqr * (0.5 * Sqr(vu) * Log(MSDX2(0) * MSDX2(1)
         / Power(scale, 4)) * (0.000625 * Power(g1p, 4) * Sqr(QHu) * Sqr(QDX +
         QDXbar) + Sqr(Kappa * Lambdax * vd * Sin2ThetaSDX / vu - Cos2ThetaSDX
         * (0.2 * Sqr(g1) - 0.025 * Sqr(g1p) * QHu * (QDX - QDXbar)))) - 0.025
         * Sqr(g1p) * QHu * (QDX + QDXbar) * Sqr(vu) * Log(MSDX2(1) / MSDX2(0))
         * (Kappa * Lambdax * vd * Sin2ThetaSDX / vu - Cos2ThetaSDX * (0.2 *
         Sqr(g1) - 0.025 * Sqr(g1p) * QHu * (QDX - QDXbar))) - 0.025 * Sqr(g1p)
         * QHu * (QDX + QDXbar) * (A0(Sqrt(MSDX2(0))) + A0(Sqrt(MSDX2(1)))) - (
         A0(Sqrt(MSDX2(1))) - A0(Sqrt(MSDX2(0)))) * (Cos2ThetaSDX * (0.2 * Sqr(
         g1) - 0.025 * Sqr(g1p) * QHu * (QDX - QDXbar)) + Sqr(vu) *
         inverse_mass_diff * (Sqr(Sin2ThetaSDX) * Sqr(0.2 * Sqr(g1) - 0.025 *
         Sqr(g1p) * QHu * (QDX - QDXbar)) + Sqr(Kappa) * Sqr(Lambdax) * Sqr(vd)
         * Sqr(Cos2ThetaSDX) / Sqr(vu) + Kappa * Lambdax * vd * Sin4ThetaSDX *
         (0.2 * Sqr(g1) - 0.025 * Sqr(g1p) * QHu * (QDX - QDXbar)) / vu )));
   }

   return result;
}

BOOST_AUTO_TEST_CASE( test_analytic_exotic_contributions )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double Delta00Pr = 0.;
   double Delta01Pr = 0.;
   double Delta11Pr = 0.;
   for (unsigned gen = 0; gen < 3; ++gen) {
      Delta00Pr += upper_bound.get_unrotated_exotic_contribution(gen, 0, 0);
      Delta01Pr += upper_bound.get_unrotated_exotic_contribution(gen, 0, 1);
      Delta11Pr += upper_bound.get_unrotated_exotic_contribution(gen, 1, 1);
   }

   const double Delta00Pr_expect = get_d2V1lp_exotic_dvd_dvd(model);
   const double Delta01Pr_expect = get_d2V1lp_exotic_dvd_dvu(model);
   const double Delta11Pr_expect = get_d2V1lp_exotic_dvu_dvu(model);

   BOOST_CHECK_CLOSE(Delta00Pr, Delta00Pr_expect, 1.0e-10);
   BOOST_CHECK_CLOSE(Delta01Pr, Delta01Pr_expect, 1.0e-10);
   BOOST_CHECK_CLOSE(Delta11Pr, Delta11Pr_expect, 1.0e-10);
}

double get_d2V1lp_inert_singlet_dvd_dvd(const CNE6SSM_mass_eigenstates& model)
{
   const double vd = model.get_vd();
   const double g1p = model.get_g1p();

   const double QHd = -3.0;
   const double QS = model.get_QS();

   CNE6SSM_higgs_upper_bound upper_bound(model);

   const double scale = model.get_scale();
   A0_fn A0(scale);

   double result = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      const double MSI02 = upper_bound.calculate_MSI02(gen);

      result += oneOver16PiSqr * 0.025 * Sqr(g1p) * QHd * QS *
         (0.025 * Sqr(g1p) * QHd * QS * Sqr(vd) * Log(MSI02 /
         Sqr(scale)) - A0(Sqrt(MSI02)));
   }

   return result;
}

double get_d2V1lp_inert_singlet_dvd_dvu(const CNE6SSM_mass_eigenstates& model)
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1p = model.get_g1p();

   const double QHd = -3.0;
   const double QHu = -2.0;
   const double QS = model.get_QS();

   const double scale = model.get_scale();

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double result = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      const double MSI02 = upper_bound.calculate_MSI02(gen);

      result += 0.000625 * oneOver16PiSqr * Power(g1p, 4) *
         QHd * QHu * Sqr(QS) * vd * vu * Log(MSI02 /
         Sqr(scale));
   }

   return result;
}

double get_d2V1lp_inert_singlet_dvu_dvu(const CNE6SSM_mass_eigenstates& model)
{
   const double vu = model.get_vu();
   const double g1p = model.get_g1p();

   const double QHu = -2.0;
   const double QS = model.get_QS();

   CNE6SSM_higgs_upper_bound upper_bound(model);

   const double scale = model.get_scale();
   A0_fn A0(scale);

   double result = 0.;

   for (unsigned gen = 0; gen < 3; ++gen) {
      const double MSI02 = upper_bound.calculate_MSI02(gen);

      result += oneOver16PiSqr * 0.025 * Sqr(g1p) * QHu * QS *
         (0.025 * Sqr(g1p) * QHu * QS * Sqr(vu) * Log(MSI02 /
         Sqr(scale)) - A0(Sqrt(MSI02)));
   }

   return result;
}

BOOST_AUTO_TEST_CASE( test_analytic_inert_singlet_contributions )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double Delta00Pr = 0.;
   double Delta01Pr = 0.;
   double Delta11Pr = 0.;
   for (unsigned gen = 0; gen < 3; ++gen) {
      Delta00Pr += upper_bound.get_unrotated_inert_singlet_contribution(gen, 0, 0);
      Delta01Pr += upper_bound.get_unrotated_inert_singlet_contribution(gen, 0, 1);
      Delta11Pr += upper_bound.get_unrotated_inert_singlet_contribution(gen, 1, 1);
   }

   const double Delta00Pr_expect = get_d2V1lp_inert_singlet_dvd_dvd(model);
   const double Delta01Pr_expect = get_d2V1lp_inert_singlet_dvd_dvu(model);
   const double Delta11Pr_expect = get_d2V1lp_inert_singlet_dvu_dvu(model);

   BOOST_CHECK_CLOSE(Delta00Pr, Delta00Pr_expect, 1.0e-10);
   BOOST_CHECK_CLOSE(Delta01Pr, Delta01Pr_expect, 1.0e-10);
   BOOST_CHECK_CLOSE(Delta11Pr, Delta11Pr_expect, 1.0e-10);
}

double get_d2V1lp_inert_neutral_higgs_dvd_dvd(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   const double gbar = Sqrt(Sqr(g2) + 0.6 * Sqr(g1));

   const double QHd = -3.0;
   const double QHu = -2.0;

   const double scale = model.get_scale();
   A0_fn A0(scale);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double result = 0.;

   for (unsigned gen = 0; gen < 2; ++gen) {
      const Eigen::Array<double,2,1> MHI02(upper_bound.calculate_MHI02(gen));
      const double Sin2ThetaHI0 = upper_bound.calculate_Sin2ThetaHI0(gen);
      const double Cos2ThetaHI0 = upper_bound.calculate_Cos2ThetaHI0(gen);
      const double Sin4ThetaHI0 = 2.0 * Sin2ThetaHI0 * Cos2ThetaHI0;

      const double Lambda = model.get_Lambda12(gen, gen);
      const double inverse_mass_diff = 1.0 / (MHI02(1) - MHI02(0));

      result += 0.5 * oneOver16PiSqr * (0.5 * Sqr(vd) * Log(MHI02(0) * MHI02(1)
         / Power(scale, 4)) * (0.000625 * Power(g1p, 4) * Sqr(QHd) * Sqr(QHd +
         QHu) + Sqr(Lambda * Lambdax * vu * Sin2ThetaHI0 / vd - Cos2ThetaHI0 *
         (0.5 * Sqr(gbar) + 0.025 * Sqr(g1p) * QHd * (QHd - QHu)))) + 0.025 *
         Sqr(g1p) * QHd * (QHd + QHu) * Sqr(vd) * Log(MHI02(1) / MHI02(0)) * (
         Lambda * Lambdax * vu * Sin2ThetaHI0 / vd - Cos2ThetaHI0 * (0.5 * Sqr(
         gbar) + 0.025 * Sqr(g1p) * QHd * (QHd - QHu))) - 0.025 * Sqr(g1p) * QHd
         * (QHd + QHu) * (A0(Sqrt(MHI02(0))) + A0(Sqrt(MHI02(1)))) + (A0(Sqrt(
         MHI02(1))) - A0(Sqrt(MHI02(0)))) * (Cos2ThetaHI0 * (0.5 * Sqr(gbar) +
         0.025 * Sqr(g1p) * QHd * (QHd - QHu)) - Sqr(vd) * inverse_mass_diff *
         (Sqr(Sin2ThetaHI0) * Sqr(0.5 * Sqr(gbar) + 0.025 * Sqr(g1p) * QHd * (
         QHd - QHu)) + Sqr(Lambda) * Sqr(Lambdax) * Sqr(vu) * Sqr(Cos2ThetaHI0)
         / Sqr(vd) + Lambda * Lambdax * vu * Sin4ThetaHI0 * (0.5 * Sqr(gbar) +
         0.025 * Sqr(g1p) * QHd * (QHd - QHu)) / vd)));
   }

   return result;
}

double get_d2V1lp_inert_neutral_higgs_dvd_dvu(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   const double gbar = Sqrt(Sqr(g2) + 0.6 * Sqr(g1));

   const double QHd = -3.0;
   const double QHu = -2.0;

   const double scale  = model.get_scale();
   A0_fn A0(scale);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double result = 0.;

   for (unsigned gen = 0; gen < 2; ++gen) {
      const Eigen::Array<double,2,1> MHI02(upper_bound.calculate_MHI02(gen));
      const double Sin2ThetaHI0 = upper_bound.calculate_Sin2ThetaHI0(gen);
      const double Cos2ThetaHI0 = upper_bound.calculate_Cos2ThetaHI0(gen);
      const double Sin4ThetaHI0 = 2.0 * Sin2ThetaHI0 * Cos2ThetaHI0;

      const double inverse_mass_diff = 1.0 / (MHI02(1) - MHI02(0));

      const double Lambda = model.get_Lambda12(gen, gen);

      result += 0.5 * oneOver16PiSqr * (0.5 * vd * vu * Log(MHI02(0) * MHI02(1)
         / Power(scale, 4)) * (0.000625 * Power(g1p, 4) * QHd * QHu * Sqr(QHd +
         QHu) + (Lambda * Lambdax * vu * Sin2ThetaHI0 / vd - Cos2ThetaHI0 * (0.5
         * Sqr(gbar) + 0.025 * Sqr(g1p) * QHd * (QHd - QHu))) * (Lambda * Lambdax
         * vd * Sin2ThetaHI0 / vu + Cos2ThetaHI0 * (0.5 * Sqr(gbar) - 0.025 * Sqr(
         g1p) * QHu * (QHd - QHu)))) + 0.5 * vd * vu * Log(MHI02(1) / MHI02(0)) *
         (0.025 * Sqr(g1p) * QHd * (QHd + QHu) * (Lambda * Lambdax * vd *
         Sin2ThetaHI0 / vu + Cos2ThetaHI0 * (0.5 * Sqr(gbar) - 0.025 * Sqr(g1p) *
         QHu * (QHd - QHu))) + 0.025 * Sqr(g1p) * QHu * (QHd + QHu) * (Lambda *
         Lambdax * vu * Sin2ThetaHI0 / vd - Cos2ThetaHI0 * (0.5 * Sqr(gbar) + 0.025
         * Sqr(g1p) * QHd * (QHd - QHu)))) - (A0(Sqrt(MHI02(1))) - A0(Sqrt(MHI02(0))))
         * (Lambda * Lambdax * Sin2ThetaHI0 - vd * vu * inverse_mass_diff * (Sqr(
         Sin2ThetaHI0) * (0.5 * Sqr(gbar) + 0.025 * Sqr(g1p) * QHd * (QHd - QHu)) *
         (0.5 * Sqr(gbar) - 0.025 * Sqr(g1p) * QHu * (QHd - QHu)) - Sqr(Lambda) *
         Sqr(Lambdax) * Sqr(Cos2ThetaHI0) - 0.5 * Lambda * Lambdax * Sin4ThetaHI0 *
         (vd * (0.5 * Sqr(gbar) + 0.025 * Sqr(g1p) * QHd * (QHd - QHu)) / vu - vu *
         (0.5 * Sqr(gbar) - 0.025 * Sqr(g1p) * QHu * (QHd - QHu)) / vd))));
   }

   return result;
}

double get_d2V1lp_inert_neutral_higgs_dvu_dvu(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   const double gbar = Sqrt(Sqr(g2) + 0.6 * Sqr(g1));

   const double QHd = -3.0;
   const double QHu = -2.0;

   const double scale  = model.get_scale();
   A0_fn A0(scale);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double result = 0.;

   for (unsigned gen = 0; gen < 2; ++gen) {
      const Eigen::Array<double,2,1> MHI02(upper_bound.calculate_MHI02(gen));
      const double Sin2ThetaHI0 = upper_bound.calculate_Sin2ThetaHI0(gen);
      const double Cos2ThetaHI0 = upper_bound.calculate_Cos2ThetaHI0(gen);
      const double Sin4ThetaHI0 = 2.0 * Sin2ThetaHI0 * Cos2ThetaHI0;

      const double inverse_mass_diff = 1.0 / (MHI02(1) - MHI02(0));

      const double Lambda = model.get_Lambda12(gen, gen);

      result += 0.5 * oneOver16PiSqr * (0.5 * Sqr(vu) * Log(MHI02(0) *
         MHI02(1) / Power(scale, 4)) * (0.000625 * Power(g1p, 4) *
         Sqr(QHu) * Sqr(QHd + QHu) + Sqr(Lambda * Lambdax * vd *
         Sin2ThetaHI0 / vu + Cos2ThetaHI0 * (0.5 * Sqr(gbar) - 0.025 *
         Sqr(g1p) * QHu * (QHd - QHu)))) + 0.025 * Sqr(g1p) * QHu * (QHd
         + QHu) * Sqr(vu) * Log(MHI02(1) / MHI02(0)) * (Lambda * Lambdax
         * vd * Sin2ThetaHI0 / vu + Cos2ThetaHI0 * (0.5 * Sqr(gbar) -
         0.025 * Sqr(g1p) * QHu * (QHd - QHu))) - 0.025 * Sqr(g1p) * QHu
         * (QHd + QHu) * (A0(Sqrt(MHI02(0))) + A0(Sqrt(MHI02(1)))) - (A0(
         Sqrt(MHI02(1))) - A0(Sqrt(MHI02(0)))) * (Cos2ThetaHI0 * (0.5 *
         Sqr(gbar) - 0.025 * Sqr(g1p) * QHu * (QHd - QHu)) + Sqr(vu) *
         inverse_mass_diff * (Sqr(Sin2ThetaHI0) * Sqr(0.5 * Sqr(gbar) -
         0.025 * Sqr(g1p) * QHu * (QHd - QHu)) + Sqr(Lambda) *
         Sqr(Lambdax) * Sqr(vd) * Sqr(Cos2ThetaHI0) / Sqr(vu) - Lambda *
         Lambdax * vd * Sin4ThetaHI0 * (0.5 * Sqr(gbar) - 0.025 * Sqr(g1p)
         * QHu * (QHd - QHu)) / vu)));
   }

   return result;
}

BOOST_AUTO_TEST_CASE( test_analytic_inert_neutral_higgs_contributions )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double Delta00Pr = 0.;
   double Delta01Pr = 0.;
   double Delta11Pr = 0.;
   for (unsigned gen = 0; gen < 2; ++gen) {
      Delta00Pr += upper_bound.get_unrotated_inert_neutral_higgs_contribution(gen, 0, 0);
      Delta01Pr += upper_bound.get_unrotated_inert_neutral_higgs_contribution(gen, 0, 1);
      Delta11Pr += upper_bound.get_unrotated_inert_neutral_higgs_contribution(gen, 1, 1);
   }

   const double Delta00Pr_expect = get_d2V1lp_inert_neutral_higgs_dvd_dvd(model);
   const double Delta01Pr_expect = get_d2V1lp_inert_neutral_higgs_dvd_dvu(model);
   const double Delta11Pr_expect = get_d2V1lp_inert_neutral_higgs_dvu_dvu(model);

   BOOST_CHECK_CLOSE(Delta00Pr, Delta00Pr_expect, 1.0e-10);
   BOOST_CHECK_CLOSE(Delta01Pr, Delta01Pr_expect, 1.0e-10);
   BOOST_CHECK_CLOSE(Delta11Pr, Delta11Pr_expect, 1.0e-10);
}

double get_d2V1lp_inert_charged_higgs_dvd_dvd(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   const double QHd = -3.0;
   const double QHu = -2.0;

   const double scale = model.get_scale();
   A0_fn A0(scale);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double result = 0.;

   for (unsigned gen = 0; gen < 2; ++gen) {
      const Eigen::Array<double,2,1> MHIPM2(upper_bound.calculate_MHIPM2(gen));
      const double Sin2ThetaHIPM = upper_bound.calculate_Sin2ThetaHIPM(gen);
      const double Cos2ThetaHIPM = upper_bound.calculate_Cos2ThetaHIPM(gen);
      const double Sin4ThetaHIPM = 2.0 * Sin2ThetaHIPM * Cos2ThetaHIPM;

      const double Lambda = model.get_Lambda12(gen, gen);
      const double inverse_mass_diff = 1.0 / (MHIPM2(1) - MHIPM2(0));

      result += 0.5 * oneOver16PiSqr * (0.5 * Sqr(vd) * Log(MHIPM2(0) * MHIPM2(1)
         / Power(scale, 4)) * (0.000625 * Power(g1p, 4) * Sqr(QHd) * Sqr(QHd + QHu)
         + Sqr(Lambda * Lambdax * vu * Sin2ThetaHIPM / vd - Cos2ThetaHIPM * (0.5
         * (Sqr(g2) - 0.6 * Sqr(g1)) - 0.025 * Sqr(g1p) * QHd * (QHd - QHu)))) -
         0.025 * Sqr(g1p) * QHd * (QHd + QHu) * Sqr(vd) * Log(MHIPM2(1) / MHIPM2(0))
         * (Lambda * Lambdax * vu * Sin2ThetaHIPM / vd - Cos2ThetaHIPM * (0.5 * (
         Sqr(g2) - 0.6 * Sqr(g1)) - 0.025 * Sqr(g1p) * QHd * (QHd - QHu))) - 0.025 *
         Sqr(g1p) * QHd * (QHd + QHu) * (A0(Sqrt(MHIPM2(0))) + A0(Sqrt(MHIPM2(1))))
         - (A0(Sqrt(MHIPM2(1))) - A0(Sqrt(MHIPM2(0)))) * (Cos2ThetaHIPM * (0.5 * (
         Sqr(g2) - 0.6 * Sqr(g1)) - 0.025 * Sqr(g1p) * QHd * (QHd - QHu)) + Sqr(vd)
         * inverse_mass_diff * (Sqr(Sin2ThetaHIPM) * Sqr(0.5 * (Sqr(g2) - 0.6 * Sqr(
         g1)) - 0.025 * Sqr(g1p) * QHd * (QHd - QHu)) + Sqr(Lambda) * Sqr(Lambdax) *
         Sqr(vu) * Sqr(Cos2ThetaHIPM) / Sqr(vd) + Lambda * Lambdax * vu *
         Sin4ThetaHIPM * (0.5 * (Sqr(g2) - 0.6 * Sqr(g1)) - 0.025 * Sqr(g1p) * QHd *
         (QHd - QHu)) / vd)));
   }

   return result;
}

double get_d2V1lp_inert_charged_higgs_dvd_dvu(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   const double QHd = -3.0;
   const double QHu = -2.0;

   const double scale  = model.get_scale();
   A0_fn A0(scale);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double result = 0.;

   for (unsigned gen = 0; gen < 2; ++gen) {
      const Eigen::Array<double,2,1> MHIPM2(upper_bound.calculate_MHIPM2(gen));
      const double Sin2ThetaHIPM = upper_bound.calculate_Sin2ThetaHIPM(gen);
      const double Cos2ThetaHIPM = upper_bound.calculate_Cos2ThetaHIPM(gen);
      const double Sin4ThetaHIPM = 2.0 * Sin2ThetaHIPM * Cos2ThetaHIPM;

      const double inverse_mass_diff = 1.0 / (MHIPM2(1) - MHIPM2(0));

      const double Lambda = model.get_Lambda12(gen, gen);

      result += 0.5 * oneOver16PiSqr * (0.5 * vd * vu * Log(MHIPM2(0) * MHIPM2(1)
         / Power(scale,4)) * (0.000625 * Power(g1p,4) * QHd * QHu * Sqr(QHd + QHu)
         + (Lambda * Lambdax * vu * Sin2ThetaHIPM / vd - Cos2ThetaHIPM * (0.5 * (
         Sqr(g2) - 0.6 * Sqr(g1)) - 0.025 * Sqr(g1p) * QHd * (QHd - QHu))) * (Lambda
         * Lambdax * vd * Sin2ThetaHIPM / vu + Cos2ThetaHIPM * (0.5 * (Sqr(g2) - 0.6
         * Sqr(g1)) + 0.025 * Sqr(g1p) * QHu * (QHd - QHu)))) - 0.5 * vd * vu * Log(
         MHIPM2(1) / MHIPM2(0)) * (0.025 * Sqr(g1p) * QHd * (QHd + QHu) * (Lambda *
         Lambdax * vd * Sin2ThetaHIPM / vu + Cos2ThetaHIPM * (0.5 * (Sqr(g2) - 0.6 *
         Sqr(g1)) + 0.025 * Sqr(g1p) * QHu * (QHd - QHu))) + 0.025 * Sqr(g1p) * QHu
         * (QHd + QHu) * (Lambda * Lambdax * vu * Sin2ThetaHIPM / vd - Cos2ThetaHIPM
         * (0.5 * (Sqr(g2) - 0.6 * Sqr(g1)) - 0.025 * Sqr(g1p) * QHd * (QHd - QHu))))
         + (A0(Sqrt(MHIPM2(1))) - A0(Sqrt(MHIPM2(0)))) * (Lambda * Lambdax *
         Sin2ThetaHIPM + vd * vu * inverse_mass_diff * (Sqr(Sin2ThetaHIPM) * (0.5 *
         (Sqr(g2) - 0.6 * Sqr(g1)) - 0.025 * Sqr(g1p) * QHd * (QHd - QHu)) * (0.5 *
         (Sqr(g2) - 0.6 * Sqr(g1)) + 0.025 * Sqr(g1p) * QHu * (QHd - QHu)) - Sqr(
         Lambda) * Sqr(Lambdax) * Sqr(Cos2ThetaHIPM) - 0.5 * Lambda * Lambdax *
         Sin4ThetaHIPM * (-vu * (0.5 * (Sqr(g2) - 0.6 * Sqr(g1)) + 0.025 * Sqr(g1p)
         * QHu * (QHd - QHu)) / vd + vd * (0.5 * (Sqr(g2) - 0.6 * Sqr(g1)) - 0.025 *
         Sqr(g1p) * QHd * (QHd - QHu)) / vu))));  
   }

   return result;
}

double get_d2V1lp_inert_charged_higgs_dvu_dvu(const CNE6SSM_mass_eigenstates& model)
{
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   const double QHd = -3.0;
   const double QHu = -2.0;

   const double scale  = model.get_scale();
   A0_fn A0(scale);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double result = 0.;

   for (unsigned gen = 0; gen < 2; ++gen) {
      const Eigen::Array<double,2,1> MHIPM2(upper_bound.calculate_MHIPM2(gen));
      const double Sin2ThetaHIPM = upper_bound.calculate_Sin2ThetaHIPM(gen);
      const double Cos2ThetaHIPM = upper_bound.calculate_Cos2ThetaHIPM(gen);
      const double Sin4ThetaHIPM = 2.0 * Sin2ThetaHIPM * Cos2ThetaHIPM;

      const double inverse_mass_diff = 1.0 / (MHIPM2(1) - MHIPM2(0));

      const double Lambda = model.get_Lambda12(gen, gen);

      result += 0.5 * oneOver16PiSqr * (0.5 * Sqr(vu) * Log(MHIPM2(0) * MHIPM2(1)
         / Power(scale, 4)) * (0.000625 * Power(g1p, 4) * Sqr(QHu) * Sqr(QHd +
         QHu) + Sqr(Lambda * Lambdax * vd * Sin2ThetaHIPM / vu + Cos2ThetaHIPM *
         (0.5 * (Sqr(g2) - 0.6 * Sqr(g1)) + 0.025 * Sqr(g1p) * QHu * (QHd - QHu))))
         - 0.025 * Sqr(g1p) * QHu * (QHd + QHu) * Sqr(vu) * Log(MHIPM2(1) / MHIPM2(0))
         * (Lambda * Lambdax * vd * Sin2ThetaHIPM / vu + Cos2ThetaHIPM * (0.5 *
         (Sqr(g2) - 0.6 * Sqr(g1)) + 0.025 * Sqr(g1p) * QHu * (QHd - QHu))) - 0.025 *
         Sqr(g1p) * QHu * (QHd + QHu) * (A0(Sqrt(MHIPM2(0))) + A0(Sqrt(MHIPM2(1)))) +
         (A0(Sqrt(MHIPM2(1))) - A0(Sqrt(MHIPM2(0)))) * (Cos2ThetaHIPM * (0.5 * (Sqr(g2)
         - 0.6 * Sqr(g1)) + 0.025 * Sqr(g1p) * QHu * (QHd - QHu)) - Sqr(vu) *
         inverse_mass_diff * (Sqr(Sin2ThetaHIPM) * Sqr(0.5 * (Sqr(g2) - 0.6 * Sqr(g1))
         + 0.025 * Sqr(g1p) * QHu * (QHd - QHu)) + Sqr(Lambda) * Sqr(Lambdax) * Sqr(vd)
         * Sqr(Cos2ThetaHIPM) / Sqr(vu) - Lambda * Lambdax * vd * Sin4ThetaHIPM * (0.5
         * (Sqr(g2) - 0.6 * Sqr(g1)) + 0.025 * Sqr(g1p) * QHu * (QHd - QHu)) / vu )));
   }

   return result;
}

BOOST_AUTO_TEST_CASE( test_analytic_inert_charged_contributions )
{
   CNE6SSM_mass_eigenstates model;

   set_test_model_parameters(model);

   CNE6SSM_higgs_upper_bound upper_bound(model);

   double Delta00Pr = 0.;
   double Delta01Pr = 0.;
   double Delta11Pr = 0.;
   for (unsigned gen = 0; gen < 2; ++gen) {
      Delta00Pr += upper_bound.get_unrotated_inert_charged_higgs_contribution(gen, 0, 0);
      Delta01Pr += upper_bound.get_unrotated_inert_charged_higgs_contribution(gen, 0, 1);
      Delta11Pr += upper_bound.get_unrotated_inert_charged_higgs_contribution(gen, 1, 1);
   }

   const double Delta00Pr_expect = get_d2V1lp_inert_charged_higgs_dvd_dvd(model);
   const double Delta01Pr_expect = get_d2V1lp_inert_charged_higgs_dvd_dvu(model);
   const double Delta11Pr_expect = get_d2V1lp_inert_charged_higgs_dvu_dvu(model);

   BOOST_CHECK_CLOSE(Delta00Pr, Delta00Pr_expect, 1.0e-10);
   BOOST_CHECK_CLOSE(Delta01Pr, Delta01Pr_expect, 1.0e-10);
   BOOST_CHECK_CLOSE(Delta11Pr, Delta11Pr_expect, 1.0e-10);
}
