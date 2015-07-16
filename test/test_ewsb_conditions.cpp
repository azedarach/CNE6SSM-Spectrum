// ====================================================================
// Test suite for implementation of EWSB conditions
// ====================================================================

#include <limits>
#include <random>
#include <chrono>
#include <Eigen/Core>

#include "wrappers.hpp"
#include "CNE6SSM_two_scale_model.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_ewsb_conditions

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace std;
using namespace flexiblesusy;

// Test that the default EWSB conditions give
// expected values for a set of parameters
BOOST_AUTO_TEST_CASE( test_default_ewsb_tree_level )
{
   CNE6SSM<Two_scale> default_model = CNE6SSM<Two_scale>();

   CNE6SSM_input_parameters<Two_scale> default_inputs;
   default_inputs.QSInput = 5.;
   default_model.set_input_parameters(default_inputs);

   default_model.set_QS(default_inputs.QSInput);
   default_model.set_vd(24.4779);
   default_model.set_vu(244.779);
   default_model.set_vs(5381.382);
   default_model.set_vsb(5919.520587);
   default_model.set_vphi(8000.);
   default_model.set_g1(0.46);
   default_model.set_g1p(0.47);
   default_model.set_g2(0.61);
   default_model.set_Lambdax(0.01);
   default_model.set_Sigmax(0.12);
   default_model.set_TLambdax(1100.);
   default_model.set_TSigmax(540.);
   default_model.set_MuPhi(100.);
   default_model.set_BMuPhi(150.);
   default_model.set_KappaPr(0.5);
   default_model.set_TKappaPr(3000.);
   default_model.set_XiF(0.4);
   default_model.set_LXiF(10.);
   default_model.set_mHd2(2.3e6);
   default_model.set_mHu2(-1.3e6);
   default_model.set_ms2(-1.9e6);
   default_model.set_msbar2(2.1e6);
   default_model.set_mphi2(1.2e5);

   BOOST_CHECK_CLOSE(default_model.get_ewsb_eq_hh_1(), -3.902254749e7, 1.0e-7);
   BOOST_CHECK_CLOSE(default_model.get_ewsb_eq_hh_2(), -1.541987751e6, 1.0e-7);
   BOOST_CHECK_CLOSE(default_model.get_ewsb_eq_hh_3(), -7.156101591e6, 1.0e-7);
   BOOST_CHECK_CLOSE(default_model.get_ewsb_eq_hh_4(), -1.393365501e6, 1.0e-7);
   BOOST_CHECK_CLOSE(default_model.get_ewsb_eq_hh_5(), 3.094450484e7, 1.0e-7);

}

void randomize_model(CNE6SSM<Two_scale>& model)
{
   static default_random_engine generator;
   static uniform_real_distribution<double> distribution(-50., 50.);

   CNE6SSM_input_parameters<Two_scale> random_inputs;
   random_inputs.QSInput = distribution(generator);
   model.set_input_parameters(random_inputs);

   model.set_QS(random_inputs.QSInput);
   model.set_vd(distribution(generator));
   model.set_vu(distribution(generator));
   model.set_vs(distribution(generator));
   model.set_vsb(distribution(generator));
   model.set_vphi(distribution(generator));
   model.set_g1(distribution(generator));
   model.set_g1p(distribution(generator));
   model.set_g2(distribution(generator));
   model.set_Lambdax(distribution(generator));
   model.set_Sigmax(distribution(generator));
   model.set_TLambdax(distribution(generator));
   model.set_TSigmax(distribution(generator));
   model.set_MuPhi(distribution(generator));
   model.set_BMuPhi(distribution(generator));
   model.set_KappaPr(distribution(generator));
   model.set_TKappaPr(distribution(generator));
   model.set_XiF(distribution(generator));
   model.set_LXiF(distribution(generator));
   model.set_mHd2(distribution(generator));
   model.set_mHu2(distribution(generator));
   model.set_ms2(distribution(generator));
   model.set_msbar2(distribution(generator));
   model.set_mphi2(distribution(generator));

}

// Test tree level solution of EWSB equations via soft masses
BOOST_AUTO_TEST_CASE( test_default_zeroes_custom )
{
   CNE6SSM<Two_scale> default_zero_model = CNE6SSM<Two_scale>();

   randomize_model(default_zero_model);

   Eigen::Array<double,5,1> tree_level_masses 
      = default_zero_model.get_ewsb_tree_level_soft_masses();

   default_zero_model.set_mHd2(tree_level_masses(0));
   default_zero_model.set_mHu2(tree_level_masses(1));
   default_zero_model.set_ms2(tree_level_masses(2));
   default_zero_model.set_msbar2(tree_level_masses(3));
   default_zero_model.set_mphi2(tree_level_masses(4));

   BOOST_CHECK_SMALL(default_zero_model.get_ewsb_eq_hh_1(), 1.0e-5);
   BOOST_CHECK_SMALL(default_zero_model.get_ewsb_eq_hh_2(), 1.0e-5);
   BOOST_CHECK_SMALL(default_zero_model.get_ewsb_eq_hh_3(), 1.0e-5);
   BOOST_CHECK_SMALL(default_zero_model.get_ewsb_eq_hh_4(), 1.0e-5);
   BOOST_CHECK_SMALL(default_zero_model.get_ewsb_eq_hh_5(), 1.0e-5);

}

// Test that parameters that zero custom EWSB conditions also zero
// default EWSB conditions
BOOST_AUTO_TEST_CASE( test_custom_zeroes_default )
{
   CNE6SSM<Two_scale> custom_zero_model = CNE6SSM<Two_scale>();

   randomize_model(custom_zero_model);

   double QS = custom_zero_model.get_QS();
   double vd = custom_zero_model.get_vd();
   double vu = custom_zero_model.get_vu();
   double vs = custom_zero_model.get_vs();
   double vsb = custom_zero_model.get_vsb();
   double vphi = custom_zero_model.get_vphi();
   double mHd2 = custom_zero_model.get_mHd2();
   double mHu2 = custom_zero_model.get_mHu2();
   double ms2 = custom_zero_model.get_ms2();
   double msbar2 = custom_zero_model.get_msbar2();
   double mphi2 = custom_zero_model.get_mphi2();
   double g1 = custom_zero_model.get_g1();
   double g1p = custom_zero_model.get_g1p();
   double g2 = custom_zero_model.get_g2();
   double Lambdax = custom_zero_model.get_Lambdax();
   double Sigmax = custom_zero_model.get_Sigmax();
   double TLambdax = custom_zero_model.get_TLambdax();
   double TSigmax = custom_zero_model.get_TSigmax();
   double MuPhi = custom_zero_model.get_MuPhi();
   double BMuPhi = custom_zero_model.get_BMuPhi();
   double XiF = custom_zero_model.get_XiF();
   double LXiF = custom_zero_model.get_LXiF();
   double KappaPr = custom_zero_model.get_KappaPr();
   double TKappaPr = custom_zero_model.get_TKappaPr();

   TLambdax = 
      0.5*(ms2*Sqr(vs) - msbar2*Sqr(vsb) + 0.5*Sqr(vs)*
           AbsSqr(Lambdax)*Sqr(vd) + 0.5*Sqr(vs)*AbsSqr(Lambdax)*Sqr(vu) + 0.5*Sqr(vphi)
           *AbsSqr(Sigmax)*Sqr(vs) - 0.5*Sqr(vphi)*AbsSqr(Sigmax)*Sqr(vsb) - 0.25*vphi*
           vsb*vd*vu*Conj(Sigmax)*Lambdax - 0.25*vphi*vsb*vd*vu*Sigmax*Conj(Lambdax) -
           0.0375*QS*Sqr(g1p)*Sqr(vs)*Sqr(vd) - 0.025*QS*Sqr(g1p)*Sqr(vs)*Sqr(vu) + 
           0.0125*Sqr(QS)*Sqr(g1p)*Power(vs,4) - 0.0375*QS*Sqr(g1p)*Sqr(vsb)*Sqr(vd) - 
           0.025*QS*Sqr(g1p)*Sqr(vsb)*Sqr(vu) - 0.0125*Sqr(QS)*Sqr(g1p)*Power(vsb,4))
      / (0.35355339059327373*vd*vu*vs);

   mHd2 = (0.0125*(28.284271247461902*vs*vu*Conj(TLambdax) -
                   20*vphi*vsb*vu*Conj(Sigmax)*Lambdax - 20*vphi*vsb*vu*Conj(Lambdax)*Sigmax -
                   6*Power(vd,3)*Sqr(g1) - 9*Power(vd,3)*Sqr(g1p) - 10*Power(vd,3)*Sqr(g2) - 40
                   *vd*AbsSqr(Lambdax)*Sqr(vs) + 3*QS*vd*Sqr(g1p)*Sqr(vs) - 3*QS*vd*Sqr(g1p)*
                   Sqr(vsb) - 40*vd*AbsSqr(Lambdax)*Sqr(vu) + 6*vd*Sqr(g1)*Sqr(vu) - 6*vd*Sqr(
                      g1p)*Sqr(vu) + 10*vd*Sqr(g2)*Sqr(vu) + 28.284271247461902*vs*vu*TLambdax))/vd;

   mHu2 = (mHd2*Sqr(vd) + 0.5*Sqr(vs)*AbsSqr(Lambdax)*
      Sqr(vd) - 0.5*Sqr(vs)*AbsSqr(Lambdax)*Sqr(vu) + 0.125*Sqr(g2)*Power(vd,4)
      + 0.075*Sqr(g1)*Power(vd,4) - 0.125*Sqr(g2)*Power(vu,4) - 0.075*Sqr(g1)*
      Power(vu,4) + 0.1125*Sqr(g1p)*Power(vd,4) - 0.05*Sqr(g1p)*Power(vu,4)
      - 0.0375*QS*Sqr(g1p)*Sqr(vd)*Sqr(vs) + 0.0375*QS*Sqr(g1p)*Sqr(vd)*Sqr(vsb)
           + 0.025*QS*Sqr(g1p)*Sqr(vu)*Sqr(vs) - 0.025*QS*Sqr(g1p)*Sqr(vu)*Sqr(vsb)) / Sqr(vu);

   XiF = (msbar2*vsb - 0.35355339059327373*MuPhi*vphi*vs*Conj(Sigmax)
          - 0.35355339059327373*vphi*vs*Conj(TSigmax) + 0.25*vd*vphi*vu*Conj(Sigmax)*
          Lambdax - 0.35355339059327373*vphi*vs*Conj(MuPhi)*
          Sigmax + 0.25*vd*vphi*vu*Conj(Lambdax)*Sigmax +
          0.0125*Power(vsb,3)*Sqr(g1p)*Sqr(QS) + 0.0375*QS*vsb*Sqr(g1p)*Sqr(vd) + 0.5*
          vsb*AbsSqr(Sigmax)*Sqr(vphi) - 0.25*vs*Conj(Sigmax)*KappaPr*Sqr(vphi) - 0.25
          *vs*Conj(KappaPr)*Sigmax*Sqr(vphi) + 0.5*vsb*AbsSqr(Sigmax)*Sqr(vs) - 0.0125
          *vsb*Sqr(g1p)*Sqr(QS)*Sqr(vs) + 0.025*QS*vsb*Sqr(g1p)*Sqr(vu) -
          0.35355339059327373*vphi*vs*TSigmax) / (Sigmax*vs);

   BMuPhi = -0.5*(mphi2*vphi + vphi*AbsSqr(MuPhi) + Power(vphi,3)*AbsSqr(
                     KappaPr) + 0.7071067811865475*
                  MuPhi*Conj(XiF) - 0.35355339059327373*MuPhi*vs*vsb*Conj(Sigmax) +
                  0.7071067811865475*Conj(LXiF) - 0.35355339059327373*vs*vsb*Conj(TSigmax) +
                  vphi*Conj(XiF)*KappaPr - 0.5*vphi*vs*vsb*Conj(Sigmax)*KappaPr + 0.25*vd*vsb*
                  vu*Conj(Sigmax)*Lambdax + 0.7071067811865475*Conj(MuPhi)*XiF + vphi*Conj(
                     KappaPr)*XiF - 0.35355339059327373*vs*vsb*Conj(MuPhi)*Sigmax - 0.5*vphi*vs*
                  vsb*Conj(KappaPr)*Sigmax + 0.25*vd*vsb*vu*Conj(Lambdax)*Sigmax +
                  0.7071067811865475*LXiF + 1.0606601717798212*MuPhi*Conj(KappaPr)*Sqr(vphi) +
                  0.35355339059327373*Conj(TKappaPr)*Sqr(vphi) + 1.0606601717798212*Conj(
                     MuPhi)*KappaPr*Sqr(vphi) + 0.5*vphi*AbsSqr(Sigmax)*Sqr(vs) + 0.5*vphi*AbsSqr
                  (Sigmax)*Sqr(vsb) + 0.35355339059327373*Sqr(vphi)*TKappaPr -
                  0.35355339059327373*vs*vsb*TSigmax) / (0.5*vphi);

   custom_zero_model.set_mHu2(mHu2);
   custom_zero_model.set_TLambdax(TLambdax);
   custom_zero_model.set_BMuPhi(BMuPhi);
   custom_zero_model.set_XiF(XiF);
   custom_zero_model.set_mHd2(mHd2);

   BOOST_CHECK_SMALL(custom_zero_model.get_ewsb_eq_hh_1(), 1.0e-5);
   BOOST_CHECK_SMALL(custom_zero_model.get_ewsb_eq_hh_2(), 1.0e-5);
   BOOST_CHECK_SMALL(custom_zero_model.get_ewsb_eq_hh_3(), 1.0e-5);
   BOOST_CHECK_SMALL(custom_zero_model.get_ewsb_eq_hh_4(), 1.0e-5);
   BOOST_CHECK_SMALL(custom_zero_model.get_ewsb_eq_hh_5(), 1.0e-5);
}
