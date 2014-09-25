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

// File generated at Mon 15 Sep 2014 17:29:10

#include "CNE6SSM_two_scale_high_scale_constraint.hpp"
#include "CNE6SSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "numerics.hpp"

#include <cassert>
#include <cmath>
#include <cerrno>
#include <cstring>

namespace flexiblesusy {

#define INPUTPARAMETER(p) inputPars.p
#define MODELPARAMETER(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define SM(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME CNE6SSM<Two_scale>

CNE6SSM_high_scale_constraint<Two_scale>::CNE6SSM_high_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , fixed_scale(0.)
   , model(0)
   , inputPars()
{
}

CNE6SSM_high_scale_constraint<Two_scale>::CNE6SSM_high_scale_constraint(const CNE6SSM_input_parameters& inputPars_)
   : Constraint<Two_scale>()
   , fixed_scale(0.)
   , model(0)
   , inputPars(inputPars_)
{
   initialize();
}

CNE6SSM_high_scale_constraint<Two_scale>::~CNE6SSM_high_scale_constraint()
{
}

void CNE6SSM_high_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: CNE6SSM_high_scale_constraint:"
          " model pointer must not be zero");

   if (std::fabs(model->get_g1()) > 3.0) {
#ifdef ENABLE_VERBOSE
      ERROR("CNE6SSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g1 = " << model->get_g1());
#endif
      model->set_g1(1.0);
   }
   if (std::fabs(model->get_g2()) > 3.0) {
#ifdef ENABLE_VERBOSE
      ERROR("CNE6SSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g2 = " << model->get_g2());
#endif
      model->set_g2(1.0);
   }
   if (std::fabs(model->get_g3()) > 3.0) {
#ifdef ENABLE_VERBOSE
      ERROR("CNE6SSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g3 = " << model->get_g3());
#endif
      model->set_g3(1.0);
   }

   update_scale();

   const auto MuPrInput = INPUTPARAMETER(MuPrInput);
   const auto BMuPrInput = INPUTPARAMETER(BMuPrInput);
   const auto MuPhiInput = INPUTPARAMETER(MuPhiInput);
   const auto BMuPhiInput = INPUTPARAMETER(BMuPhiInput);
   const auto SigmaLInput = INPUTPARAMETER(SigmaLInput);
   const auto KappaPrInput = INPUTPARAMETER(KappaPrInput);
   const auto SigmaxInput = INPUTPARAMETER(SigmaxInput);
   const auto gDInput = INPUTPARAMETER(gDInput);
   const auto hEInput = INPUTPARAMETER(hEInput);
   const auto KappaInput = INPUTPARAMETER(KappaInput);
   const auto Lambda12Input = INPUTPARAMETER(Lambda12Input);
   const auto LambdaxInput = INPUTPARAMETER(LambdaxInput);
   const auto fuInput = INPUTPARAMETER(fuInput);
   const auto fdInput = INPUTPARAMETER(fdInput);
   const auto Azero = INPUTPARAMETER(Azero);
   const auto m0 = INPUTPARAMETER(m0);
   const auto m12 = INPUTPARAMETER(m12);
   const auto g1 = MODELPARAMETER(g1);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto KappaPr = MODELPARAMETER(KappaPr);
   const auto Sigmax = MODELPARAMETER(Sigmax);
   const auto hE = MODELPARAMETER(hE);
   const auto SigmaL = MODELPARAMETER(SigmaL);
   const auto gD = MODELPARAMETER(gD);
   const auto fu = MODELPARAMETER(fu);
   const auto fd = MODELPARAMETER(fd);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   MODEL->set_g1p(g1);
   MODEL->set_MuPr(MuPrInput);
   MODEL->set_BMuPr(BMuPrInput);
   MODEL->set_MuPhi(MuPhiInput);
   MODEL->set_BMuPhi(BMuPhiInput);
   MODEL->set_SigmaL(SigmaLInput);
   MODEL->set_KappaPr(KappaPrInput);
   MODEL->set_Sigmax(SigmaxInput);
   MODEL->set_gD(gDInput);
   MODEL->set_hE(hEInput);
   MODEL->set_Kappa(KappaInput);
   MODEL->set_Lambda12(Lambda12Input);
   //MODEL->set_Lambdax(LambdaxInput);
   MODEL->set_fu(fuInput);
   MODEL->set_fd(fdInput);
   MODEL->set_TYe(Azero*Ye);
   MODEL->set_TYd(Azero*Yd);
   MODEL->set_TYu(Azero*Yu);
   MODEL->set_TKappaPr(Azero*KappaPr);
   MODEL->set_TSigmax(Azero*Sigmax);
   MODEL->set_ThE(Azero*hE);
   MODEL->set_TSigmaL(Azero*SigmaL);
   MODEL->set_TgD(Azero*gD);
   MODEL->set_Tfu(Azero*fu);
   MODEL->set_Tfd(Azero*fd);
   MODEL->set_TKappa(Azero*Kappa);
   MODEL->set_TLambda12(Azero*Lambda12);
   MODEL->set_TLambdax(Azero*Lambdax);
   MODEL->set_mHd2(Sqr(m0));
   MODEL->set_mHu2(Sqr(m0));
   MODEL->set_ms2(Sqr(m0));
   MODEL->set_msbar2(Sqr(m0));
   MODEL->set_mphi2(Sqr(m0));
   MODEL->set_mHp2(Sqr(m0));
   MODEL->set_mHpbar2(Sqr(m0));
   MODEL->set_mH1I2(Sqr(m0)*UNITMATRIX(2));
   MODEL->set_mH2I2(Sqr(m0)*UNITMATRIX(2));
   MODEL->set_mSI2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_mq2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_ml2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_md2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_mu2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_me2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_mDx2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_mDxbar2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_MassB(m12);
   MODEL->set_MassWB(m12);
   MODEL->set_MassG(m12);
   MODEL->set_MassBp(m12);
   std::cout << "Q = " << model->get_scale() << ", ";
   std::cout << "g1 = " << model->get_g1() << ", ";
   std::cout << "g1p = " << model->get_g1p() << ", ";
   std::cout << "g2 = " << model->get_g2() << ", ";
   std::cout << "g3 = " << model->get_g3() << ", ";
   std::cout << "yt = " << model->get_Yu(2,2) << ", ";
   std::cout << "yb = " << model->get_Yd(2,2) << ", ";
   std::cout << "ytau = " << model->get_Ye(2,2) << ", ";
   std::cout << "vd = " << model->get_vd() << ", ";
   std::cout << "vu = " << model->get_vu() << ", ";
   std::cout << "vs = " << model->get_vs() << ", ";
   std::cout << "vsb = " << model->get_vsb() << ", ";
   std::cout << "vphi = " << model->get_vphi() << "\n";
}

double CNE6SSM_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double CNE6SSM_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void CNE6SSM_high_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<CNE6SSM<Two_scale> >(model_);
}

void CNE6SSM_high_scale_constraint<Two_scale>::set_input_parameters(const CNE6SSM_input_parameters& inputPars_)
{
   inputPars = inputPars_;
}

void CNE6SSM_high_scale_constraint<Two_scale>::set_scale(double s)
{
   fixed_scale = s;
}

void CNE6SSM_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   fixed_scale = 0.;
   model = NULL;
}

void CNE6SSM_high_scale_constraint<Two_scale>::initialize()
{
   initial_scale_guess = 2.e16;

   scale = initial_scale_guess;
}

void CNE6SSM_high_scale_constraint<Two_scale>::update_scale()
{
   if (!is_zero(fixed_scale)) {
      scale = fixed_scale;
      return;
   }

   const double currentScale = model->get_scale();
   const CNE6SSM_soft_parameters beta_functions(model->calc_beta());

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto beta_g1 = BETAPARAMETER(g1);
   const auto beta_g2 = BETAPARAMETER(g2);

   scale = currentScale*exp((-g1 + g2)/(BETA(g1) - BETA(g2)));


   if (errno == ERANGE) {
#ifdef ENABLE_VERBOSE
      ERROR("CNE6SSM_high_scale_constraint<Two_scale>: Overflow error"
            " during calculation of high scale: " << strerror(errno) << '\n'
            << "   current scale = " << currentScale << '\n'
            << "   new scale = " << scale << '\n'
            << "   resetting scale to " << get_initial_scale_guess());
#endif
      scale = get_initial_scale_guess();
      errno = 0;
   }


}

} // namespace flexiblesusy
