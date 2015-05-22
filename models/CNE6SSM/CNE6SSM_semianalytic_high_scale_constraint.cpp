// ====================================================================
// Test implementation of high scale constraint for semianalytic
// algorithm
// ====================================================================

#include "CNE6SSM_semianalytic_high_scale_constraint.hpp"
#include "CNE6SSM_semianalytic_model.hpp"
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

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define SM(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME CNE6SSM<Semianalytic>

CNE6SSM_high_scale_constraint<Semianalytic>::CNE6SSM_high_scale_constraint()
   : Constraint<Semianalytic>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
{
}

CNE6SSM_high_scale_constraint<Semianalytic>::CNE6SSM_high_scale_constraint(
   CNE6SSM<Semianalytic>* model_)
   : Constraint<Semianalytic>()
   , model(model_)
{
   initialize();
}

CNE6SSM_high_scale_constraint<Semianalytic>::~CNE6SSM_high_scale_constraint()
{
}

void CNE6SSM_high_scale_constraint<Semianalytic>::apply()
{
   assert(model && "Error: CNE6SSM_high_scale_constraint::apply():"
          " model pointer must not be zero");

   if (std::fabs(model->get_g1()) > 3.54491) {
#ifdef ENABLE_VERBOSE
      ERROR("CNE6SSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g1 = " << model->get_g1());
#endif
      model->set_g1(3.54491);
   }
   if (std::fabs(model->get_g2()) > 3.54491) {
#ifdef ENABLE_VERBOSE
      ERROR("CNE6SSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g2 = " << model->get_g2());
#endif
      model->set_g2(3.54491);
   }
   if (std::fabs(model->get_g3()) > 3.54491) {
#ifdef ENABLE_VERBOSE
      ERROR("CNE6SSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g3 = " << model->get_g3());
#endif
      model->set_g3(3.54491);
   }
   // Yukawa couplings should probably be reset, as in softsusy
//    if (MaxAbsValue(model->get_Yu()) > 3.0) {
// #ifdef ENABLE_VERBOSE
//       ERROR("CNE6SSM_high_scale_constraint: Non-perturbative gauge "
//             "coupling Yu = " << model->get_Yu());
// #endif
//       for (std::size_t i = 0; i < 3; ++i) {
//          if (model->get_Yu(i,i) > 3.0)
//             model->set_Yu(i,i,1.0);
//       }
//    }
//    if (MaxAbsValue(model->get_Yd()) > 3.0) {
// #ifdef ENABLE_VERBOSE
//       ERROR("CNE6SSM_high_scale_constraint: Non-perturbative gauge "
//             "coupling Yd = " << model->get_Yd());
// #endif
//       for (std::size_t i = 0; i < 3; ++i) {
//          if (model->get_Yd(i,i) > 3.0)
//             model->set_Yd(i,i,1.0);
//       }
//    }
//    if (MaxAbsValue(model->get_Ye()) > 3.0) {
// #ifdef ENABLE_VERBOSE
//       ERROR("CNE6SSM_high_scale_constraint: Non-perturbative gauge "
//             "coupling Ye = " << model->get_Yu());
// #endif
//       for (std::size_t i = 0; i < 3; ++i) {
//          if (model->get_Ye(i,i) > 3.0)
//             model->set_Ye(i,i,1.0);
//       }
//    }

   update_scale();

   const auto MuPrInput = INPUTPARAMETER(MuPrInput);
   const auto MuPhiInput = INPUTPARAMETER(MuPhiInput);
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
   const auto g1 = MODELPARAMETER(g1);

   MODEL->set_g1p(g1);
   MODEL->set_MuPr(MuPrInput);
   MODEL->set_MuPhi(MuPhiInput);
   MODEL->set_SigmaL(SigmaLInput);
   MODEL->set_KappaPr(KappaPrInput);
   MODEL->set_Sigmax(SigmaxInput);
   MODEL->set_gD(gDInput);
   MODEL->set_hE(hEInput);
   MODEL->set_Kappa(KappaInput);
   MODEL->set_Lambda12(Lambda12Input);
   MODEL->set_Lambdax(LambdaxInput);
   MODEL->set_fu(fuInput);
   MODEL->set_fd(fdInput);

   {
      const auto g1 = MODELPARAMETER(g1);
      const auto g2 = MODELPARAMETER(g2);
      const auto g3 = MODELPARAMETER(g3);
      const auto g1p = MODELPARAMETER(g1p);
      const auto Yd = MODELPARAMETER(Yd);
      const auto hE = MODELPARAMETER(hE);
      const auto Ye = MODELPARAMETER(Ye);
      const auto SigmaL = MODELPARAMETER(SigmaL);
      const auto KappaPr = MODELPARAMETER(KappaPr);
      const auto Sigmax = MODELPARAMETER(Sigmax);
      const auto gD = MODELPARAMETER(gD);
      const auto Kappa = MODELPARAMETER(Kappa);
      const auto Lambda12 = MODELPARAMETER(Lambda12);
      const auto Lambdax = MODELPARAMETER(Lambdax);
      const auto fu = MODELPARAMETER(fu);
      const auto fd = MODELPARAMETER(fd);
      const auto Yu = MODELPARAMETER(Yu);

      if (MaxAbsValue(g1) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("g1", MaxAbsValue(g1), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("g1");
      if (MaxAbsValue(g2) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("g2", MaxAbsValue(g2), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("g2");
      if (MaxAbsValue(g3) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("g3", MaxAbsValue(g3), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("g3");
      if (MaxAbsValue(g1p) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("g1p", MaxAbsValue(g1p), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("g1p");
      if (MaxAbsValue(Yd) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Yd", MaxAbsValue(Yd), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Yd");
      if (MaxAbsValue(hE) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("hE", MaxAbsValue(hE), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("hE");
      if (MaxAbsValue(Ye) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Ye", MaxAbsValue(Ye), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Ye");
      if (MaxAbsValue(SigmaL) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("SigmaL", MaxAbsValue(SigmaL), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("SigmaL");
      if (MaxAbsValue(KappaPr) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("KappaPr", MaxAbsValue(KappaPr), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("KappaPr");
      if (MaxAbsValue(Sigmax) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Sigmax", MaxAbsValue(Sigmax), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Sigmax");
      if (MaxAbsValue(gD) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("gD", MaxAbsValue(gD), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("gD");
      if (MaxAbsValue(Kappa) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Kappa", MaxAbsValue(Kappa), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Kappa");
      if (MaxAbsValue(Lambda12) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Lambda12", MaxAbsValue(Lambda12), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Lambda12");
      if (MaxAbsValue(Lambdax) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Lambdax", MaxAbsValue(Lambdax), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Lambdax");
      if (MaxAbsValue(fu) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("fu", MaxAbsValue(fu), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("fu");
      if (MaxAbsValue(fd) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("fd", MaxAbsValue(fd), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("fd");
      if (MaxAbsValue(Yu) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter("Yu", MaxAbsValue(Yu), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter("Yu");

   }
}

double CNE6SSM_high_scale_constraint<Semianalytic>::get_scale() const
{
   return scale;
}

double CNE6SSM_high_scale_constraint<Semianalytic>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const CNE6SSM_input_parameters<Semianalytic>& CNE6SSM_high_scale_constraint<Semianalytic>::get_input_parameters() const
{
   return model->get_input();
}

CNE6SSM<Semianalytic>* CNE6SSM_high_scale_constraint<Semianalytic>::get_model() const
{
   return model;
}

void CNE6SSM_high_scale_constraint<Semianalytic>::set_model(Semianalytic_model* model_)
{
   model = cast_model<CNE6SSM<Semianalytic>*>(model_);
}

void CNE6SSM_high_scale_constraint<Semianalytic>::set_scale(double s)
{
   scale = s;
}

void CNE6SSM_high_scale_constraint<Semianalytic>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void CNE6SSM_high_scale_constraint<Semianalytic>::initialize()
{
   assert(model && "CNE6SSM_high_scale_constraint<Semianalytic>::"
          "initialize(): model pointer is zero.");

   initial_scale_guess = 2.e16;

   scale = initial_scale_guess;
}

void CNE6SSM_high_scale_constraint<Semianalytic>::update_scale()
{
   assert(model && "CNE6SSM_high_scale_constraint<Semianalytic>::"
          "update_scale(): model pointer is zero.");

   const double currentScale = model->get_scale();
   const CNE6SSM_soft_parameters beta_functions(model->calc_beta());

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto beta_g1 = BETAPARAMETER(g1);
   const auto beta_g2 = BETAPARAMETER(g2);

   scale = currentScale*exp((-g1 + g2)/(BETA(g1) - BETA(g2)));


   if (errno == ERANGE) {
#ifdef ENABLE_VERBOSE
      ERROR("CNE6SSM_high_scale_constraint<Semianalytic>: Overflow error"
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
