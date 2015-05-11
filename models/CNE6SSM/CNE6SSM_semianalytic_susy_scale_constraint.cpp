// ====================================================================
// Test implementation of SUSY scale constraint for semianalytic
// algorithm
// ====================================================================

#include "CNE6SSM_semianalytic_susy_scale_constraint.hpp"
#include "CNE6SSM_semianalytic_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"

#include <cassert>
#include <cmath>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME CNE6SSM<Semianalytic>

CNE6SSM_susy_scale_constraint<Semianalytic>::CNE6SSM_susy_scale_constraint()
   : Constraint<Semianalytic>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
{
}

CNE6SSM_susy_scale_constraint<Semianalytic>::CNE6SSM_susy_scale_constraint(
   CNE6SSM<Semianalytic>* model_)
   : Constraint<Semianalytic>()
   , model(model_)
{
   initialize();
}

CNE6SSM_susy_scale_constraint<Semianalytic>::~CNE6SSM_susy_scale_constraint()
{
}

void CNE6SSM_susy_scale_constraint<Semianalytic>::apply()
{
   assert(model && "Error: CNE6SSM_susy_scale_constraint::apply():"
          " model pointer must not be zero");

   model->calculate_coefficients(high_constraint->get_scale());
   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints

   // the parameters, which are fixed by the EWSB eqs., will now be
   // defined at this scale (at the EWSB loop level defined in the
   // model)
   model->solve_ewsb();
}

double CNE6SSM_susy_scale_constraint<Semianalytic>::get_scale() const
{
   return scale;
}

double CNE6SSM_susy_scale_constraint<Semianalytic>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const CNE6SSM_input_parameters<Semianalytic>& CNE6SSM_susy_scale_constraint<Semianalytic>::get_input_parameters() const
{
   assert(model && "Error: CNE6SSM_susy_scale_constraint::"
          "get_input_parameters(): model pointer is zero.");

   return model->get_input();
}

CNE6SSM<Semianalytic>* CNE6SSM_susy_scale_constraint<Semianalytic>::get_model() const
{
   return model;
}

void CNE6SSM_susy_scale_constraint<Semianalytic>::set_model(Semianalytic_model* model_)
{
   model = cast_model<CNE6SSM<Semianalytic>*>(model_);
}

void CNE6SSM_susy_scale_constraint<Semianalytic>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void CNE6SSM_susy_scale_constraint<Semianalytic>::initialize()
{
   assert(model && "CNE6SSM_susy_scale_constraint<Semianalytic>::"
          "initialize(): model pointer is zero.");

   const auto sInput = INPUTPARAMETER(sInput);
   const auto SigmaxInput = INPUTPARAMETER(SigmaxInput);

   initial_scale_guess = SigmaxInput * sInput;

   scale = initial_scale_guess;
}

void CNE6SSM_susy_scale_constraint<Semianalytic>::update_scale()
{
   assert(model && "CNE6SSM_susy_scale_constraint<Semianalytic>::"
          "update_scale(): model pointer is zero.");

   const auto ZU = MODELPARAMETER(ZU);
   const auto MSu = MODELPARAMETER(MSu);

   scale = Sqrt(Power(MSu(0),Sqr(Abs(ZU(0,2))) + Sqr(Abs(ZU(0,5))))*Power(MSu(1
      ),Sqr(Abs(ZU(1,2))) + Sqr(Abs(ZU(1,5))))*Power(MSu(2),Sqr(Abs(ZU(2,2))) +
      Sqr(Abs(ZU(2,5))))*Power(MSu(3),Sqr(Abs(ZU(3,2))) + Sqr(Abs(ZU(3,5))))*Power
      (MSu(4),Sqr(Abs(ZU(4,2))) + Sqr(Abs(ZU(4,5))))*Power(MSu(5),Sqr(Abs(ZU(5,2))
      ) + Sqr(Abs(ZU(5,5)))));


}

} // namespace flexiblesusy
