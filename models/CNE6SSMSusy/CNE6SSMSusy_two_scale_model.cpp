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

// File generated at Sun 19 Apr 2015 20:37:13

/**
 * @file CNE6SSMSusy_two_scale_model.cpp
 * @brief implementation of the CNE6SSMSusy two scale model class
 *
 * This file was generated at Sun 19 Apr 2015 20:37:13 with FlexibleSUSY
 * 1.0.4 (git commit: v1.0.3-cpc6-749-g227a308) and SARAH 4.5.3 .
 */

#include "CNE6SSMSusy_two_scale_model.hpp"
#include "numerics.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "root_finder.hpp"
#include "fixed_point_iterator.hpp"
#include "gsl_utils.hpp"
#include "config.h"
#include "functors.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

#include <gsl/gsl_multiroots.h>

namespace flexiblesusy {

using namespace CNE6SSMSusy_info;

#define CLASSNAME CNE6SSMSusy<Two_scale>

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model->get_input().parameter
#define LOCALINPUT(parameter) input.parameter

CLASSNAME::CNE6SSMSusy(const CNE6SSMSusy_input_parameters<Two_scale>& input_)
   : Two_scale_model()
   , CNE6SSMSusy_susy_parameters()
   , input(input_)
   , precision(1.0e-3)
   , problems(CNE6SSMSusy_info::particle_names)
{
}

CLASSNAME::~CNE6SSMSusy()
{
}

void CLASSNAME::set_precision(double precision_)
{
   precision = precision_;
}

const CNE6SSMSusy_input_parameters<Two_scale>& CLASSNAME::get_input() const
{
   return input;
}

void CLASSNAME::set_input_parameters(const CNE6SSMSusy_input_parameters<Two_scale>& input_)
{
   input = input_;
}

const Problems<CNE6SSMSusy_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems() const
{
   return problems;
}

Problems<CNE6SSMSusy_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems()
{
   return problems;
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "CNE6SSMSusy (solver type: two_scale)\n"
           "========================================\n";
   CNE6SSMSusy_susy_parameters::print(ostr);
}

/**
 * calculates spectrum for model once the DRbar parameters at
 * at low energies are known
 */
void CLASSNAME::calculate_spectrum()
{

}

void CLASSNAME::clear_problems()
{

}

void CLASSNAME::clear()
{
   CNE6SSMSusy_susy_parameters::clear();
}

std::string CLASSNAME::name() const
{
   return "CNE6SSMSusy";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   CNE6SSMSusy_susy_parameters::run_to(scale, eps);
}

double CLASSNAME::get_parameter(unsigned parameter) const
{
   if (parameter >= CNE6SSMSusy_info::NUMBER_OF_PARAMETERS)
      throw UnknownModelParameterError(parameter);

   switch (parameter) {

   case CNE6SSMSusy_info::Yd00:
      return Yd(0,0);
   case CNE6SSMSusy_info::Yd01:
      return Yd(0,1);
   case CNE6SSMSusy_info::Yd02:
      return Yd(0,2);
   case CNE6SSMSusy_info::Yd10:
      return Yd(1,0); 
   case CNE6SSMSusy_info::Yd11:
      return Yd(1,1); 
   case CNE6SSMSusy_info::Yd12:
      return Yd(1,2);
   case CNE6SSMSusy_info::Yd20:
      return Yd(2,0);
   case CNE6SSMSusy_info::Yd21:
      return Yd(2,1);
   case CNE6SSMSusy_info::Yd22:
      return Yd(2,2);
   case CNE6SSMSusy_info::hE00:
      return hE(0,0);
   case CNE6SSMSusy_info::hE01:
      return hE(0,1);
   case CNE6SSMSusy_info::hE10:
      return hE(1,0); 
   case CNE6SSMSusy_info::hE11:
      return hE(1,1); 
   case CNE6SSMSusy_info::hE20:
      return hE(2,0); 
   case CNE6SSMSusy_info::hE21:
      return hE(2,1);
   case CNE6SSMSusy_info::Ye00:
      return Ye(0,0); 
   case CNE6SSMSusy_info::Ye01:
      return Ye(0,1);
   case CNE6SSMSusy_info::Ye02:
      return Ye(0,2); 
   case CNE6SSMSusy_info::Ye10:
      return Ye(1,0); 
   case CNE6SSMSusy_info::Ye11:
      return Ye(1,1); 
   case CNE6SSMSusy_info::Ye12:
      return Ye(1,2);
   case CNE6SSMSusy_info::Ye20:
      return Ye(2,0); 
   case CNE6SSMSusy_info::Ye21:
      return Ye(2,1); 
   case CNE6SSMSusy_info::Ye22:
      return Ye(2,2); 
   case CNE6SSMSusy_info::SigmaL:
      return SigmaL; 
   case CNE6SSMSusy_info::KappaPr:
      return KappaPr; 
   case CNE6SSMSusy_info::Sigmax:
      return Sigmax; 
   case CNE6SSMSusy_info::gD00:
      return gD(0,0); 
   case CNE6SSMSusy_info::gD01:
      return gD(0,1); 
   case CNE6SSMSusy_info::gD02:
      return gD(0,2); 
   case CNE6SSMSusy_info::gD10:
      return gD(1,0); 
   case CNE6SSMSusy_info::gD11:
      return gD(1,1);
   case CNE6SSMSusy_info::gD12:
      return gD(1,2); 
   case CNE6SSMSusy_info::gD20:
      return gD(2,0);
   case CNE6SSMSusy_info::gD21:
      return gD(2,1); 
   case CNE6SSMSusy_info::gD22:
      return gD(2,2); 
   case CNE6SSMSusy_info::Kappa00:
      return Kappa(0,0); 
   case CNE6SSMSusy_info::Kappa01:
      return Kappa(0,1); 
   case CNE6SSMSusy_info::Kappa02:
      return Kappa(0,2); 
   case CNE6SSMSusy_info::Kappa10:
      return Kappa(1,0); 
   case CNE6SSMSusy_info::Kappa11:
      return Kappa(1,1); 
   case CNE6SSMSusy_info::Kappa12:
      return Kappa(1,2);
   case CNE6SSMSusy_info::Kappa20:
      return Kappa(2,0); 
   case CNE6SSMSusy_info::Kappa21:
      return Kappa(2,1); 
   case CNE6SSMSusy_info::Kappa22:
      return Kappa(2,2); 
   case CNE6SSMSusy_info::Lambda1200:
      return Lambda12(0,0); 
   case CNE6SSMSusy_info::Lambda1201:
      return Lambda12(0,1); 
   case CNE6SSMSusy_info::Lambda1210:
      return Lambda12(1,0); 
   case CNE6SSMSusy_info::Lambda1211:
      return Lambda12(1,1);
   case CNE6SSMSusy_info::Lambdax:
      return Lambdax;
   case CNE6SSMSusy_info::fu00:
      return fu(0,0);
   case CNE6SSMSusy_info::fu01:
      return fu(0,1); 
   case CNE6SSMSusy_info::fu10:
      return fu(1,0);
   case CNE6SSMSusy_info::fu11:
      return fu(1,1); 
   case CNE6SSMSusy_info::fu20:
      return fu(2,0); 
   case CNE6SSMSusy_info::fu21:
      return fu(2,1); 
   case CNE6SSMSusy_info::fd00:
      return fd(0,0); 
   case CNE6SSMSusy_info::fd01:
      return fd(0,1);
   case CNE6SSMSusy_info::fd10:
      return fd(1,0); 
   case CNE6SSMSusy_info::fd11:
      return fd(1,1); 
   case CNE6SSMSusy_info::fd20:
      return fd(2,0);
   case CNE6SSMSusy_info::fd21:
      return fd(2,1); 
   case CNE6SSMSusy_info::Yu00:
      return Yu(0,0); 
   case CNE6SSMSusy_info::Yu01:
      return Yu(0,1); 
   case CNE6SSMSusy_info::Yu02:
      return Yu(0,2); 
   case CNE6SSMSusy_info::Yu10:
      return Yu(1,0); 
   case CNE6SSMSusy_info::Yu11:
      return Yu(1,1); 
   case CNE6SSMSusy_info::Yu12:
      return Yu(1,2); 
   case CNE6SSMSusy_info::Yu20:
      return Yu(2,0); 
   case CNE6SSMSusy_info::Yu21:
      return Yu(2,1); 
   case CNE6SSMSusy_info::Yu22:
      return Yu(2,2); 
   case CNE6SSMSusy_info::MuPr:
      return MuPr;
   case CNE6SSMSusy_info::MuPhi:
      return MuPhi; 
   case CNE6SSMSusy_info::XiF:
      return XiF;
   case CNE6SSMSusy_info::g1:
      return g1; 
   case CNE6SSMSusy_info::g2:
      return g2;
   case CNE6SSMSusy_info::g3:
      return g3; 
   case CNE6SSMSusy_info::g1p:
      return g1p; 
   case CNE6SSMSusy_info::vd:
      return vd; 
   case CNE6SSMSusy_info::vu:
      return vu; 
   case CNE6SSMSusy_info::vs:
      return vs; 
   case CNE6SSMSusy_info::vsb:
      return vsb; 
   case CNE6SSMSusy_info::vphi:
      return vphi; 
   case CNE6SSMSusy_info::QS:
      return QS;

   default:
      throw UnknownModelParameterError(parameter);
   }
}

void CLASSNAME::set_parameter(unsigned parameter, double x)
{
   if (parameter >= CNE6SSMSusy_info::NUMBER_OF_PARAMETERS)
      throw UnknownModelParameterError(parameter);

   switch (parameter) {

   case CNE6SSMSusy_info::Yd00:
      Yd(0,0) = x;
      break;
   case CNE6SSMSusy_info::Yd01:
      Yd(0,1) = x;
      break;
   case CNE6SSMSusy_info::Yd02:
      Yd(0,2) = x;
      break;
   case CNE6SSMSusy_info::Yd10:
      Yd(1,0) = x;
      break;
   case CNE6SSMSusy_info::Yd11:
      Yd(1,1) = x;
      break;
   case CNE6SSMSusy_info::Yd12:
      Yd(1,2) = x;
      break;
   case CNE6SSMSusy_info::Yd20:
      Yd(2,0) = x;
      break;
   case CNE6SSMSusy_info::Yd21:
      Yd(2,1) = x;
      break;
   case CNE6SSMSusy_info::Yd22:
      Yd(2,2) = x;
      break;
   case CNE6SSMSusy_info::hE00:
      hE(0,0) = x;
      break;
   case CNE6SSMSusy_info::hE01:
      hE(0,1) = x;
      break;
   case CNE6SSMSusy_info::hE10:
      hE(1,0) = x;
      break;
   case CNE6SSMSusy_info::hE11:
      hE(1,1) = x;
      break;
   case CNE6SSMSusy_info::hE20:
      hE(2,0) = x;
      break;
   case CNE6SSMSusy_info::hE21:
      hE(2,1) = x;
      break;
   case CNE6SSMSusy_info::Ye00:
      Ye(0,0) = x; 
      break;
   case CNE6SSMSusy_info::Ye01:
      Ye(0,1) = x;
      break;
   case CNE6SSMSusy_info::Ye02:
      Ye(0,2) = x;
      break;
   case CNE6SSMSusy_info::Ye10:
      Ye(1,0) = x;
      break;
   case CNE6SSMSusy_info::Ye11:
      Ye(1,1) = x;
      break;
   case CNE6SSMSusy_info::Ye12:
      Ye(1,2) = x;
      break;
   case CNE6SSMSusy_info::Ye20:
      Ye(2,0) = x;
      break;
   case CNE6SSMSusy_info::Ye21:
      Ye(2,1) = x;
      break;
   case CNE6SSMSusy_info::Ye22:
      Ye(2,2) = x;
      break;
   case CNE6SSMSusy_info::SigmaL:
      SigmaL = x;
      break;
   case CNE6SSMSusy_info::KappaPr:
      KappaPr = x;
      break;
   case CNE6SSMSusy_info::Sigmax:
      Sigmax = x;
      break;
   case CNE6SSMSusy_info::gD00:
      gD(0,0) = x;
      break;
   case CNE6SSMSusy_info::gD01:
      gD(0,1) = x;
      break;
   case CNE6SSMSusy_info::gD02:
      gD(0,2) = x;
      break;
   case CNE6SSMSusy_info::gD10:
      gD(1,0) = x;
      break;
   case CNE6SSMSusy_info::gD11:
      gD(1,1) = x;
      break;
   case CNE6SSMSusy_info::gD12:
      gD(1,2) = x;
      break;
   case CNE6SSMSusy_info::gD20:
      gD(2,0) = x;
      break;
   case CNE6SSMSusy_info::gD21:
      gD(2,1) = x;
      break;
   case CNE6SSMSusy_info::gD22:
      gD(2,2) = x;
      break;
   case CNE6SSMSusy_info::Kappa00:
      Kappa(0,0) = x;
      break;
   case CNE6SSMSusy_info::Kappa01:
      Kappa(0,1) = x;
      break;
   case CNE6SSMSusy_info::Kappa02:
      Kappa(0,2) = x;
      break;
   case CNE6SSMSusy_info::Kappa10:
      Kappa(1,0) = x;
      break;
   case CNE6SSMSusy_info::Kappa11:
      Kappa(1,1) = x;
      break;
   case CNE6SSMSusy_info::Kappa12:
      Kappa(1,2) = x;
      break;
   case CNE6SSMSusy_info::Kappa20:
      Kappa(2,0) = x;
      break;
   case CNE6SSMSusy_info::Kappa21:
      Kappa(2,1) = x;
      break;
   case CNE6SSMSusy_info::Kappa22:
      Kappa(2,2) = x;
      break;
   case CNE6SSMSusy_info::Lambda1200:
      Lambda12(0,0) = x;
      break;
   case CNE6SSMSusy_info::Lambda1201:
      Lambda12(0,1) = x;
      break;
   case CNE6SSMSusy_info::Lambda1210:
      Lambda12(1,0) = x;
      break;
   case CNE6SSMSusy_info::Lambda1211:
      Lambda12(1,1) = x;
      break;
   case CNE6SSMSusy_info::Lambdax:
      Lambdax = x;
      break;
   case CNE6SSMSusy_info::fu00:
      fu(0,0) = x;
      break;
   case CNE6SSMSusy_info::fu01:
      fu(0,1) = x;
      break;
   case CNE6SSMSusy_info::fu10:
      fu(1,0) = x;
      break;
   case CNE6SSMSusy_info::fu11:
      fu(1,1) = x;
      break;
   case CNE6SSMSusy_info::fu20:
      fu(2,0) = x;
      break;
   case CNE6SSMSusy_info::fu21:
      fu(2,1) = x;
      break;
   case CNE6SSMSusy_info::fd00:
      fd(0,0) = x;
      break;
   case CNE6SSMSusy_info::fd01:
      fd(0,1) = x;
      break;
   case CNE6SSMSusy_info::fd10:
      fd(1,0) = x;
      break;
   case CNE6SSMSusy_info::fd11:
      fd(1,1) = x;
      break;
   case CNE6SSMSusy_info::fd20:
      fd(2,0) = x;
      break;
   case CNE6SSMSusy_info::fd21:
      fd(2,1) = x;
      break;
   case CNE6SSMSusy_info::Yu00:
      Yu(0,0) = x;
      break;
   case CNE6SSMSusy_info::Yu01:
      Yu(0,1) = x;
      break;
   case CNE6SSMSusy_info::Yu02:
      Yu(0,2) = x;
      break;
   case CNE6SSMSusy_info::Yu10:
      Yu(1,0) = x;
      break;
   case CNE6SSMSusy_info::Yu11:
      Yu(1,1) = x;
      break;
   case CNE6SSMSusy_info::Yu12:
      Yu(1,2) = x;
      break;
   case CNE6SSMSusy_info::Yu20:
      Yu(2,0) = x;
      break;
   case CNE6SSMSusy_info::Yu21:
      Yu(2,1) = x;
      break;
   case CNE6SSMSusy_info::Yu22:
      Yu(2,2) = x;
      break;
   case CNE6SSMSusy_info::MuPr:
      MuPr = x;
      break;
   case CNE6SSMSusy_info::MuPhi:
      MuPhi = x;
      break;
   case CNE6SSMSusy_info::XiF:
      XiF = x;
      break;
   case CNE6SSMSusy_info::g1:
      g1 = x;
      break;
   case CNE6SSMSusy_info::g2:
      g2 = x;
      break;
   case CNE6SSMSusy_info::g3:
      g3 = x;
      break;
   case CNE6SSMSusy_info::g1p:
      g1p = x;
      break;
   case CNE6SSMSusy_info::vd:
      vd = x;
      break;
   case CNE6SSMSusy_info::vu:
      vu = x;
      break;
   case CNE6SSMSusy_info::vs:
      vs = x;
      break;
   case CNE6SSMSusy_info::vsb:
      vsb = x;
      break;
   case CNE6SSMSusy_info::vphi:
      vphi = x;
      break;
   case CNE6SSMSusy_info::QS:
      QS = x;
      break;

   default:
      throw UnknownModelParameterError(parameter);
   }
}

std::ostream& operator<<(std::ostream& ostr, const CNE6SSMSusy<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
