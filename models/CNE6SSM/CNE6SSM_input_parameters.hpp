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

// File generated at Mon 15 Sep 2014 17:29:07

#ifndef CNE6SSM_INPUT_PARAMETERS_H
#define CNE6SSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct CNE6SSM_input_parameters {
   double m0;
   double m12;
   double TanBeta;
   double Azero;
   double ssumInput;
   double QS;
   Eigen::Matrix<double,3,2> hEInput;
   double SigmaLInput;
   double KappaPrInput;
   double SigmaxInput;
   Eigen::Matrix<double,3,3> gDInput;
   Eigen::Matrix<double,3,3> KappaInput;
   Eigen::Matrix<double,2,2> Lambda12Input;
   double LambdaxInput;
   Eigen::Matrix<double,3,2> fuInput;
   Eigen::Matrix<double,3,2> fdInput;
   double MuPrInput;
   double MuPhiInput;
   double BMuPrInput;
   double BMuPhiInput;

   CNE6SSM_input_parameters()
      : m0(0), m12(0), TanBeta(0), Azero(0), ssumInput(0), QS(0), hEInput(
   Eigen::Matrix<double,3,2>::Zero()), SigmaLInput(0), KappaPrInput(0),
   SigmaxInput(0), gDInput(Eigen::Matrix<double,3,3>::Zero()), KappaInput(
   Eigen::Matrix<double,3,3>::Zero()), Lambda12Input(Eigen::Matrix<double,2,2>
   ::Zero()), LambdaxInput(0), fuInput(Eigen::Matrix<double,3,2>::Zero()),
   fdInput(Eigen::Matrix<double,3,2>::Zero()), MuPrInput(0), MuPhiInput(0),
   BMuPrInput(0), BMuPhiInput(0)

   {}
};

} // namespace flexiblesusy

#endif
