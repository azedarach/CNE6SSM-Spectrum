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

// File generated at Sun 19 Apr 2015 20:31:35

#ifndef CNE6SSM_INPUT_PARAMETERS_H
#define CNE6SSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct CNE6SSM_input_parameters {
   double m0;
   double m12;
   double TanBeta;
   int SignLambdax;
   double Azero;
   double ssumInput;
   double QSInput;
   Eigen::Matrix<double,3,2> hEInput;
   double SigmaLInput;
   double KappaPrInput;
   double SigmaxInput;
   Eigen::Matrix<double,3,3> gDInput;
   Eigen::Matrix<double,3,3> KappaInput;
   Eigen::Matrix<double,2,2> Lambda12Input;
   Eigen::Matrix<double,3,2> fuInput;
   Eigen::Matrix<double,3,2> fdInput;
   double MuPrInput;
   double MuPhiInput;
   double BMuPrInput;
   double BMuPhiInput;

   CNE6SSM_input_parameters()
      : m0(0), m12(0), TanBeta(0), SignLambdax(1), Azero(0), ssumInput(0), QSInput(0),
   hEInput(Eigen::Matrix<double,3,2>::Zero()), SigmaLInput(0), KappaPrInput(0),
   SigmaxInput(0), gDInput(Eigen::Matrix<double,3,3>::Zero()), KappaInput(
   Eigen::Matrix<double,3,3>::Zero()), Lambda12Input(Eigen::Matrix<double,2,2>
   ::Zero()), fuInput(Eigen::Matrix<double,3,2>::Zero()), fdInput(Eigen::Matrix
   <double,3,2>::Zero()), MuPrInput(0), MuPhiInput(0), BMuPrInput(0),
   BMuPhiInput(0)

   {}
};

std::ostream& operator<<(std::ostream&, const CNE6SSM_input_parameters&);

} // namespace flexiblesusy

#endif
