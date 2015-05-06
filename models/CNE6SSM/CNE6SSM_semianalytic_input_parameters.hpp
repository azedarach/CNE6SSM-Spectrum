// ====================================================================
// Specialisation of input parameters for models solved using the
// semianalytic algorithm
// ====================================================================

#ifndef CNE6SSM_SEMIANALYTIC_INPUT_PARAMETERS_H
#define CNE6SSM_SEMIANALYTIC_INPUT_PARAMETERS_H

#include "CNE6SSM_input_parameters.hpp"

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

class Semianalytic;

template<>
struct CNE6SSM_input_parameters<Semianalytic> {
   double m12;
   double TanBeta;
   double sInput;
   double TanTheta;
   double QSInput;
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
      : m12(0), TanBeta(0), sInput(0), TanTheta(0), QSInput(0)
      , hEInput(Eigen::Matrix<double,3,2>::Zero()), SigmaLInput(0)
      , KappaPrInput(0), SigmaxInput(0)
      , gDInput(Eigen::Matrix<double,3,3>::Zero())
      , KappaInput(Eigen::Matrix<double,3,3>::Zero())
      , Lambda12Input(Eigen::Matrix<double,2,2>::Zero()), LambdaxInput(0)
      , fuInput(Eigen::Matrix<double,3,2>::Zero())
      , fdInput(Eigen::Matrix<double,3,2>::Zero()), MuPrInput(0)
      , MuPhiInput(0), BMuPrInput(0), BMuPhiInput(0)
      {}
};

std::ostream& operator<<(std::ostream&, const CNE6SSM_input_parameters<Semianalytic>&);

} // namespace flexiblesusy

#endif
