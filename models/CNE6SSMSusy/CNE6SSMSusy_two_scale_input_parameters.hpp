// ====================================================================
// Specialisation of input parameters for models solved using the
// two scale algorithm
// ====================================================================

#ifndef CNE6SSMSusy_TWO_SCALE_INPUT_PARAMETERS_H
#define CNE6SSMSusy_TWO_SCALE_INPUT_PARAMETERS_H

#include "CNE6SSMSusy_input_parameters.hpp"

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

class Two_scale;

template<>
struct CNE6SSMSusy_input_parameters<Two_scale> {
   double TanBeta;
   double sInput;
   double TanTheta;
   double vphiInput;
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
   double XiFInput;

   CNE6SSMSusy_input_parameters()
      : TanBeta(0), sInput(0), TanTheta(0), vphiInput(0), QSInput(0),
   hEInput(Eigen::Matrix<double,3,2>::Zero()), SigmaLInput(0), KappaPrInput(0),
   SigmaxInput(0), gDInput(Eigen::Matrix<double,3,3>::Zero()), KappaInput(
   Eigen::Matrix<double,3,3>::Zero()), Lambda12Input(Eigen::Matrix<double,2,2>
   ::Zero()), LambdaxInput(0), fuInput(Eigen::Matrix<double,3,2>::Zero()),
        fdInput(Eigen::Matrix<double,3,2>::Zero()), MuPrInput(0), MuPhiInput(0), XiFInput(0)

   {}
};

std::ostream& operator<<(std::ostream&, const CNE6SSMSusy_input_parameters<Two_scale>&);

} // namespace flexiblesusy

#endif
