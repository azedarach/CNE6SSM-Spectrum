// ====================================================================
// Specialisation of input parameters for models solved using the
// two scale algorithm
// ====================================================================

#include "CNE6SSMSusy_two_scale_input_parameters.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

std::ostream& operator<<(std::ostream& ostr, const CNE6SSMSusy_input_parameters<Two_scale>& input)
{
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "sInput = " << INPUT(sInput) << ", ";
   ostr << "TanTheta = " << INPUT(TanTheta) << ", ";
   ostr << "vphiInput = " << INPUT(vphiInput) << ", ";
   ostr << "QSInput = " << INPUT(QSInput) << ", ";
   ostr << "hEInput = " << INPUT(hEInput) << ", ";
   ostr << "SigmaLInput = " << INPUT(SigmaLInput) << ", ";
   ostr << "KappaPrInput = " << INPUT(KappaPrInput) << ", ";
   ostr << "SigmaxInput = " << INPUT(SigmaxInput) << ", ";
   ostr << "gDInput = " << INPUT(gDInput) << ", ";
   ostr << "KappaInput = " << INPUT(KappaInput) << ", ";
   ostr << "Lambda12Input = " << INPUT(Lambda12Input) << ", ";
   ostr << "LambdaxInput = " << INPUT(LambdaxInput) << ", ";
   ostr << "fuInput = " << INPUT(fuInput) << ", ";
   ostr << "fdInput = " << INPUT(fdInput) << ", ";
   ostr << "MuPrInput = " << INPUT(MuPrInput) << ", ";
   ostr << "MuPhiInput = " << INPUT(MuPhiInput) << ", ";
   ostr << "XiFInput = " << INPUT(XiFInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy
