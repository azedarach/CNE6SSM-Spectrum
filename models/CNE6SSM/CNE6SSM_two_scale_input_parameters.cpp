// ====================================================================
// Specialisation of input parameters for models solved using the
// two scale algorithm
// ====================================================================

#include "CNE6SSM_two_scale_input_parameters.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

std::ostream& operator<<(std::ostream& ostr, const CNE6SSM_input_parameters<Two_scale>& input)
{
   ostr << "m0 = " << INPUT(m0) << ", ";
   ostr << "m12 = " << INPUT(m12) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "SignLambdax = " << INPUT(SignLambdax) << ", ";
   ostr << "Azero = " << INPUT(Azero) << ", ";
   ostr << "ssumInput = " << INPUT(ssumInput) << ", ";
   ostr << "QS = " << INPUT(QSInput) << ", ";
   ostr << "hEInput = " << INPUT(hEInput) << ", ";
   ostr << "SigmaLInput = " << INPUT(SigmaLInput) << ", ";
   ostr << "KappaPrInput = " << INPUT(KappaPrInput) << ", ";
   ostr << "SigmaxInput = " << INPUT(SigmaxInput) << ", ";
   ostr << "gDInput = " << INPUT(gDInput) << ", ";
   ostr << "KappaInput = " << INPUT(KappaInput) << ", ";
   ostr << "Lambda12Input = " << INPUT(Lambda12Input) << ", ";
   ostr << "fuInput = " << INPUT(fuInput) << ", ";
   ostr << "fdInput = " << INPUT(fdInput) << ", ";
   ostr << "MuPrInput = " << INPUT(MuPrInput) << ", ";
   ostr << "MuPhiInput = " << INPUT(MuPhiInput) << ", ";
   ostr << "BMuPrInput = " << INPUT(BMuPrInput) << ", ";
   ostr << "BMuPhiInput = " << INPUT(BMuPhiInput) << ", ";

   return ostr;
}

CNE6SSM_input_parameters_slha CNE6SSM_input_parameters<Two_scale>::get_slha_format() const
{
   CNE6SSM_input_parameters_slha slha_format;

   slha_format.minpar.push_back(SLHA_input_entry("MINPAR", "m0", 1));
   slha_format.minpar.push_back(SLHA_input_entry("MINPAR", "m12", 2));
   slha_format.minpar.push_back(SLHA_input_entry("MINPAR", "TanBeta", 3));
   slha_format.minpar.push_back(SLHA_input_entry("MINPAR", "SignLambdax", 4));
   slha_format.minpar.push_back(SLHA_input_entry("MINPAR", "Azero", 5));

   slha_format.extpar.push_back(SLHA_input_entry("EXTPAR", "ssumInput", 65));
   slha_format.extpar.push_back(SLHA_input_entry("EXTPAR", "QSInput", 72));

   slha_format.entries.push_back(SLHA_input_entry("HMIXIN", "MuPhiInput", 31));
   slha_format.entries.push_back(SLHA_input_entry("HMIXIN", "KappaPrInput", 32));
   slha_format.entries.push_back(SLHA_input_entry("HMIXING", "SigmaxInput", 33));
   slha_format.blocks.push_back(SLHA_input_block("ESIXHEYUKIN", "hEInput"));
   slha_format.entries.push_back(SLHA_input_entry("ESIXRUNIN", "SigmaLInput", 42));
   slha_format.blocks.push_back(SLHA_input_block("ESIXGDYUKIN", "gDInput"));
   slha_format.blocks.push_back(SLHA_input_block("ESIXFUYUKIN", "fuInput"));
   slha_format.blocks.push_back(SLHA_input_block("ESIXFDYUKIN", "fdInput"));
   slha_format.entries.push_back(SLHA_input_entry("ESIXRUNIN", "BMuPhiInput", 30));
   slha_format.blocks.push_back(SLHA_input_block("ESIXKAPPAIN", "KappaInput"));
   slha_format.blocks.push_back(SLHA_input_block("ESIXLAMBDAIN", "Lambda12Input"));
   slha_format.entries.push_back(SLHA_input_entry("ESIXRUNIN", "MuPrInput", 0));
   slha_format.entries.push_back(SLHA_input_entry("ESIXRUNIN", "BMuPrInput", 101));

   return slha_format;
}

} // namespace flexiblesusy
