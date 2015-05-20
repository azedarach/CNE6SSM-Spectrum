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

// File generated at Sun 19 Apr 2015 20:31:38

#include "CNE6SSMSusy_info.hpp"

#include <iostream>

namespace flexiblesusy {

namespace CNE6SSMSusy_info {
   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;
   const double normalization_g1p = 0.15811388300841897;

   const unsigned particle_multiplicities[NUMBER_OF_PARTICLES] = {1, 1, 3, 1, 1
      , 1, 1, 6, 3, 6, 6, 6, 5, 5, 2, 8, 2, 3, 3, 3, 3, 7, 4, 2, 7, 2, 2, 2, 1};

   const char* particle_names[NUMBER_OF_PARTICLES] = {"VG", "Glu", "Fv", "ChaP"
      , "VP", "VZ", "VZp", "Sd", "Sv", "Su", "Se", "SDX", "hh", "Ah", "Hpm", "Chi"
      , "Cha", "Fe", "Fd", "Fu", "FDX", "SHI0", "SHIPM", "ChaI", "ChiI", "SHp0",
      "SHpp", "ChiP", "VWm"};

   const char* particle_latex_names[NUMBER_OF_PARTICLES] = {   "g",
      "\\tilde{g}", "\\nu", "\\tilde{\\chi'}^{-}", "\\gamma", "Z", "{Z\\ {}^\\prime}",
      "\\tilde{d}", "\\tilde{\\nu}", "\\tilde{u}", "\\tilde{e}", "\\tilde{D}", "h"
      , "A^0", "H^{\\pm}", "\\tilde{\\chi}^0", "\\tilde{\\chi}^-", "e", "d", "u", "D",
      "h_I^{0}", "h_I^{-}", "\\tilde{\\chi}_I^{-}",
      "\\tilde{\\chi}_I^{0}", "L_4^{0}", "L_4^{-}", "\\tilde{\\chi'}^{0}", "W^-"
      };

   const char* parameter_names[NUMBER_OF_PARAMETERS] = {"Yd(0,0)", "Yd(0,1)",
      "Yd(0,2)", "Yd(1,0)", "Yd(1,1)", "Yd(1,2)", "Yd(2,0)", "Yd(2,1)", "Yd(2,2)",
      "hE(0,0)", "hE(0,1)", "hE(1,0)", "hE(1,1)", "hE(2,0)", "hE(2,1)", "Ye(0,0)"
      , "Ye(0,1)", "Ye(0,2)", "Ye(1,0)", "Ye(1,1)", "Ye(1,2)", "Ye(2,0)",
      "Ye(2,1)", "Ye(2,2)", "SigmaL", "KappaPr", "Sigmax", "gD(0,0)", "gD(0,1)",
      "gD(0,2)", "gD(1,0)", "gD(1,1)", "gD(1,2)", "gD(2,0)", "gD(2,1)", "gD(2,2)",
      "Kappa(0,0)", "Kappa(0,1)", "Kappa(0,2)", "Kappa(1,0)", "Kappa(1,1)",
      "Kappa(1,2)", "Kappa(2,0)", "Kappa(2,1)", "Kappa(2,2)", "Lambda12(0,0)",
      "Lambda12(0,1)", "Lambda12(1,0)", "Lambda12(1,1)", "Lambdax", "fu(0,0)",
      "fu(0,1)", "fu(1,0)", "fu(1,1)", "fu(2,0)", "fu(2,1)", "fd(0,0)", "fd(0,1)",
      "fd(1,0)", "fd(1,1)", "fd(2,0)", "fd(2,1)", "Yu(0,0)", "Yu(0,1)", "Yu(0,2)"
      , "Yu(1,0)", "Yu(1,1)", "Yu(1,2)", "Yu(2,0)", "Yu(2,1)", "Yu(2,2)", "MuPr",
      "MuPhi", "XiF", "g1", "g2", "g3", "g1p", "vd", "vu", "vs", "vsb", "vphi",
      "QS"};

   const int parameter_mass_dimensions[NUMBER_OF_PARAMETERS] = {0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 1, 1, 2, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0};

   const char* parameter_latex_names[NUMBER_OF_PARAMETERS] = {"y_{11}^D",
      "y_{12}^D", "y_{13}^D", "y_{21}^D", "y_{22}^D", "y_{23}^D", "y_{31}^D",
      "y_{32}^D", "y_{33}^D", "h_{11}^E", "h_{12}^E", "h_{21}^E", "h_{22}^E",
      "h_{31}^E", "h_{32}^E", "y_{11}^E", "y_{12}^E", "y_{13}^E", "y_{21}^E",
      "y_{22}^E", "y_{23}^E", "y_{31}^E", "y_{32}^E", "y_{33}^E", "\\sigma_L",
      "\\kappa", "\\sigma", "g_{11}^D", "g_{12}^D", "g_{13}^D", "g_{21}^D",
      "g_{22}^D", "g_{23}^D", "g_{31}^D", "g_{32}^D", "g_{33}^D", "\\kappa_{11}",
      "\\kappa_{12}", "\\kappa_{13}", "\\kappa_{21}", "\\kappa_{22}", 
      "\\kappa_{23}", "\\kappa_{31}", "\\kappa_{32}", "\\kappa_{33}", "\\lambda_{11}",
      "\\lambda_{12}", "\\lambda_{21}", "\\lambda_{22}", "\\lambda", "\\tilde{f}_{11}",
      "\\tilde{f}_{12}", "\\tilde{f}_{21}", "\\tilde{f}_{22}", "\\tilde{f}_{31}",
      "\\tilde{f}_{32}", "f_{11}", "f_{12}", "f_{21}", "f_{22}", "f_{31}", "f_{32}",
      "y_{11}^U", "y_{12}^U", "y_{13}^U", "y_{21}^U", "y_{22}^U", "y_{23}^U",
      "y_{31}^U", "y_{32}^U", "y_{33}^U", "\\mu_L", "\\mu", "\\Lambda", "g_1", 
      "g_2", "g_3", "g_1'", "v_1", "v_2", "s_1", "s_2", "\\varphi", "Q_S"};

   const char* input_names[NUMBER_OF_INPUTS] = {"TanBeta", "sInput", "TanTheta",
      "vphiInput", "QSInput", "hEInput(0,0)", "hEInput(0,1)", "hEInput(1,0)", "hEInput(1,1)", 
      "hEInput(2,0)", "hEInput(2,1)", "SigmaLInput", "KappaPrInput", "SigmaxInput", 
      "gDInput(0,0)", "gDInput(0,1)", "gDInput(0,2)", "gDInput(1,0)", "gDInput(1,1)", 
      "gDInput(1,2)", "gDInput(2,0)", "gDInput(2,1)", "gDInput(2,2)", "KappaInput(0,0)",
      "KappaInput(0,1)", "KappaInput(0,2)", "KappaInput(1,0)", "KappaInput(1,1)", 
      "KappaInput(1,2)", "KappaInput(2,0)", "KappaInput(2,1)", "KappaInput(2,2)", 
      "Lambda12Input(0,0)", "Lambda12Input(0,1)", "Lambda12Input(1,0)", "Lambda12Input(1,1)",
      "LambdaxInput", "fuInput(0,0)", "fuInput(0,1)", "fuInput(1,0)", "fuInput(1,1)", "fuInput(2,0)",
      "fuInput(2,1)", "fdInput(0,0)", "fdInput(0,1)", "fdInput(1,0)", "fdInput(1,1)",
      "fdInput(2,0)", "fdInput(2,1)", "MuPrInput", "MuPhiInput", "XiFInput"};

   const int input_mass_dimensions[NUMBER_OF_INPUTS] = {0, 1, 0, 1, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2};

   const char* input_latex_names[NUMBER_OF_INPUTS] = {"\\tan\\beta", "s", "\\tan\\theta",
      "\\varphi", "\\tilde{Q}_S", "h_{11}^E", "h_{12}^E",
      "h_{21}^E", "h_{22}^E", "h_{31}^E", "h_{32}^E", "\\sigma_L", "\\kappa", "\\sigma",
      "g_{11}^D", "g_{12}^D", "g_{13}^D", "g_{21}^D", "g_{22}^D", "g_{23}^D", "g_{31}^D",
      "g_{32}^D", "g_{33}^D", "\\kappa_{11}", "\\kappa_{12}", "\\kappa_{13}", 
      "\\kappa_{21}", "\\kappa_{22}", "\\kappa_{23}", "\\kappa_{31}", "\\kappa_{32}",
      "\\kappa_{33}", "\\lambda_{11}", "\\lambda_{12}", "\\lambda_{21}", "\\lambda_{22}",
      "\\lambda", "\\tilde{f}_{11}", "\\tilde{f}_{12}", "\\tilde{f}_{21}", "\\tilde{f}_{22}", 
      "\\tilde{f}_{31}", "\\tilde{f}_{32}", "f_{11}", "f_{12}", "f_{21}", "f_{22}", 
      "f_{31}", "f_{32}", "\\mu_L", "\\mu", "\\Lambda"};

   const char* model_name = "CNE6SSMSusy";
   const bool is_low_energy_model = false;
   const bool is_supersymmetric_model = true;

void print(std::ostream& ostr)
{
   ostr
      << "Model information\n"
      << "=================\n"
      << "Model name:                " << model_name << '\n'
      << "Is a low-energy model:     "
      << (is_low_energy_model ? "yes" : "no") << '\n'
      << "Is a supersymmetric model: "
      << (is_supersymmetric_model ? "yes" : "no") << '\n'
      << "Number of multiplets:      " << NUMBER_OF_PARTICLES << '\n'
      << "Number of parameters:      " << NUMBER_OF_PARAMETERS << '\n'
      ;

   ostr << "\n"
      "Multiplets:                ";
   for (unsigned i = 0; i < NUMBER_OF_PARTICLES; i++) {
      ostr << particle_names[i]
           << '[' << particle_multiplicities[i] << ']';
      if (i + 1 < NUMBER_OF_PARTICLES)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Parameters:                ";
   for (unsigned i = 0; i < NUMBER_OF_PARAMETERS; i++) {
      ostr << parameter_names[i];
      if (i + 1 < NUMBER_OF_PARAMETERS)
         ostr << ", ";
   }
   ostr << '\n';
}

} // namespace CNE6SSM_info

} // namespace flexiblesusy

