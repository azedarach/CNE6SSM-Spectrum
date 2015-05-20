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

#ifndef CNE6SSMSusy_INFO_H
#define CNE6SSMSusy_INFO_H

#include <iosfwd>

namespace flexiblesusy {

namespace CNE6SSMSusy_info {
   enum Particles : unsigned {VG, Glu, Fv, ChaP, VP, VZ, VZp, Sd, Sv, Su, Se,
      SDX, hh, Ah, Hpm, Chi, Cha, Fe, Fd, Fu, FDX, SHI0, SHIPM, ChaI, ChiI, SHp0,
      SHpp, ChiP, VWm, NUMBER_OF_PARTICLES};

   enum Parameters : unsigned {Yd00, Yd01, Yd02, Yd10, Yd11, Yd12, Yd20, Yd21,
      Yd22, hE00, hE01, hE10, hE11, hE20, hE21, Ye00, Ye01, Ye02, Ye10, Ye11, Ye12
      , Ye20, Ye21, Ye22, SigmaL, KappaPr, Sigmax, gD00, gD01, gD02, gD10, gD11,
      gD12, gD20, gD21, gD22, Kappa00, Kappa01, Kappa02, Kappa10, Kappa11, Kappa12
      , Kappa20, Kappa21, Kappa22, Lambda1200, Lambda1201, Lambda1210, Lambda1211,
      Lambdax, fu00, fu01, fu10, fu11, fu20, fu21, fd00, fd01, fd10, fd11, fd20,
      fd21, Yu00, Yu01, Yu02, Yu10, Yu11, Yu12, Yu20, Yu21, Yu22, MuPr, MuPhi, XiF
      , g1, g2, g3, g1p, vd, vu, vs, vsb, vphi, QS, NUMBER_OF_PARAMETERS};

   // DH:: added enum for input parameters as well
   enum Inputs : unsigned {TanBeta, sInput, TanTheta, vphiInput,
         QSInput, hEInput00, hEInput01, hEInput10, hEInput11, hEInput20,
         hEInput21, SigmaLInput, KappaPrInput, SigmaxInput, gDInput00,
         gDInput01, gDInput02, gDInput10, gDInput11, gDInput12, gDInput20,
         gDInput21, gDInput22, KappaInput00, KappaInput01, KappaInput02,
         KappaInput10, KappaInput11, KappaInput12, KappaInput20, KappaInput21,
         KappaInput22, Lambda12Input00, Lambda12Input01, Lambda12Input10,
         Lambda12Input11, LambdaxInput, fuInput00, fuInput01, fuInput10, fuInput11,
         fuInput20, fuInput21, fdInput00, fdInput01, fdInput10, fdInput11,
         fdInput20, fdInput21, MuPrInput, MuPhiInput, XiFInput,
         NUMBER_OF_INPUTS};

   extern const double normalization_g1;
   extern const double normalization_g2;
   extern const double normalization_g3;
   extern const double normalization_g1p;

   extern const unsigned particle_multiplicities[NUMBER_OF_PARTICLES];
   extern const char* particle_names[NUMBER_OF_PARTICLES];
   extern const char* particle_latex_names[NUMBER_OF_PARTICLES];
   extern const char* parameter_names[NUMBER_OF_PARAMETERS];
   extern const int parameter_mass_dimensions[NUMBER_OF_PARAMETERS];
   extern const char* parameter_latex_names[NUMBER_OF_PARAMETERS];
   extern const char* input_names[NUMBER_OF_INPUTS];
   extern const int input_mass_dimensions[NUMBER_OF_INPUTS];
   extern const char* input_latex_names[NUMBER_OF_INPUTS];

   extern const char* model_name;
   extern const bool is_low_energy_model;
   extern const bool is_supersymmetric_model;

   void print(std::ostream&);
}

} // namespace flexiblesusy

#endif
