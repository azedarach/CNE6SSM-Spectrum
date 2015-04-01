// ====================================================================
// Reads in an SLHA file to define a CNE6SSM model at the GUT scale,
// then gets the coefficients of m0, m12, Azero etc. in mHd2, mHu2
// ====================================================================

#include "CNE6SSM_input_parameters.hpp"
#include "CNE6SSM_slha_io.hpp"
#include "CNE6SSM_spectrum_generator.hpp"

#include "error.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "command_line_options.hpp"

#include <gsl/gsl_deriv.h>

#include <iostream>
#include <chrono>
#include <cstdlib>
#include <map>
#include <vector>
#include <random>

using namespace flexiblesusy;
using namespace softsusy;

double get_parameter_from_inputs(CNE6SSM<Two_scale> model, CNE6SSM_info::Parameters p, 
                                 double m0, double m12, double Azero, double output_scale)
{
   const auto Ye = model.get_Ye();
   const auto Yd = model.get_Yd();
   const auto Yu = model.get_Yu();
   const auto KappaPr = model.get_KappaPr();
   const auto Sigmax = model.get_Sigmax();
   const auto hE = model.get_hE();
   const auto SigmaL = model.get_SigmaL();
   const auto gD = model.get_gD();
   const auto fu = model.get_fu();
   const auto fd = model.get_fd();
   const auto Kappa = model.get_Kappa();
   const auto Lambda12 = model.get_Lambda12();
   const auto Lambdax = model.get_Lambdax();

   model.set_TYe(Azero * Ye);
   model.set_TYd(Azero * Yd);
   model.set_TYu(Azero * Yu);
   model.set_TKappaPr(Azero * KappaPr);
   model.set_TSigmax(Azero * Sigmax);
   model.set_ThE(Azero * hE);
   model.set_TSigmaL(Azero * SigmaL);
   model.set_TgD(Azero * gD);
   model.set_Tfu(Azero * fu);
   model.set_Tfd(Azero * fd);
   model.set_TKappa(Azero * Kappa);
   model.set_TLambda12(Azero * Lambda12);
   model.set_TLambdax(Azero * Lambdax);
   model.set_mHd2(Sqr(m0));
   model.set_mHu2(Sqr(m0));
   model.set_ms2(Sqr(m0));
   model.set_msbar2(Sqr(m0));
   model.set_mphi2(Sqr(m0));
   model.set_mHp2(Sqr(m0));
   model.set_mHpbar2(Sqr(m0));
   model.set_mH1I2(Sqr(m0)*UNITMATRIX(2));
   model.set_mH2I2(Sqr(m0)*UNITMATRIX(2));
   model.set_mSI2(Sqr(m0)*UNITMATRIX(3));
   model.set_mq2(Sqr(m0)*UNITMATRIX(3));
   model.set_ml2(Sqr(m0)*UNITMATRIX(3));
   model.set_md2(Sqr(m0)*UNITMATRIX(3));
   model.set_mu2(Sqr(m0)*UNITMATRIX(3));
   model.set_me2(Sqr(m0)*UNITMATRIX(3));
   model.set_mDx2(Sqr(m0)*UNITMATRIX(3));
   model.set_mDxbar2(Sqr(m0)*UNITMATRIX(3));
   model.set_MassB(m12);
   model.set_MassWB(m12);
   model.set_MassG(m12);
   model.set_MassBp(m12);

   model.run_to(output_scale);

   return model.get_parameter(p);
}

double get_tree_level_Lambdax_soft_term(CNE6SSM<Two_scale> model, double mHd2_coeff, double mHu2_coeff)
{
   const auto vd = model.get_vd();
   const auto vu = model.get_vu();
   const auto vs = model.get_vs();

   double coeff = 2.0 * (Sqr(vd) * mHd2_coeff - Sqr(vu) * mHu2_coeff) / (Sqr(vs) * (Sqr(vu) - Sqr(vd)));

   return coeff;
}

double get_tree_level_Lambdax_constant_term(CNE6SSM<Two_scale> model)
{
   const auto g1 = model.get_g1();
   const auto g2 = model.get_g2();
   const auto g1p = model.get_g1p();
   const auto vd = model.get_vd();
   const auto vu = model.get_vu();
   const auto vs = model.get_vs();
   const auto vsb = model.get_vsb();
   const auto QS = model.get_input().QS;

   double coeff = - 0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) * (Sqr(vd) + Sqr(vu)) / Sqr(vs)
      + 0.025 * Sqr(g1p) * (2.0 * Sqr(vu) - 3.0 * Sqr(vd)) * 
      (-3.0 * Sqr(vd) - 2.0 * Sqr(vu) + QS * (Sqr(vs) - Sqr(vsb))) / (Sqr(vs) * (Sqr(vu) - Sqr(vd)));

   return coeff;
}

int main(int argc, const char * argv[])
{
   typedef Two_scale algorithm_type;

   Command_line_options options(argc, argv);
   if (options.must_print_model_info())
      CNE6SSM_info::print(std::cout);
   if (options.must_exit())
      return options.status();

   const std::string rgflow_file(options.get_rgflow_file());
   const std::string slha_input_file(options.get_slha_input_file());
   const std::string slha_output_file(options.get_slha_output_file());
   const std::string spectrum_file(options.get_spectrum_file());
   CNE6SSM_slha_io slha_io;
   Spectrum_generator_settings spectrum_generator_settings;
   QedQcd oneset;
   CNE6SSM_input_parameters input;

   if (slha_input_file.empty()) {
      ERROR("No SLHA input file given!\n"
            "   Please provide one via the option --slha-input-file=");
      return EXIT_FAILURE;
   }

   try {
      slha_io.read_from_file(slha_input_file);
      slha_io.fill(oneset);
      slha_io.fill(input);
      slha_io.fill(spectrum_generator_settings);
   } catch (const Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }

   oneset.toMz(); // run SM fermion masses to MZ

      CNE6SSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(
      spectrum_generator_settings.get(Spectrum_generator_settings::precision));
   spectrum_generator.set_max_iterations(
      spectrum_generator_settings.get(Spectrum_generator_settings::max_iterations));
   spectrum_generator.set_calculate_sm_masses(
      spectrum_generator_settings.get(Spectrum_generator_settings::calculate_sm_masses) >= 1.0);
   spectrum_generator.set_alternate_ewsb(
      spectrum_generator_settings.get(Spectrum_generator_settings::alternate_ewsb) >= 1.0);
   spectrum_generator.set_input_scale(
      slha_io.get_input_scale());
   spectrum_generator.set_parameter_output_scale(
      slha_io.get_parameter_output_scale());
   spectrum_generator.set_pole_mass_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   spectrum_generator.set_ewsb_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::ewsb_loop_order));
   spectrum_generator.set_beta_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::beta_loop_order));
   spectrum_generator.set_threshold_corrections(
      spectrum_generator_settings.get(Spectrum_generator_settings::threshold_corrections));

   spectrum_generator.run(oneset, input);

   const CNE6SSM<algorithm_type>& model
      = spectrum_generator.get_model();
   const Problems<CNE6SSM_info::NUMBER_OF_PARTICLES>& problems
      = spectrum_generator.get_problems();

   // output
   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(oneset);
   slha_io.set_minpar(input);
   slha_io.set_extpar(input);
   if (!problems.have_serious_problem())
      slha_io.set_spectrum(model);

   if (!slha_output_file.empty()) {
      slha_io.write_to_file(slha_output_file);
   }

   if (!spectrum_file.empty())
      spectrum_generator.write_spectrum(spectrum_file);

   if (!rgflow_file.empty())
      spectrum_generator.write_running_couplings(rgflow_file);

   const int exit_code = spectrum_generator.get_exit_code();

   // if point is valid, calculate coefficients in RG running
   if (exit_code == EXIT_SUCCESS) {
      // soft scalar masses to calculate coefficients for
      std::vector<CNE6SSM_info::Parameters> soft_scalar_masses 
         = {CNE6SSM_info::mHd2, CNE6SSM_info::mHu2, CNE6SSM_info::ms2,
            CNE6SSM_info::msbar2, CNE6SSM_info::mphi2, CNE6SSM_info::mq200,
            CNE6SSM_info::mu200};

      // soft gaugino masses to calculate coefficients for
      std::vector<CNE6SSM_info::Parameters> soft_gaugino_masses
         = {CNE6SSM_info::MassB, CNE6SSM_info::MassWB, CNE6SSM_info::MassG, 
            CNE6SSM_info::MassBp};

      // soft trilinears to calculate coefficients for
      std::vector<CNE6SSM_info::Parameters> soft_trilinears
         = {CNE6SSM_info::TYu22, CNE6SSM_info::TYu00, CNE6SSM_info::TSigmax,
            CNE6SSM_info::TLambdax};

      // get high scale
      double high_scale = spectrum_generator.get_high_scale();

      // get SUSY scale - this will be fixed scale
      // for comparison at
      double susy_scale = spectrum_generator.get_susy_scale();

      // calculate coefficients and percentage errors the original way
      std::map<CNE6SSM_info::Parameters, std::vector<double> > soft_scalar_mass_coeffs;
      std::map<CNE6SSM_info::Parameters, double> soft_scalar_mass_errors;
      for (std::vector<CNE6SSM_info::Parameters>::const_iterator it = soft_scalar_masses.begin(),
              end = soft_scalar_masses.end(); it != end; ++it) {
         const Eigen::Array<double,4,1> coeffs = model.get_soft_scalar_mass_coeffs(*it, susy_scale, high_scale);
         soft_scalar_mass_coeffs[*it] = {coeffs(0), coeffs(1), coeffs(2), coeffs(3)};
         const double pred_value = coeffs(0) * Sqr(input.m0) + coeffs(1) * Sqr(input.m12)
            + coeffs(2) * input.m12 * input.Azero + coeffs(3) * Sqr(input.Azero);
         soft_scalar_mass_errors[*it] = 100.0 * Abs((model.get_parameter(*it) - pred_value) / 
                                                    (0.5 * (model.get_parameter(*it) + pred_value)));
      }

      std::map<CNE6SSM_info::Parameters, std::vector<double> > soft_gaugino_mass_coeffs;
      std::map<CNE6SSM_info::Parameters, double> soft_gaugino_mass_errors;
      for (std::vector<CNE6SSM_info::Parameters>::const_iterator it = soft_gaugino_masses.begin(),
              end = soft_gaugino_masses.end(); it != end; ++it) {
         const Eigen::Array<double,2,1> coeffs = model.get_soft_gaugino_mass_coeffs(*it, susy_scale, high_scale);
         soft_gaugino_mass_coeffs[*it] = {coeffs(0), coeffs(1)};
         const double pred_value = coeffs(0) * input.Azero + coeffs(1) * input.m12;
         soft_gaugino_mass_errors[*it] = 100.0 * Abs((model.get_parameter(*it) - pred_value) / 
                                                     (0.5 * (model.get_parameter(*it) + pred_value)));
      }

      std::map<CNE6SSM_info::Parameters, std::vector<double> > soft_trilinear_coeffs;
      std::map<CNE6SSM_info::Parameters, double> soft_trilinear_errors;
      for (std::vector<CNE6SSM_info::Parameters>::const_iterator it = soft_trilinears.begin(),
              end = soft_trilinears.end(); it != end; ++it) {
         const Eigen::Array<double,2,1> coeffs = model.get_soft_trilinear_coeffs(*it, susy_scale, high_scale);
         soft_trilinear_coeffs[*it] = {coeffs(0), coeffs(1)};
         const double pred_value = coeffs(0) * input.Azero + coeffs(1) * input.m12;
         soft_trilinear_errors[*it] = 100.0 * Abs((model.get_parameter(*it) - pred_value) /
                                                 (0.5 * (model.get_parameter(*it) + pred_value)));
      }

      // print results
      std::cout << "Coefficients for requested parameters:\n";
      for (std::vector<CNE6SSM_info::Parameters>::const_iterator it = soft_scalar_masses.begin(),
              end = soft_scalar_masses.end(); it != end; ++it) {
         std::cout << "a(" << CNE6SSM_info::parameter_names[*it] 
                   << ", " << susy_scale << " GeV) = " << soft_scalar_mass_coeffs[*it][0] << "\n";
         std::cout << "b(" << CNE6SSM_info::parameter_names[*it] 
                   << ", " << susy_scale << " GeV) = " << soft_scalar_mass_coeffs[*it][1] << "\n";
         std::cout << "c(" << CNE6SSM_info::parameter_names[*it] 
                   << ", " << susy_scale << " GeV) = " << soft_scalar_mass_coeffs[*it][2] << "\n";
         std::cout << "d(" << CNE6SSM_info::parameter_names[*it] 
                   << ", " << susy_scale << " GeV) = " << soft_scalar_mass_coeffs[*it][3] << "\n";
         std::cout << "% error = " << soft_scalar_mass_errors[*it] << "\n";
      }

      for (std::vector<CNE6SSM_info::Parameters>::const_iterator it = soft_gaugino_masses.begin(),
              end = soft_gaugino_masses.end(); it != end; ++it) {
         std::cout << "p(" << CNE6SSM_info::parameter_names[*it] 
                   << ", " << susy_scale << " GeV) = " << soft_gaugino_mass_coeffs[*it][0] << "\n";
         std::cout << "q(" << CNE6SSM_info::parameter_names[*it] 
                   << ", " << susy_scale << " GeV) = " << soft_gaugino_mass_coeffs[*it][1] << "\n";
         std::cout << "% error = " << soft_scalar_mass_errors[*it] << "\n";
      }

      for (std::vector<CNE6SSM_info::Parameters>::const_iterator it = soft_trilinears.begin(),
              end = soft_trilinears.end(); it != end; ++it) {
         std::cout << "e(" << CNE6SSM_info::parameter_names[*it] 
                   << ", " << susy_scale << " GeV) = " << soft_trilinear_coeffs[*it][0] << "\n";
         std::cout << "f(" << CNE6SSM_info::parameter_names[*it] 
                   << ", " << susy_scale << " GeV) = " << soft_trilinear_coeffs[*it][1] << "\n";
         std::cout << "% error = " << soft_trilinear_errors[*it] << "\n";
      }

   }

   return exit_code;
}
