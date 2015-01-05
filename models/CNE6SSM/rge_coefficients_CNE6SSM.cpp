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
      // get high scale
      double high_scale = spectrum_generator.get_high_scale();

      // get SUSY scale - this will be fixed scale
      // for comparison at
      double susy_scale = spectrum_generator.get_susy_scale();

      CNE6SSM<algorithm_type> running_model(model);

      running_model.run_to(high_scale);

      // all though in principle 6 terms, 2 are known to vanish
      const std::size_t num_terms = 4;
      
      Eigen::Matrix<double, num_terms, num_terms> input_values;
      Eigen::VectorXd mHu2_values(num_terms);
      Eigen::VectorXd mHd2_values(num_terms);
      Eigen::VectorXd ms2_values(num_terms);

      double m0_centre = model.get_input().m0;
      double m12_centre = model.get_input().m12;
      double Azero_centre = model.get_input().Azero;

      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::default_random_engine generator(seed);
      std::normal_distribution<double> m0_distribution(m0_centre, 0.1 * m0_centre);
      std::normal_distribution<double> m12_distribution(m12_centre, 0.1 * m12_centre);
      std::normal_distribution<double> Azero_distribution(Azero_centre, 0.1 * Azero_centre);

      for (std::size_t i = 0; i < num_terms; ++i) {
         // generate random values from a Gaussian centred
         // on the initial values
         double m0_tmp = m0_distribution(generator);
         double m12_tmp = m12_distribution(generator);
         double Azero_tmp = Azero_distribution(generator);

         input_values(i, 0) = Sqr(m0_tmp);
         input_values(i, 1) = Sqr(m12_tmp);
         input_values(i, 2) = m12_tmp * Azero_tmp;
         input_values(i, 3) = Sqr(Azero_tmp);

         // get the soft mass values
         mHu2_values(i) = get_parameter_from_inputs(running_model, CNE6SSM_info::mHu2, m0_tmp, m12_tmp, Azero_tmp, susy_scale);
         mHd2_values(i) = get_parameter_from_inputs(running_model, CNE6SSM_info::mHd2, m0_tmp, m12_tmp, Azero_tmp, susy_scale);
         ms2_values(i) = get_parameter_from_inputs(running_model, CNE6SSM_info::ms2, m0_tmp, m12_tmp, Azero_tmp, susy_scale);
      }

      // solve for the coefficients in the expansion of mHd2 and mHu2
      Eigen::VectorXd mHd2_coeffs = input_values.fullPivHouseholderQr().solve(mHd2_values);
      Eigen::VectorXd mHu2_coeffs = input_values.fullPivHouseholderQr().solve(mHu2_values);
      Eigen::VectorXd ms2_coeffs = input_values.fullPivHouseholderQr().solve(ms2_values);

      std::cout << "Estimate for coefficients:\n";
      std::cout << "aHd(" << susy_scale << " GeV) = " << mHd2_coeffs(0) << "\n";
      std::cout << "aHu(" << susy_scale << " GeV) = " << mHu2_coeffs(0) << "\n";
      std::cout << "bHd(" << susy_scale << " GeV) = " << mHd2_coeffs(1) << "\n";
      std::cout << "bHu(" << susy_scale << " GeV) = " << mHu2_coeffs(1) << "\n";
      std::cout << "cHd(" << susy_scale << " GeV) = " << mHd2_coeffs(2) << "\n";
      std::cout << "cHu(" << susy_scale << " GeV) = " << mHu2_coeffs(2) << "\n";
      std::cout << "dHd(" << susy_scale << " GeV) = " << mHd2_coeffs(3) << "\n";
      std::cout << "dHu(" << susy_scale << " GeV) = " << mHu2_coeffs(3) << "\n";
      std::cout << "aS1(" << susy_scale << " GeV) = " << ms2_coeffs(0) << "\n";
      std::cout << "bS1(" << susy_scale << " GeV) = " << ms2_coeffs(1) << "\n";
      std::cout << "cS1(" << susy_scale << " GeV) = " << ms2_coeffs(2) << "\n";
      std::cout << "dS1(" << susy_scale << " GeV) = " << ms2_coeffs(3) << "\n";
      std::cout << "aLambdax(" << susy_scale << " GeV) = " 
                << get_tree_level_Lambdax_soft_term(model, mHd2_coeffs(0), mHu2_coeffs(0)) << "\n";
      std::cout << "bLambdax(" << susy_scale << " GeV) = " 
                << get_tree_level_Lambdax_soft_term(model, mHd2_coeffs(1), mHu2_coeffs(1)) << "\n";
      std::cout << "cLambdax(" << susy_scale << " GeV) = " 
                << get_tree_level_Lambdax_soft_term(model, mHd2_coeffs(2), mHu2_coeffs(2)) << "\n";
      std::cout << "dLambdax(" << susy_scale << " GeV) = " 
                << get_tree_level_Lambdax_soft_term(model, mHd2_coeffs(3), mHu2_coeffs(3)) << "\n";
      std::cout << "lLambdax(" << susy_scale << " GeV) = " << get_tree_level_Lambdax_constant_term(model) << "\n";

   }

   return exit_code;
}
