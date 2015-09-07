// ====================================================================
// Program to calculate pole masses and mixings of a CNE6SSM parameter
// space point using a semianalytic solver
// ====================================================================

#include "CNE6SSM_semi_two_scale_input_parameters.hpp"
#include "CNE6SSM_slha_io.hpp"
#include "CNE6SSM_semianalytic_spectrum_generator.hpp"
// note
#include "CNE6SSM_higgs_upper_bound.hpp"

#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "command_line_options.hpp"

#include <iostream>
#include <cstdlib>

int main(int argc, const char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
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
   CNE6SSM_semianalytic_input_parameters<algorithm_type> input;

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

   CNE6SSM_semianalytic_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(
      spectrum_generator_settings.get(Spectrum_generator_settings::precision));
   spectrum_generator.set_max_iterations(
      spectrum_generator_settings.get(Spectrum_generator_settings::max_iterations));
   spectrum_generator.set_calculate_sm_masses(
      spectrum_generator_settings.get(Spectrum_generator_settings::calculate_sm_masses) >= 1.0);
   spectrum_generator.set_force_output(
      spectrum_generator_settings.get(Spectrum_generator_settings::force_output) >= 1.0);
   spectrum_generator.set_parameter_output_scale(
      slha_io.get_parameter_output_scale());
   spectrum_generator.set_pole_mass_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   spectrum_generator.set_ewsb_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::ewsb_loop_order));
   spectrum_generator.set_beta_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::beta_loop_order));
   spectrum_generator.set_threshold_corrections_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::threshold_corrections_loop_order));
   spectrum_generator.set_two_loop_corrections(
      spectrum_generator_settings.get_two_loop_corrections());
   // note
   spectrum_generator.set_ewsb_iteration_precision(0.5);

   spectrum_generator.run(oneset, input);

   const CNE6SSM_semianalytic_slha<algorithm_type> model(spectrum_generator.get_model());
   const Problems<CNE6SSM_info::NUMBER_OF_PARTICLES>& problems
      = spectrum_generator.get_problems();

   CNE6SSM_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale = spectrum_generator.get_low_scale();

   // output
   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(oneset);
   slha_io.set_minpar(input);
   slha_io.set_extpar(input);
   if (!problems.have_problem() ||
       spectrum_generator_settings.get(Spectrum_generator_settings::force_output)) {
      slha_io.set_spectrum(model);
      slha_io.set_extra(model, scales);
   }

   if (slha_output_file.empty()) {
      slha_io.write_to_stream(std::cout);
   } else {
      slha_io.write_to_file(slha_output_file);
   }

   if (!spectrum_file.empty())
      spectrum_generator.write_spectrum(spectrum_file);

   if (!rgflow_file.empty())
      spectrum_generator.write_running_couplings(rgflow_file);

   // note
   CNE6SSM_higgs_upper_bound upper_bound(model);

   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   std::cout << "tree level upper bound = "
             << upper_bound.calculate_tree_level_upper_bound() << '\n';

   std::cout << "1-lp upper bound (only stops) = "
             << upper_bound.calculate_one_loop_upper_bound() << '\n';

   upper_bound.set_include_down_tadpoles(true);

   std::cout << "1-lp upper bound (stops + sbottoms) = "
             << upper_bound.calculate_one_loop_upper_bound() << '\n';

   upper_bound.set_include_all_SM_generations(true);
   upper_bound.set_include_down_tadpoles(false);

   std::cout << "1-lp upper bound (all up squarks) = "
             << upper_bound.calculate_one_loop_upper_bound() << '\n';

   upper_bound.set_include_down_tadpoles(true);

   std::cout << "1-lp upper bound (all up squarks + all down squarks) = "
             << upper_bound.calculate_one_loop_upper_bound() << '\n';

   upper_bound.set_include_exotic_tadpoles(true);

   std::cout << "1-lp upper bound (3Su + 3Sd + 3SDX) = "
             << upper_bound.calculate_one_loop_upper_bound() << '\n';

   upper_bound.set_include_inert_singlet_tadpoles(true);

   std::cout << "1-lp upper bound (3Su + 3Sd + 3SDX + 3S) = "
             << upper_bound.calculate_one_loop_upper_bound() << '\n';

   upper_bound.set_include_inert_neutral_higgs_tadpoles(true);

   std::cout << "1-lp upper bound (3Su + 3Sd + 3SDX + 3S + 2HI0) = "
             << upper_bound.calculate_one_loop_upper_bound() << '\n';

   upper_bound.set_include_inert_charged_higgs_tadpoles(true);

   std::cout << "1-lp upper bound (3Su + 3Sd + 3SDX + 3S + 2HI0 + 2HIPM) = "
             << upper_bound.calculate_one_loop_upper_bound() << '\n';

   double up_contribution = 0.;
   double unrotated_up_00 = 0.;
   double unrotated_up_01 = 0.;
   double unrotated_up_11 = 0.;
   double unrotated_down_00 = 0.;
   double unrotated_down_01 = 0.;
   double unrotated_down_11 = 0.;
   double exotic_contribution = 0.;
   double unrotated_exotic_00 = 0.;
   double unrotated_exotic_01 = 0.;
   double unrotated_exotic_11 = 0.;
   for (unsigned gen = 0; gen < 3; ++gen) {
      up_contribution += upper_bound.get_up_contribution(gen);
      unrotated_up_00 += upper_bound.get_unrotated_up_contribution(gen, 0, 0);
      unrotated_up_01 += upper_bound.get_unrotated_up_contribution(gen, 0, 1);
      unrotated_up_11 += upper_bound.get_unrotated_up_contribution(gen, 1, 1);
      unrotated_down_00 += upper_bound.get_unrotated_down_contribution(gen, 0, 0);
      unrotated_down_01 += upper_bound.get_unrotated_down_contribution(gen, 0, 1);
      unrotated_down_11 += upper_bound.get_unrotated_down_contribution(gen, 1, 1);
      exotic_contribution += upper_bound.get_exotic_contribution(gen);
      unrotated_exotic_00 += upper_bound.get_unrotated_exotic_contribution(gen, 0, 0);
      unrotated_exotic_01 += upper_bound.get_unrotated_exotic_contribution(gen, 0, 1);
      unrotated_exotic_11 += upper_bound.get_unrotated_exotic_contribution(gen, 1, 1);
   }

   std::cout << "up contribution to upper bound = "
             << up_contribution << '\n';

   std::cout << "rotated 0p Su self-energy(0,0) = "
             << upper_bound.get_up_self_energy(0., 0, 0) << '\n';

   std::cout << "unrotated up contribution (0,0) = "
             << unrotated_up_00 << '\n';
   std::cout << "unrotated up contribution (0,1) = "
             << unrotated_up_01 << '\n';
   std::cout << "unrotated up contribution (1,1) = "
             << unrotated_up_11 << '\n';

   std::cout << "unrotated up 0p self-energy(0,0) = "
             << upper_bound.get_unrotated_up_self_energy(0., 0, 0) << '\n';
   std::cout << "unrotated up 0p self-energy(0,1) = "
             << upper_bound.get_unrotated_up_self_energy(0., 0, 1) << '\n';
   std::cout << "unrotated up 0p self-energy(1,1) = "
             << upper_bound.get_unrotated_up_self_energy(0., 1, 1) << '\n';

   std::cout << "unrotated down contribution (0,0) = "
             << unrotated_down_00 << '\n';
   std::cout << "unrotated down contribution (0,1) = "
             << unrotated_down_01 << '\n';
   std::cout << "unrotated down contribution (1,1) = "
             << unrotated_down_11 << '\n';

   std::cout << "unrotated down 0p self-energy(0,0) = "
             << upper_bound.get_unrotated_down_self_energy(0., 0, 0) << '\n';
   std::cout << "unrotated down 0p self-energy(0,1) = "
             << upper_bound.get_unrotated_down_self_energy(0., 0, 1) << '\n';
   std::cout << "unrotated down 0p self-energy(1,1) = "
             << upper_bound.get_unrotated_down_self_energy(0., 1, 1) << '\n';

   std::cout << "unrotated exotic contribution (0,0) = "
             << unrotated_exotic_00 << '\n';
   std::cout << "unrotated exotic contribution (0,1) = "
             << unrotated_exotic_01 << '\n';
   std::cout << "unrotated exotic contribution (1,1) = "
             << unrotated_exotic_11 << '\n';

   std::cout << "unrotated exotic 0p self-energy(0,0) = "
             << upper_bound.get_unrotated_exotic_self_energy(0., 0, 0) << '\n';
   std::cout << "unrotated exotic 0p self-energy(0,1) = "
             << upper_bound.get_unrotated_exotic_self_energy(0., 0, 1) << '\n';
   std::cout << "unrotated exotic 0p self-energy(1,1) = "
             << upper_bound.get_unrotated_exotic_self_energy(0., 1, 1) << '\n';

   std::cout << "exotic contribution to upper bound = "
             << exotic_contribution << '\n';

   std::cout << "rotated 0p SDX self-energy(0,0) = "
             << upper_bound.get_exotic_self_energy(0., 0, 0) << '\n';

   std::cout << "unrotated Su self-energy (0,0) = "
             << upper_bound.get_unrotated_up_self_energy(model.get_Mhh()(0), 0, 0) << '\n';

   std::cout << "unrotated Su self-energy (1,1) = "
             << upper_bound.get_unrotated_up_self_energy(model.get_Mhh()(0), 1, 1) << '\n';

   std::cout << "unrotated SDX self-energy (0,0) = "
             << upper_bound.get_unrotated_exotic_self_energy(model.get_Mhh()(0), 0, 0) << '\n';

   std::cout << "unrotated SDX self-energy (1,1) = "
             << upper_bound.get_unrotated_exotic_self_energy(model.get_Mhh()(0), 1, 1) << '\n';

   std::cout << "unrotated self-energy (0,0) = "
             << upper_bound.get_unrotated_full_self_energy(model.get_Mhh()(0), 0, 0) << '\n';

   std::cout << "unrotated self-energy (1,1) = "
             << upper_bound.get_unrotated_full_self_energy(model.get_Mhh()(0), 1, 1) << '\n';

   std::cout << "full self-energy (1,1) element = "
             << upper_bound.get_full_self_energy(model.get_Mhh()(0), 0, 0) << '\n';

   std::cout << "self-energy (1,1) element = "
             << upper_bound.get_self_energy(model.get_Mhh()(0), 0, 0) << '\n';

   upper_bound.set_include_all_SM_generations(false);
   upper_bound.set_include_up_tadpoles(false);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(true);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   std::cout << "t_1 = " << upper_bound.get_tadpole_vd() << '\n';
   std::cout << "t_2 = " << upper_bound.get_tadpole_vu() << '\n';
   std::cout << "tadpole contribution = " << (upper_bound.get_tadpole_vd() / model.get_vd())
      * Sqr(Cos(ArcTan(model.get_vu() / model.get_vd()))) + (upper_bound.get_tadpole_vu() / model.get_vu())
      * Sqr(Sin(ArcTan(model.get_vu() / model.get_vd()))) << '\n';
   std::cout << "Delta_{11}^D = " << upper_bound.get_exotic_contribution(0) +
      upper_bound.get_exotic_contribution(1) + upper_bound.get_exotic_contribution(2) << '\n';
   std::cout << "sum = " << (upper_bound.get_tadpole_vd() / model.get_vd())
      * Sqr(Cos(ArcTan(model.get_vu() / model.get_vd()))) + (upper_bound.get_tadpole_vu() / model.get_vu())
      * Sqr(Sin(ArcTan(model.get_vu() / model.get_vd()))) + upper_bound.get_exotic_contribution(0) +
      upper_bound.get_exotic_contribution(1) + upper_bound.get_exotic_contribution(2) << '\n';

   upper_bound.set_include_all_SM_generations(false);
   upper_bound.set_include_up_tadpoles(true);
   upper_bound.set_include_down_tadpoles(false);
   upper_bound.set_include_exotic_tadpoles(false);
   upper_bound.set_include_inert_singlet_tadpoles(false);
   upper_bound.set_include_inert_neutral_higgs_tadpoles(false);
   upper_bound.set_include_inert_charged_higgs_tadpoles(false);

   std::cout << "t_1 = " << upper_bound.get_tadpole_vd() << '\n';
   std::cout << "t_2 = " << upper_bound.get_tadpole_vu() << '\n';
   std::cout << "tadpole contribution = " << (upper_bound.get_tadpole_vd() / model.get_vd())
      * Sqr(Cos(ArcTan(model.get_vu() / model.get_vd()))) + (upper_bound.get_tadpole_vu() / model.get_vu())
      * Sqr(Sin(ArcTan(model.get_vu() / model.get_vd()))) << '\n';
   std::cout << "Delta_{11}^D = " << upper_bound.get_up_contribution(2) << '\n';
   std::cout << "sum = " << (upper_bound.get_tadpole_vd() / model.get_vd())
      * Sqr(Cos(ArcTan(model.get_vu() / model.get_vd()))) + (upper_bound.get_tadpole_vu() / model.get_vu())
      * Sqr(Sin(ArcTan(model.get_vu() / model.get_vd()))) + upper_bound.get_up_contribution(2) << '\n';


   const int exit_code = spectrum_generator.get_exit_code();

   return exit_code;
}
