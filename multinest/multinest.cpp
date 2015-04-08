// ====================================================================
// Test implementation of using MultiNest with FlexibleSUSY
// ====================================================================

#include "CNE6SSM_input_parameters.hpp"
#include "CNE6SSM_slha_io.hpp"
#include "CNE6SSM_spectrum_generator.hpp"

#include "error.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "command_line_options.hpp"
#include "wrappers.hpp"

#include <iostream>
#include <cstdlib>

#include "multinest.h"

using namespace flexiblesusy;

void print_usage()
{
   std::cout << 
      "Usage: multinest.x [options]\n"
      "Options:\n"
      "  --root-file-name=<name>           base file name for output\n"
      "  --help,-h                         print this help message"
             << std::endl;
}

void get_command_line_arguments(int argc, char* argv[], std::string& base_name)
{
   for (int i = 1; i < argc; ++i) {
      const std::string option(argv[i]);

      if (Command_line_options::starts_with(option,"--root-file-name=")) {
         base_name = Command_line_options::get_option_value(option, "=");
         if (base_name.empty())
            WARNING("no root file name given");
         continue;
      }

      if (option == "--help" || option == "-h") {
         print_usage();
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option);
      exit(EXIT_FAILURE);
   }
}

void set_default_parameter_values(CNE6SSM_input_parameters& input)
{   
   input.ssumInput = 40000.0; // GeV
   input.QS = 5.;
   
   input.hEInput(0,0) = 0.;
   input.hEInput(0,1) = 0.;
   input.hEInput(1,0) = 0.;
   input.hEInput(1,1) = 0.;
   input.hEInput(2,0) = 0.;
   input.hEInput(2,1) = 0.;
   
   input.SigmaLInput = 3.0e-1;
   input.KappaPrInput = 2.0e-2;
   input.SigmaxInput = 1.0e-1;
   
   input.gDInput(0,0) = 0.;
   input.gDInput(0,1) = 0.;
   input.gDInput(0,2) = 0.;
   input.gDInput(1,0) = 0.;
   input.gDInput(1,1) = 0.;
   input.gDInput(1,2) = 0.;
   input.gDInput(2,0) = 0.;
   input.gDInput(2,1) = 0.;
   input.gDInput(2,2) = 0.;
   
   input.KappaInput(0,1) = 0.;
   input.KappaInput(0,2) = 0.;
   input.KappaInput(1,0) = 0.;
   input.KappaInput(1,2) = 0.;
   input.KappaInput(2,0) = 0.;
   input.KappaInput(2,1) = 0.;
   
   input.Lambda12Input(0,1) = 0.;
   input.Lambda12Input(1,0) = 0.;
   
   input.fuInput(0,0) = 1.0e-7;
   input.fuInput(0,1) = 0.;
   input.fuInput(1,0) = 0.;
   input.fuInput(1,1) = 1.0e-7;
   input.fuInput(2,0) = 1.0e-7;
   input.fuInput(2,1) = 0.;
   
   input.fdInput(0,0) = 1.0e-7;
   input.fdInput(0,1) = 0.;
   input.fdInput(1,0) = 0.;
   input.fdInput(1,1) = 1.0e-7;
   input.fdInput(2,0) = 0.;
   input.fdInput(2,1) = 1.0e-7;
   
   input.MuPrInput = 1.0e4;
   input.MuPhiInput = 0.;
   input.BMuPrInput = 1.0e4;
   input.BMuPhiInput = 0.;
   
}

// Map the unit hypercube to the required prior distribution of the physical
// parameters
void transform_hypercube_to_priors(double *cube, int num_dims, int num_params, void *parameters)
{
   const double m0_lower = 0.;
   const double m0_upper = 5000.;
   const double m12_lower = 0.;
   const double m12_upper = 5000.;
   const double TanBeta_lower = 5.;
   const double TanBeta_upper = 40.;
   const double SignLambdax_lower = -1.;
   const double SignLambdax_upper = 1.;
   const double Azero_lower = -15000.;
   const double Azero_upper = 15000.;
   const double Kappa_lower = -3.0;
   const double Kappa_upper = 3.0;
   const double Lambda12_lower = -3.0;
   const double Lambda12_upper = 3.0;
   
   cube[0] = m0_lower + cube[0] * (m0_upper - m0_lower);
   cube[1] = m12_lower + cube[1] * (m12_upper - m12_lower);
   cube[2] = TanBeta_lower + cube[2] * (TanBeta_upper - TanBeta_lower);
   cube[3] = SignLambdax_lower + cube[3] * (SignLambdax_upper - SignLambdax_lower);
   cube[4] = Azero_lower + cube[4] * (Azero_upper - Azero_lower);
   cube[5] = Kappa_lower + cube[5] * (Kappa_upper - Kappa_lower);
   cube[6] = Lambda12_lower + cube[6] * (Lambda12_upper - Lambda12_lower);
}

// Calculate the log likelihood using the physical parameters found in the first
// argument
double calculate_log_likelihood(double *cube, int &num_dims, int &num_params, void *parameters)
{
   CNE6SSM_input_parameters input;
   set_default_parameter_values(input);
   input.m0 = cube[0];
   input.m12 = cube[1];
   input.TanBeta = cube[2];
   input.SignLambdax = Sign(cube[3]);
   input.Azero = cube[4];
   input.KappaInput(0,0) = cube[5];
   input.KappaInput(1,1) = cube[5];
   input.KappaInput(2,2) = cube[5];
   input.Lambda12Input(0,0) = cube[6];
   input.Lambda12Input(1,1) = cube[6];

   QedQcd oneset;
   oneset.toMz();

   CNE6SSM_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_precision_goal(1.0e-3);
   spectrum_generator.set_max_iterations(0);
   spectrum_generator.set_calculate_sm_masses(0);
   spectrum_generator.set_parameter_output_scale(0);
   spectrum_generator.set_alternate_ewsb(1);
   spectrum_generator.set_pole_mass_loop_order(2);
   spectrum_generator.set_ewsb_loop_order(2);
   spectrum_generator.set_beta_loop_order(2);
   spectrum_generator.set_threshold_corrections(1);

   spectrum_generator.run(oneset, input);

   const CNE6SSM<Two_scale>& model = spectrum_generator.get_model();
   const Problems<CNE6SSM_info::NUMBER_OF_PARTICLES>& problems
      = spectrum_generator.get_problems();

   const double flag_log_likelihood = -2.0e90;
   const double inverse_width = 0.5; // (standard deviation)^{-1}

   if (problems.have_serious_problem()) {
      cube[7] = 0.;
      cube[8] = 0.;
      cube[9] = 1.0e3;
      cube[10] = 0.;

      return flag_log_likelihood;
   } else {
      const CNE6SSM_physical& pole_masses = model.get_physical();

      // output ratio and masses, and Higgs mass
      cube[7] = pole_masses.MSu(5);
      cube[8] = pole_masses.MGlu;
      cube[9] = pole_masses.MSu(5) / pole_masses.MGlu;
      cube[10] = pole_masses.Mhh(0);

      // use MSu(5) as an upper bound on the ratio
      return -Sqr(inverse_width * pole_masses.MSu(5) / pole_masses.MGlu);
   }
}

// Function that is passed to MultiNest for computing the likelihood
void get_log_likelihood(double *unit_hypercube, int &num_dims, int &num_params, double &log_likelihood, void *parameters)
{
   transform_hypercube_to_priors(unit_hypercube, num_dims, num_params, parameters);
   log_likelihood = calculate_log_likelihood(unit_hypercube, num_dims, num_params, parameters);
}

// Routine called by MultiNest useful for displaying statistics
void dump_statistics(int &num_posterior_samples, int &total_num_live, int &num_params, 
                     double** physical_live_points, double** posterior_points, 
                     double** parameter_stats, double &max_log_likelihood, double &log_evidence,
                     double &log_evidence_INS, double &log_evidence_error, void* parameters)
{
   // currently unused
}

struct Multinest_options {
   bool do_nested_importance_samp;
   bool do_mode_separation;
   bool use_const_efficiency;
   int num_live_points;
   double efficiency;
   double tolerance;
   int num_iters_per_update;
   double min_log_evidence;
   int expected_max_num_modes;
   int random_seed;
   bool print_to_stdout;
   bool resume_from_previous;
   bool do_write_output_files;
   bool multinest_init_mpi;
   double min_log_likelihood;
   int max_num_iters;
};

void run_multinest(int num_free_params, int num_derived_params, int num_mode_sep_params,
                   int* periodic_bcs, const Multinest_options& options, const std::string& base_name,
                   void* parameters = 0)
{
   const int num_params = num_free_params + num_derived_params;

   nested::run(options.do_nested_importance_samp, options.do_mode_separation, 
               options.use_const_efficiency, options.num_live_points,
               options.tolerance, options.efficiency, num_free_params, 
               num_params, num_mode_sep_params, options.expected_max_num_modes, 
               options.num_iters_per_update, options.min_log_evidence, base_name,
               options.random_seed, periodic_bcs, options.print_to_stdout, 
               options.resume_from_previous, options.do_write_output_files, 
               options.multinest_init_mpi, options.min_log_likelihood, 
               options.max_num_iters, get_log_likelihood, dump_statistics, 
               parameters);
}

int main(int argc, char* argv[])
{
   std::string base_name;
   get_command_line_arguments(argc, argv, base_name);

   if (base_name.empty()) {
      WARNING("No root file name given!\n"
              "   Using default root: ./chains/multinest-\n"
              "   You can provide one via the option --root-file-name=");
      base_name = "./chains/multinest-";
   }

   // number of parameters of interest in the problem,
   // specifying which have periodic boundary conditions
   const int num_free_params = 7;
   const int num_derived_params = 4;
   const int num_mode_sep_params = num_free_params;
   
   int periodic_bcs[num_free_params];
   for (int i = 0; i < num_free_params; ++i) periodic_bcs[i] = 0;

   // configuration options for MultiNest
   Multinest_options options;
   options.do_nested_importance_samp = true;
   options.do_mode_separation = false;
   options.use_const_efficiency = false;
   options.num_live_points = 1000;
   options.efficiency = 0.5;
   options.tolerance = 0.5;
   options.num_iters_per_update = 1000;
   options.min_log_evidence = -1.0e90;
   options.expected_max_num_modes = 100;
   options.random_seed = -1;
   options.print_to_stdout = true;
   options.resume_from_previous = true;
   options.do_write_output_files = true;
   options.multinest_init_mpi = true;
   options.min_log_likelihood = -1.0e90;
   options.max_num_iters = 0;
   
   run_multinest(num_free_params, num_derived_params, num_mode_sep_params,
                 periodic_bcs, options, base_name);

   return 0;
}
