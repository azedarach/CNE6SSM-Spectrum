// ====================================================================
// Does a grid scan of the CNE6SSM parameter space
//   - trial version before generalising it to arbitrary input
//     parameters (by which I mean arbitrary scans over the 
//     parameters defined in CNE6SSM_input_parameters)
// ====================================================================

#include "CNE6SSM_input_parameters.hpp"
#include "CNE6SSM_spectrum_generator.hpp"

#include "command_line_options.hpp"
#include "error.hpp"
#include "grid_scanner.hpp"
#include "scan.hpp"
#include "lowe.h"

#include <iostream>

namespace flexiblesusy {

   void print_usage()
   {
      std::cout <<
         "Usage: gridscan_CNE6SSM.x [options]\n"
         "Options:\n"
         "  --m0=<value>\n"
         "  --m12=<value>\n"
         "  --TanBeta=<value>\n"
         "  --SignLambdax=<value>\n"
         "  --Azero=<value>\n"
         
         "  --help,-h                         print this help message"
                << std::endl;
   }
   
   void set_command_line_parameters(int argc, char* argv[],
                                    CNE6SSM_input_parameters& input)
   {
      for (int i = 1; i < argc; ++i) {
         const std::string option(argv[i]);
         
         if (Command_line_options::get_parameter_value(option, "--m0=", input.m0))
            continue;
         
         if(Command_line_options::get_parameter_value(option, "--m12=", input.m12))
            continue;
         
         if(Command_line_options::get_parameter_value(option, "--TanBeta=", input.TanBeta))
            continue;
         
         if(Command_line_options::get_parameter_value(option, "--SignLambdax=", input.SignLambdax))
            continue;
         
         if(Command_line_options::get_parameter_value(option, "--Azero=", input.Azero))
            continue;
         
         
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
      if (is_zero(input.TanBeta))
         input.TanBeta = 10.0;

      if (is_zero(input.SignLambdax))
         input.SignLambdax = 1;

      if (is_zero(input.Azero))
         input.Azero = 1000.;

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
      
      input.KappaInput(0,0) = 2.0e-1;
      input.KappaInput(0,1) = 0.;
      input.KappaInput(0,2) = 0.;
      input.KappaInput(1,0) = 0.;
      input.KappaInput(1,1) = 2.0e-1;
      input.KappaInput(1,2) = 0.;
      input.KappaInput(2,0) = 0.;
      input.KappaInput(2,1) = 0.;
      input.KappaInput(2,2) = 2.0e-1;
      
      input.Lambda12Input(0,0) = 5.0e-1;
      input.Lambda12Input(0,1) = 0.;
      input.Lambda12Input(1,0) = 0.;
      input.Lambda12Input(1,1) = 5.0e-1;
      
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

   void set_soft_mass_values(const std::vector<std::size_t>& posn, const std::vector<std::size_t>& dims, CNE6SSM_input_parameters& input)
   {
      const double m0_lower = 0.; // GeV
      const double m0_upper = 3000.; // GeV
      const double m12_lower = 0.; // GeV
      const double m12_upper = 3000.; // GeV
      
      double m0_incr = 0.0;
      double m12_incr = 0.0;
      
      if (dims.at(0) > 1) 
         m0_incr = (m0_upper - m0_lower) / (dims.at(0) - 1.0);
      if (dims.at(1) > 1) 
         m12_incr = (m12_upper - m12_lower) / (dims.at(1) - 1.0);
      
      input.m0 = m0_lower + m0_incr * posn.at(0);
      input.m12 = m12_lower + m12_incr * posn.at(1);
      
   }

} // namespace flexiblesusy

int main(int argc, char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;

   CNE6SSM_input_parameters input;
   set_command_line_parameters(argc, argv, input);
   set_default_parameter_values(input);

   QedQcd oneset;
   oneset.toMz();

   CNE6SSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(5.0e-3);
   spectrum_generator.set_max_iterations(0);   // 0 == automatic
   spectrum_generator.set_calculate_sm_masses(0); // 0 == no
   spectrum_generator.set_alternate_ewsb(1); // 1 == yes
   spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale
   
   std::size_t m0_npts = 45;
   std::size_t m12_npts = 45;

   std::vector<std::size_t> scan_dimensions = {m0_npts, m12_npts};

   Grid_scanner scan(scan_dimensions);

   cout << "# "
        << std::setw(12) << std::left << "m0" << ' '
        << std::setw(12) << std::left << "m12" << ' '
        << std::setw(12) << std::left << "TanBeta" << ' '
        << std::setw(12) << std::left << "Azero" << ' '
        << std::setw(12) << std::left << "SignLambdax" << ' '
        << std::setw(12) << std::left << "Mhh(1)/GeV" << ' '
        << std::setw(12) << std::left << "error"
        << '\n';

   while (!scan.has_finished()) {
      set_soft_mass_values(scan.get_position(), scan.get_dimensions(), input);

      spectrum_generator.run(oneset, input);

      const CNE6SSM<algorithm_type>& model = spectrum_generator.get_model();
      const CNE6SSM_physical& pole_masses = model.get_physical();
      const Problems<CNE6SSM_info::NUMBER_OF_PARTICLES>& problems
         = spectrum_generator.get_problems();
      const double higgs = pole_masses.Mhh(0);
      const bool error = problems.have_serious_problem();

      cout << " "
           << std::setw(12) << std::left << input.m0 << ' '
           << std::setw(12) << std::left << input.m12 << ' '
           << std::setw(12) << std::left << input.TanBeta << ' '
           << std::setw(12) << std::left << input.Azero << ' '
           << std::setw(12) << std::left << input.SignLambdax << ' '
           << std::setw(12) << std::left << higgs << ' '
           << std::setw(12) << std::left << error;
      if (error) {
         cout << "\t# " << problems;
      }
      cout << '\n';
      scan.step_forward();
   }

   return 0;
}

