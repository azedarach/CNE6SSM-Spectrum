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
#include <chrono>
#include <sys/time.h>

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

double get_wall_time()
{
   struct timeval time;
   if (gettimeofday(&time,NULL)) {
      return 0;
   }
   return (double)time.tv_sec + (double)time.tv_usec*0.000001;
}

double get_cpu_time()
{
   return (double)clock() / CLOCKS_PER_SEC;
}

int main(int argc, char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;
   typedef std::chrono::duration<int,std::micro> microseconds_t;
   std::chrono::high_resolution_clock::time_point start_point = std::chrono::high_resolution_clock::now();
   double wall_start = get_wall_time();
   double cpu_start = get_cpu_time();

   CNE6SSM_input_parameters input;
   set_command_line_parameters(argc, argv, input);
   set_default_parameter_values(input);

   QedQcd oneset;
   oneset.toMz();

   CNE6SSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(1.0e-3);
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
        << std::setw(12) << std::left << "MGlu/GeV" << ' '
        << std::setw(12) << std::left << "MChaP/GeV" << ' '
        << std::setw(12) << std::left << "MVZp/GeV" << ' '
        << std::setw(12) << std::left << "MSd(1)/GeV" << ' '
        << std::setw(12) << std::left << "MSd(2)/GeV" << ' '
        << std::setw(12) << std::left << "MSd(3)/GeV" << ' '
        << std::setw(12) << std::left << "MSd(4)/GeV" << ' '
        << std::setw(12) << std::left << "MSd(5)/GeV" << ' '
        << std::setw(12) << std::left << "MSd(6)/GeV" << ' '
        << std::setw(12) << std::left << "MSv(1)/GeV" << ' '
        << std::setw(12) << std::left << "MSv(2)/GeV" << ' '
        << std::setw(12) << std::left << "MSv(3)/GeV" << ' '
        << std::setw(12) << std::left << "MSu(1)/GeV" << ' '
        << std::setw(12) << std::left << "MSu(2)/GeV" << ' '
        << std::setw(12) << std::left << "MSu(3)/GeV" << ' '
        << std::setw(12) << std::left << "MSu(4)/GeV" << ' '
        << std::setw(12) << std::left << "MSu(5)/GeV" << ' '
        << std::setw(12) << std::left << "MSu(6)/GeV" << ' '
        << std::setw(12) << std::left << "MSe(1)/GeV" << ' '
        << std::setw(12) << std::left << "MSe(2)/GeV" << ' '
        << std::setw(12) << std::left << "MSe(3)/GeV" << ' '
        << std::setw(12) << std::left << "MSe(4)/GeV" << ' '
        << std::setw(12) << std::left << "MSe(5)/GeV" << ' '
        << std::setw(12) << std::left << "MSe(6)/GeV" << ' '
        << std::setw(12) << std::left << "MSDX(1)/GeV" << ' '
        << std::setw(12) << std::left << "MSDX(2)/GeV" << ' '
        << std::setw(12) << std::left << "MSDX(3)/GeV" << ' '
        << std::setw(12) << std::left << "MSDX(4)/GeV" << ' '
        << std::setw(12) << std::left << "MSDX(5)/GeV" << ' '
        << std::setw(12) << std::left << "MSDX(6)/GeV" << ' '
        << std::setw(12) << std::left << "Mhh(1)/GeV" << ' '
        << std::setw(12) << std::left << "Mhh(2)/GeV" << ' '
        << std::setw(12) << std::left << "Mhh(3)/GeV" << ' '
        << std::setw(12) << std::left << "Mhh(4)/GeV" << ' '
        << std::setw(12) << std::left << "Mhh(5)/GeV" << ' '
        << std::setw(12) << std::left << "MAh(1)/GeV" << ' '
        << std::setw(12) << std::left << "MAh(2)/GeV" << ' '
        << std::setw(12) << std::left << "MAh(3)/GeV" << ' '
        << std::setw(12) << std::left << "MAh(4)/GeV" << ' '
        << std::setw(12) << std::left << "MAh(5)/GeV" << ' '
        << std::setw(12) << std::left << "MHpm(1)/GeV" << ' '
        << std::setw(12) << std::left << "MHpm(2)/GeV" << ' '
        << std::setw(12) << std::left << "MChi(1)/GeV" << ' '
        << std::setw(12) << std::left << "MChi(2)/GeV" << ' '
        << std::setw(12) << std::left << "MChi(3)/GeV" << ' '
        << std::setw(12) << std::left << "MChi(4)/GeV" << ' '
        << std::setw(12) << std::left << "MChi(5)/GeV" << ' '
        << std::setw(12) << std::left << "MChi(6)/GeV" << ' '
        << std::setw(12) << std::left << "MChi(7)/GeV" << ' '
        << std::setw(12) << std::left << "MChi(8)/GeV" << ' '
        << std::setw(12) << std::left << "MCha(1)/GeV" << ' '
        << std::setw(12) << std::left << "MCha(2)/GeV" << ' '
        << std::setw(12) << std::left << "MFDX(1)/GeV" << ' '
        << std::setw(12) << std::left << "MFDX(2)/GeV" << ' '
        << std::setw(12) << std::left << "MFDX(3)/GeV" << ' '
        << std::setw(12) << std::left << "MSHI0(1)/GeV" << ' '
        << std::setw(12) << std::left << "MSHI0(2)/GeV" << ' '
        << std::setw(12) << std::left << "MSHI0(3)/GeV" << ' '
        << std::setw(12) << std::left << "MSHI0(4)/GeV" << ' '
        << std::setw(12) << std::left << "MSHI0(5)/GeV" << ' '
        << std::setw(12) << std::left << "MSHI0(6)/GeV" << ' '
        << std::setw(12) << std::left << "MSHI0(7)/GeV" << ' '
        << std::setw(12) << std::left << "MSHIPM(1)/GeV" << ' '
        << std::setw(12) << std::left << "MSHIPM(2)/GeV" << ' '
        << std::setw(12) << std::left << "MSHIPM(3)/GeV" << ' '
        << std::setw(12) << std::left << "MSHIPM(4)/GeV" << ' '
        << std::setw(12) << std::left << "MChaI(1)/GeV" << ' '
        << std::setw(12) << std::left << "MChaI(2)/GeV" << ' '
        << std::setw(12) << std::left << "MChiI(1)/GeV" << ' '
        << std::setw(12) << std::left << "MChiI(2)/GeV" << ' '
        << std::setw(12) << std::left << "MChiI(3)/GeV" << ' '
        << std::setw(12) << std::left << "MChiI(4)/GeV" << ' '
        << std::setw(12) << std::left << "MChiI(5)/GeV" << ' '
        << std::setw(12) << std::left << "MChiI(6)/GeV" << ' '
        << std::setw(12) << std::left << "MChiI(7)/GeV" << ' '
        << std::setw(12) << std::left << "MSHp0(1)/GeV" << ' '
        << std::setw(12) << std::left << "MSHp0(2)/GeV" << ' '
        << std::setw(12) << std::left << "MSHpp(1)/GeV" << ' '
        << std::setw(12) << std::left << "MSHpp(2)/GeV" << ' '
        << std::setw(12) << std::left << "MChiP(1)/GeV" << ' '
        << std::setw(12) << std::left << "MChiP(2)/GeV" << ' '
        << std::setw(12) << std::left << "error"
        << '\n';

   while (!scan.has_finished()) {
      set_soft_mass_values(scan.get_position(), scan.get_dimensions(), input);

      spectrum_generator.run(oneset, input);

      const CNE6SSM<algorithm_type>& model = spectrum_generator.get_model();
      const CNE6SSM_physical& pole_masses = model.get_physical();
      const Problems<CNE6SSM_info::NUMBER_OF_PARTICLES>& problems
         = spectrum_generator.get_problems();
      const bool error = problems.have_serious_problem();

      cout << " "
           << std::setw(12) << std::left << input.m0 << ' '
           << std::setw(12) << std::left << input.m12 << ' '
           << std::setw(12) << std::left << input.TanBeta << ' '
           << std::setw(12) << std::left << input.Azero << ' '
           << std::setw(12) << std::left << input.SignLambdax << ' '
           << std::setw(12) << std::left << pole_masses.MGlu << ' '
           << std::setw(12) << std::left << pole_masses.MChaP << ' '
           << std::setw(12) << std::left << pole_masses.MVZp << ' '
           << std::setw(12) << std::left << pole_masses.MSd(0) << ' '
           << std::setw(12) << std::left << pole_masses.MSd(1) << ' '
           << std::setw(12) << std::left << pole_masses.MSd(2) << ' '
           << std::setw(12) << std::left << pole_masses.MSd(3) << ' '
           << std::setw(12) << std::left << pole_masses.MSd(4) << ' '
           << std::setw(12) << std::left << pole_masses.MSd(5) << ' '
           << std::setw(12) << std::left << pole_masses.MSv(0) << ' '
           << std::setw(12) << std::left << pole_masses.MSv(1) << ' '
           << std::setw(12) << std::left << pole_masses.MSv(2) << ' '
           << std::setw(12) << std::left << pole_masses.MSu(0) << ' '
           << std::setw(12) << std::left << pole_masses.MSu(1) << ' '
           << std::setw(12) << std::left << pole_masses.MSu(2) << ' '
           << std::setw(12) << std::left << pole_masses.MSu(3) << ' '
           << std::setw(12) << std::left << pole_masses.MSu(4) << ' '
           << std::setw(12) << std::left << pole_masses.MSu(5) << ' '
           << std::setw(12) << std::left << pole_masses.MSe(0) << ' '
           << std::setw(12) << std::left << pole_masses.MSe(1) << ' '
           << std::setw(12) << std::left << pole_masses.MSe(2) << ' '
           << std::setw(12) << std::left << pole_masses.MSe(3) << ' '
           << std::setw(12) << std::left << pole_masses.MSe(4) << ' '
           << std::setw(12) << std::left << pole_masses.MSe(5) << ' '
           << std::setw(12) << std::left << pole_masses.MSDX(0) << ' '
           << std::setw(12) << std::left << pole_masses.MSDX(1) << ' '
           << std::setw(12) << std::left << pole_masses.MSDX(2) << ' '
           << std::setw(12) << std::left << pole_masses.MSDX(3) << ' '
           << std::setw(12) << std::left << pole_masses.MSDX(4) << ' '
           << std::setw(12) << std::left << pole_masses.MSDX(5) << ' '
           << std::setw(12) << std::left << pole_masses.Mhh(0) << ' '
           << std::setw(12) << std::left << pole_masses.Mhh(1) << ' '
           << std::setw(12) << std::left << pole_masses.Mhh(2) << ' '
           << std::setw(12) << std::left << pole_masses.Mhh(3) << ' '
           << std::setw(12) << std::left << pole_masses.Mhh(4) << ' '
           << std::setw(12) << std::left << pole_masses.MAh(0) << ' '
           << std::setw(12) << std::left << pole_masses.MAh(1) << ' '
           << std::setw(12) << std::left << pole_masses.MAh(2) << ' '
           << std::setw(12) << std::left << pole_masses.MAh(3) << ' '
           << std::setw(12) << std::left << pole_masses.MAh(4) << ' '
           << std::setw(12) << std::left << pole_masses.MHpm(0) << ' '
           << std::setw(12) << std::left << pole_masses.MHpm(1) << ' '
           << std::setw(12) << std::left << pole_masses.MChi(0) << ' '
           << std::setw(12) << std::left << pole_masses.MChi(1) << ' '
           << std::setw(12) << std::left << pole_masses.MChi(2) << ' '
           << std::setw(12) << std::left << pole_masses.MChi(3) << ' '
           << std::setw(12) << std::left << pole_masses.MChi(4) << ' '
           << std::setw(12) << std::left << pole_masses.MChi(5) << ' '
           << std::setw(12) << std::left << pole_masses.MChi(6) << ' '
           << std::setw(12) << std::left << pole_masses.MChi(7) << ' '
           << std::setw(12) << std::left << pole_masses.MCha(0) << ' '
           << std::setw(12) << std::left << pole_masses.MCha(1) << ' '
           << std::setw(12) << std::left << pole_masses.MFDX(0) << ' '
           << std::setw(12) << std::left << pole_masses.MFDX(1) << ' '
           << std::setw(12) << std::left << pole_masses.MFDX(2) << ' '
           << std::setw(12) << std::left << pole_masses.MSHI0(0) << ' '
           << std::setw(12) << std::left << pole_masses.MSHI0(1) << ' '
           << std::setw(12) << std::left << pole_masses.MSHI0(2) << ' '
           << std::setw(12) << std::left << pole_masses.MSHI0(3) << ' '
           << std::setw(12) << std::left << pole_masses.MSHI0(4) << ' '
           << std::setw(12) << std::left << pole_masses.MSHI0(5) << ' '
           << std::setw(12) << std::left << pole_masses.MSHI0(6) << ' '
           << std::setw(12) << std::left << pole_masses.MSHIPM(0) << ' '
           << std::setw(12) << std::left << pole_masses.MSHIPM(1) << ' '
           << std::setw(12) << std::left << pole_masses.MSHIPM(2) << ' '
           << std::setw(12) << std::left << pole_masses.MSHIPM(3) << ' '
           << std::setw(12) << std::left << pole_masses.MChaI(0) << ' '
           << std::setw(12) << std::left << pole_masses.MChaI(1) << ' '
           << std::setw(12) << std::left << pole_masses.MChiI(0) << ' '
           << std::setw(12) << std::left << pole_masses.MChiI(1) << ' '
           << std::setw(12) << std::left << pole_masses.MChiI(2) << ' '
           << std::setw(12) << std::left << pole_masses.MChiI(3) << ' '
           << std::setw(12) << std::left << pole_masses.MChiI(4) << ' '
           << std::setw(12) << std::left << pole_masses.MChiI(5) << ' '
           << std::setw(12) << std::left << pole_masses.MChiI(6) << ' '
           << std::setw(12) << std::left << pole_masses.MSHp0(0) << ' '
           << std::setw(12) << std::left << pole_masses.MSHp0(1) << ' '
           << std::setw(12) << std::left << pole_masses.MSHpp(0) << ' '
           << std::setw(12) << std::left << pole_masses.MSHpp(1) << ' '
           << std::setw(12) << std::left << pole_masses.MChiP(0) << ' '
           << std::setw(12) << std::left << pole_masses.MChiP(1) << ' '
           << std::setw(12) << std::left << error;
      if (error) {
         cout << "\t# " << problems;
      }
      cout << '\n';
      scan.step_forward();
   }

   std::chrono::high_resolution_clock::time_point end_point = std::chrono::high_resolution_clock::now();
   microseconds_t duration(std::chrono::duration_cast<microseconds_t>(end_point - start_point));
   double time_in_seconds = duration.count() * 0.000001;
   double wall_end = get_wall_time();
   double cpu_end = get_cpu_time();

   cout << "# Scan completed in " << time_in_seconds << " seconds\n";
   cout << "# Wall time = " << wall_end - wall_start << " seconds\n";
   cout << "# CPU time  = " << cpu_end - cpu_start << " seconds\n";

   return 0;
}

