// ====================================================================
// Does a grid scan of the CNE6SSM parameter space
// ====================================================================

#include "CNE6SSM_input_parameters.hpp"
#include "CNE6SSM_spectrum_generator.hpp"

#include "error.hpp"
#include "grid_scanner.hpp"
#include "scan.hpp"
#include "lowe.h"

#include <iostream>

void get_current_inputs(const std::vector<std::size_t>&, flexiblesusy::CNE6SSM_input_parameters&);

int main()
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;

   CNE6SSM_input_parameters input;
   QedQcd oneset;
   oneset.toMz();

   CNE6SSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(1.0e-4);
   spectrum_generator.set_max_iterations(0);   // 0 == automatic
   spectrum_generator.set_calculate_sm_masses(0); // 0 == no
   spectrum_generator.set_alternate_ewsb(1); // 1 == yes
   spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale

   std::size_t m0_npts = 5;
   std::size_t m12_npts = 2;
   std::size_t A0_npts = 1;

   std::vector<std::size_t> scan_dimensions = {m0_npts, m12_npts, A0_npts};

   Grid_scanner scan(scan_dimensions);

   cout << "# "
        << std::setw(12) << std::left << "m0" << ' '
        << std::setw(12) << std::left << "m12" << ' '
        << std::setw(12) << std::left << "TanBeta" << ' '
        << std::setw(12) << std::left << "Azero" << ' '
        << std::setw(12) << std::left << "Mhh(1)/GeV" << ' '
        << std::setw(12) << std::left << "error"
        << '\n';

   while (!scan.has_finished()) {
      get_current_inputs(scan.get_position(), input);

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

void get_current_inputs(const std::vector<std::size_t>& posn, flexiblesusy::CNE6SSM_input_parameters&)
{
   
}
