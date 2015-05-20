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

// File generated at Sun 19 Apr 2015 20:37:13

/**
 * @file CNE6SSMSusy_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver
 *
 * This file was generated at Sun 19 Apr 2015 20:37:13 with FlexibleSUSY
 * 1.0.4 (git commit: v1.0.3-cpc6-749-g227a308) and SARAH 4.5.3 .
 */

#ifndef CNE6SSMSusy_TWO_SCALE_H
#define CNE6SSMSusy_TWO_SCALE_H

#include "CNE6SSMSusy_info.hpp"
#include "CNE6SSMSusy_model.hpp"
#include "CNE6SSMSusy_susy_parameters.hpp"
#include "CNE6SSMSusy_two_scale_input_parameters.hpp"
#include "two_scale_model.hpp"
#include "problems.hpp"
#include "config.h"

#include <iosfwd>
#include <string>

#include <gsl/gsl_vector.h>
#include <Eigen/Core>

namespace flexiblesusy {

class EWSB_solver;
class Two_scale;
/**
 * @class CNE6SSMSusy<Two_scale>
 * @brief model class with routines for solving boundary value problem
 */
template<>
class CNE6SSMSusy<Two_scale> : public Two_scale_model, public CNE6SSMSusy_susy_parameters {
public:
   explicit CNE6SSMSusy(const CNE6SSMSusy_input_parameters<Two_scale>& input_ = CNE6SSMSusy_input_parameters<Two_scale>());
   virtual ~CNE6SSMSusy();

   const Problems<CNE6SSMSusy_info::NUMBER_OF_PARTICLES>& get_problems() const;
   Problems<CNE6SSMSusy_info::NUMBER_OF_PARTICLES>& get_problems();

   virtual void clear();
   const CNE6SSMSusy_input_parameters<Two_scale>& get_input() const;
   void set_input_parameters(const CNE6SSMSusy_input_parameters<Two_scale>&);

   // interface functions
   virtual void calculate_spectrum();
   virtual void clear_problems();
   virtual std::string name() const;
   virtual void run_to(double scale, double eps = -1.0);
   virtual void print(std::ostream&) const;
   virtual void set_precision(double);

   // model interface
   double get_parameter(unsigned) const;
   void set_parameter(unsigned, double);

private:
   // input parameters
   CNE6SSMSusy_input_parameters<Two_scale> input;
   double precision;              ///< RG running precision
   Problems<CNE6SSMSusy_info::NUMBER_OF_PARTICLES> problems;
};

std::ostream& operator<<(std::ostream&, const CNE6SSMSusy<Two_scale>&);

} // namespace flexiblesusy

#endif
