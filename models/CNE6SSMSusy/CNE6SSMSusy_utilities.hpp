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

#ifndef CNE6SSMSusy_UTILITIES_H
#define CNE6SSMSusy_UTILITIES_H

#include "CNE6SSMSusy_susy_parameters.hpp"
#include "CNE6SSMSusy_info.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>
#include <string>
#include <vector>
#include <valarray>
#include <utility>

#define MODELPARAMETER(p) model.get_##p()

namespace flexiblesusy {

class CNE6SSMSusy_parameter_getter {
public:
   Eigen::ArrayXd get_parameters(const CNE6SSMSusy_susy_parameters& model) {
      return model.get();
   }
   std::vector<std::string> get_parameter_names(const CNE6SSMSusy_susy_parameters&) const {
      using namespace CNE6SSMSusy_info;
      return std::vector<std::string>(parameter_names,
                                      parameter_names + NUMBER_OF_PARAMETERS);
   }
};

} // namespace flexiblesusy

#endif
