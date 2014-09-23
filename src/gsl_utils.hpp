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

#ifndef GSL_UTILS_H
#define GSL_UTILS_H

#include <cstddef>
#include <vector>
#include <gsl/gsl_vector.h>

namespace flexiblesusy {

/// Returns true if GSL vector contains a Nan, false otherwise
bool contains_nan(const gsl_vector*, std::size_t);

/// DH:: overloaded version doing the same for a std::vector
bool contains_nan(const std::vector<double>&, std::size_t);

}

#endif
