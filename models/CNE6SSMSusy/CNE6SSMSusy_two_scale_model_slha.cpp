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

// File generated at Sun 19 Apr 2015 20:31:41

/**
 * @file CNE6SSMSusy_two_scale_model_slha.cpp
 * @brief CNE6SSMSusy model class wrapper for SLHA conversion
 */

#include "CNE6SSMSusy_two_scale_model_slha.hpp"
#include "slha_io.hpp"
#include "linalg2.hpp"

namespace flexiblesusy {

#define CLASSNAME CNE6SSMSusy_slha<Two_scale>
#define LOCALPHYSICAL(p) physical.p

CLASSNAME::CNE6SSMSusy_slha(const CNE6SSMSusy_input_parameters<Two_scale>& input_)
   : CNE6SSMSusy<Two_scale>(input_)
   , ZUR_slha(Eigen::Matrix<double,3,3>::Zero())
   , ZUL_slha(Eigen::Matrix<double,3,3>::Zero())
   , ZDR_slha(Eigen::Matrix<double,3,3>::Zero())
   , ZDL_slha(Eigen::Matrix<double,3,3>::Zero())
   , ZER_slha(Eigen::Matrix<double,3,3>::Zero())
   , ZEL_slha(Eigen::Matrix<double,3,3>::Zero())
{
}

/**
 * Copy constructor.  Copies from base class (two-scale model class in
 * BPMZ convention) and converts parameters to SLHA.
 *
 * @param model_ model class in BPMZ convention
 */
CLASSNAME::CNE6SSMSusy_slha(const CNE6SSMSusy<Two_scale>& model_)
   : CNE6SSMSusy<Two_scale>(model_)
{
   convert_to_slha();
}

CLASSNAME::~CNE6SSMSusy_slha()
{
}

void CLASSNAME::clear()
{
   CNE6SSMSusy<Two_scale>::clear();
}

void CLASSNAME::calculate_spectrum()
{
   CNE6SSMSusy<Two_scale>::calculate_spectrum();
   convert_to_slha();
}

void CLASSNAME::convert_to_slha()
{
   convert_yukawa_couplings_to_slha();
}

/**
 * Convert Yukawa couplings to SLHA convention
 */
void CLASSNAME::convert_yukawa_couplings_to_slha()
{
   fs_svd(Yu, Yu_slha, ZUR_slha, ZUL_slha);
   fs_svd(Yd, Yd_slha, ZDR_slha, ZDL_slha);
   fs_svd(Ye, Ye_slha, ZER_slha, ZEL_slha);

}

void CLASSNAME::print(std::ostream& ostr) const
{
   CNE6SSMSusy<Two_scale>::print(ostr);
}

} // namespace flexiblesusy
