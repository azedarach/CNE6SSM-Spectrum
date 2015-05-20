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

/**
 * @file CNE6SSMSusy_two_scale_model_slha.hpp
 * @brief contains wrapper class for model class in SLHA convention
 */

// File generated at Sun 19 Apr 2015 20:31:41

#ifndef CNE6SSMSusy_TWO_SCALE_SLHA_H
#define CNE6SSMSusy_TWO_SCALE_SLHA_H

#include "CNE6SSMSusy_two_scale_model.hpp"
#include "CNE6SSMSusy_model_slha.hpp"

namespace flexiblesusy {

class Two_scale;

/**
 * @class CNE6SSMSusy_slha<Two_scale>
 * @brief model class wrapper in SLHA convention
 */

template<>
class CNE6SSMSusy_slha<Two_scale> : public CNE6SSMSusy<Two_scale> {
public:
   explicit CNE6SSMSusy_slha(const CNE6SSMSusy_input_parameters<Two_scale>& input_
                             = CNE6SSMSusy_input_parameters<Two_scale>());
   explicit CNE6SSMSusy_slha(const CNE6SSMSusy<Two_scale>&);
   virtual ~CNE6SSMSusy_slha();

   virtual void clear();
   void convert_to_slha(); ///< converts couplings to SLHA convention

   // interface functions
   virtual void calculate_spectrum();
   virtual void print(std::ostream&) const;

   const Eigen::Array<double,3,1>& get_Yu_slha() const { return Yu_slha; }
   double get_Yu_slha(int i) const { return Yu_slha(i); }
   const Eigen::Array<double,3,1>& get_Yd_slha() const { return Yd_slha; }
   double get_Yd_slha(int i) const { return Yd_slha(i); }
   const Eigen::Array<double,3,1>& get_Ye_slha() const { return Ye_slha; }
   double get_Ye_slha(int i) const { return Ye_slha(i); }

private:
   Eigen::Array<double,3,1> Yu_slha;
   Eigen::Array<double,3,1> Yd_slha;
   Eigen::Array<double,3,1> Ye_slha;

   Eigen::Matrix<double,3,3> ZUR_slha;
   Eigen::Matrix<double,3,3> ZUL_slha;
   Eigen::Matrix<double,3,3> ZDR_slha;
   Eigen::Matrix<double,3,3> ZDL_slha;
   Eigen::Matrix<double,3,3> ZER_slha;
   Eigen::Matrix<double,3,3> ZEL_slha;

   void convert_yukawa_couplings_to_slha();
};

} // namespace flexiblesusy

#endif
