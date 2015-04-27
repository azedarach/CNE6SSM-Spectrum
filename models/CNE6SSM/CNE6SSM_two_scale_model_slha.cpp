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
 * @file CNE6SSM_two_scale_model_slha.cpp
 * @brief CNE6SSM model class wrapper for SLHA conversion
 */

#include "CNE6SSM_two_scale_model_slha.hpp"
#include "slha_io.hpp"
#include "ckm.hpp"
#include "pmns.hpp"

namespace flexiblesusy {

#define CLASSNAME CNE6SSM_slha<Two_scale>
#define LOCALPHYSICAL(p) physical.p

CLASSNAME::CNE6SSM_slha(const CNE6SSM_input_parameters<Two_scale>& input_)
   : CNE6SSM<Two_scale>(input_)
   , physical_slha()
   , ckm(Eigen::Matrix<std::complex<double>,3,3>::Identity())
   , pmns(Eigen::Matrix<std::complex<double>,3,3>::Identity())
{
}

/**
 * Copy constructor.  Copies from base class (two-scale model class in
 * BPMZ convention) and converts parameters to SLHA.
 *
 * @param model_ model class in BPMZ convention
 */
CLASSNAME::CNE6SSM_slha(const CNE6SSM<Two_scale>& model_)
   : CNE6SSM<Two_scale>(model_)
{
   convert_to_slha();
}

CLASSNAME::~CNE6SSM_slha()
{
}

void CLASSNAME::clear()
{
   CNE6SSM<Two_scale>::clear();
   physical_slha.clear();
}

void CLASSNAME::calculate_spectrum()
{
   CNE6SSM<Two_scale>::calculate_spectrum();
   convert_to_slha();
}

void CLASSNAME::convert_to_slha()
{
   physical_slha = get_physical();

   convert_to_slha(physical_slha);

   convert_yukawa_couplings_to_slha();
   calculate_ckm_matrix();
   calculate_pmns_matrix();
   convert_trilinear_couplings_to_slha();
   convert_soft_squared_masses_to_slha();
}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 *
 * @param physical struct of physical parameters to convert
 */
void CLASSNAME::convert_to_slha(CNE6SSM_physical& physical)
{
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChiI), LOCALPHYSICAL(ZNI));
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChiP), LOCALPHYSICAL(ZNp));

}

void CLASSNAME::calculate_ckm_matrix()
{
   ckm = ZUL_slha * ZDL_slha.adjoint();
   CKM_parameters::to_pdg_convention(ckm, ZUL_slha, ZDL_slha, ZUR_slha, ZDR_slha);

}

void CLASSNAME::calculate_pmns_matrix()
{
   pmns << 1, 0, 0, 0, 1, 0, 0, 0, 1;

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

/**
 * Convert trilinear couplings to SLHA convention
 */
void CLASSNAME::convert_trilinear_couplings_to_slha()
{
   TYu_slha = (ZUR_slha.conjugate() * TYu * ZUL_slha.adjoint()).real();
   TYd_slha = (ZDR_slha.conjugate() * TYd * ZDL_slha.adjoint()).real();
   TYe_slha = (ZER_slha.conjugate() * TYe * ZEL_slha.adjoint()).real();

}

/**
 * Convert trilinear couplings to SLHA convention
 */
void CLASSNAME::convert_soft_squared_masses_to_slha()
{
   mq2_slha = (ZDL_slha * mq2 * ZDL_slha.adjoint()).real();
   mu2_slha = (ZUR_slha.conjugate() * mu2 * ZUR_slha.transpose()).real();
   md2_slha = (ZDR_slha.conjugate() * md2 * ZDR_slha.transpose()).real();
   ml2_slha = (ZEL_slha * ml2 * ZEL_slha.adjoint()).real();
   me2_slha = (ZER_slha.conjugate() * me2 * ZER_slha.transpose()).real();

}

const CNE6SSM_physical& CLASSNAME::get_physical_slha() const
{
   return physical_slha;
}

CNE6SSM_physical& CLASSNAME::get_physical_slha()
{
   return physical_slha;
}

void CLASSNAME::print(std::ostream& ostr) const
{
   CNE6SSM<Two_scale>::print(ostr);

   ostr << "----------------------------------------\n"
           "SLHA convention:\n"
           "----------------------------------------\n";
   physical_slha.print(ostr);
}

} // namespace flexiblesusy
