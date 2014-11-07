// ====================================================================
// Routines used in scan class
// ====================================================================

#ifndef CNE6SSM_SCAN_UTILITIES_H
#define CNE6SSM_SCAN_UTILITIES_H

#include "CNE6SSM_two_scale_model.hpp"
#include "CNE6SSM_info.hpp"
#include "CNE6SSM_utilities.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>
#include <string>
#include <vector>
#include <valarray>
#include <utility>

#define PHYSICAL(p) model.get_physical().p
#define MODELPARAMETER(p) model.get_##p()


namespace flexiblesusy {

void write_CNE6SSM_inputs(const CNE6SSM_input_parameters&, std::ostream &, std::size_t);
void write_CNE6SSM_inputs_list(std::ostream &, std::size_t);
std::valarray<double> to_valarray(double);

template <class Scalar, int M, int N>
std::valarray<double> to_valarray(const Eigen::Array<Scalar, M, N>& v)
{
   return std::valarray<double>(v.data(), v.size());
}

template <int M, int N>
std::valarray<double> to_valarray(const Eigen::Matrix<double, M, N>& v)
{
   return std::valarray<double>(v.data(), v.size());
}

// this needs to extract real and imaginary parts, then
// interleave them in the valarray
template <int M, int N>
std::valarray<double> to_valarray(const Eigen::Matrix<std::complex<double>, M, N> & v)
{
   std::size_t num_entries = v.size();
   Eigen::Matrix<double,M,N> re_v = v.real();
   std::valarray<double> real_parts(re_v.data(), num_entries);
   Eigen::Matrix<double,M,N> im_v = v.imag();
   std::valarray<double>imag_parts(im_v.data(), num_entries);

   // combine into a single valarray
   std::valarray<double> result(2 * v.size());
   result[std::slice(0, num_entries, 1)] = real_parts;
   result[std::slice(num_entries, num_entries, 1)] = imag_parts;
   return result;
}

class CNE6SSM_pole_mass_writer {
public:
   CNE6SSM_pole_mass_writer();
   ~CNE6SSM_pole_mass_writer() {}

   template <class T>
   void extract_pole_masses(const CNE6SSM<T>&);
   void write_pole_masses_comment_line(std::ostream &) const;
   void write_pole_masses_line(std::ostream &) const;

private:
   struct TPoleMass {
      std::string name;
      std::valarray<double> masses;
      TPoleMass(const std::string& name_, const std::valarray<double>& masses_)
         : name(name_)
         , masses(masses_)
         {}
   };
   typedef std::vector<TPoleMass> TPoleMasses;
   TPoleMasses pole_masses;

   CNE6SSM_input_parameters pole_masses_inputs;
   Problems<CNE6SSM_info::NUMBER_OF_PARTICLES> pole_masses_problems;

   double pole_masses_scale;
   unsigned width;
};

class CNE6SSM_drbar_values_writer {
public:
   CNE6SSM_drbar_values_writer();
   ~CNE6SSM_drbar_values_writer() {}

   template <class T>
   void extract_drbar_masses(const CNE6SSM<T>&);
   template <class T>
   void extract_drbar_susy_pars(const CNE6SSM<T>&);
   template <class T>
   void extract_drbar_soft_pars(const CNE6SSM<T>&);
   template <class T>
   void extract_drbar_mixings(const CNE6SSM<T>&);

   void write_drbar_masses_comment_line(std::ostream &) const;
   void write_drbar_susy_pars_comment_line(std::ostream &) const;
   void write_drbar_soft_pars_comment_line(std::ostream &) const;
   void write_drbar_mixings_comment_line(std::ostream &) const;

   void write_drbar_masses_line(std::ostream &) const;
   void write_drbar_susy_pars_line(std::ostream &) const;
   void write_drbar_soft_pars_line(std::ostream &) const;
   void write_drbar_mixings_line(std::ostream &) const;

private:
   struct TMass {
      std::string name;
      std::valarray<double> masses;
      TMass(const std::string& name_, const std::valarray<double>& masses_)
         : name(name_)
         , masses(masses_)
         {}
   };

   // note: stored in column-major format, i.e.
   // elements are (0,0), (1,0), (2,0),...
   struct TParameter {
      std::string name;
      std::size_t mass_dimension;
      std::size_t rows;
      std::size_t cols;
      std::valarray<double> values;
      TParameter(const std::string& name_, 
                 std::size_t mass_dimension_,
                 std::size_t rows_, std::size_t cols_,
                 const std::valarray<double>& values_)
         : name(name_)
         , mass_dimension(mass_dimension_)
         , rows(rows_)
         , cols(cols_)
         , values(values_)
         {}
   };

   // note: stored in column-major format, i.e.
   // elements are (0,0), (1,0), (2,0),...
   struct TMixing {
      std::string name;
      std::size_t dimension;
      bool is_real;
      std::valarray<double> mixings;
      TMixing(const std::string name_, std::size_t dimension_,
              bool is_real_,
              const std::valarray<double>& mixings_)
         : name(name_)
         , dimension(dimension_)
         , is_real(is_real_)
         , mixings(mixings_)
         {}
   };

   typedef std::vector<TMass> TRunningMasses;
   typedef std::vector<TParameter> TSusyPars;
   typedef std::vector<TParameter> TSoftPars;
   typedef std::vector<TMixing> TMixings;

   TRunningMasses drbar_masses;
   TSusyPars drbar_susy_pars;
   TSoftPars drbar_soft_pars;
   TMixings drbar_mixings;

   CNE6SSM_input_parameters drbar_masses_inputs;
   CNE6SSM_input_parameters drbar_susy_pars_inputs;
   CNE6SSM_input_parameters drbar_soft_pars_inputs;
   CNE6SSM_input_parameters drbar_mixings_inputs;

   Problems<CNE6SSM_info::NUMBER_OF_PARTICLES> drbar_masses_problems;
   Problems<CNE6SSM_info::NUMBER_OF_PARTICLES> drbar_susy_pars_problems;
   Problems<CNE6SSM_info::NUMBER_OF_PARTICLES> drbar_soft_pars_problems;
   Problems<CNE6SSM_info::NUMBER_OF_PARTICLES> drbar_mixings_problems;

   double drbar_masses_scale;
   double drbar_susy_pars_scale;
   double drbar_soft_pars_scale;
   double drbar_mixings_scale;
   unsigned width;
};

template <class T>
void CNE6SSM_pole_mass_writer::extract_pole_masses(const CNE6SSM<T>& model)
{
   pole_masses.clear();
   pole_masses_scale = model.get_scale();
   pole_masses_inputs = model.get_input();
   pole_masses_problems = model.get_problems();

   pole_masses.push_back(TPoleMass("MGlu", to_valarray(PHYSICAL(MGlu))));
   pole_masses.push_back(TPoleMass("MChaP", to_valarray(PHYSICAL(MChaP))));
   pole_masses.push_back(TPoleMass("MVZp", to_valarray(PHYSICAL(MVZp))));
   pole_masses.push_back(TPoleMass("MSd", to_valarray(PHYSICAL(MSd))));
   pole_masses.push_back(TPoleMass("MSv", to_valarray(PHYSICAL(MSv))));
   pole_masses.push_back(TPoleMass("MSu", to_valarray(PHYSICAL(MSu))));
   pole_masses.push_back(TPoleMass("MSe", to_valarray(PHYSICAL(MSe))));
   pole_masses.push_back(TPoleMass("MSDX", to_valarray(PHYSICAL(MSDX))));
   pole_masses.push_back(TPoleMass("Mhh", to_valarray(PHYSICAL(Mhh))));
   pole_masses.push_back(TPoleMass("MAh", to_valarray(PHYSICAL(MAh))));
   pole_masses.push_back(TPoleMass("MHpm", to_valarray(PHYSICAL(MHpm))));
   pole_masses.push_back(TPoleMass("MChi", to_valarray(PHYSICAL(MChi))));
   pole_masses.push_back(TPoleMass("MCha", to_valarray(PHYSICAL(MCha))));
   pole_masses.push_back(TPoleMass("MFDX", to_valarray(PHYSICAL(MFDX))));
   pole_masses.push_back(TPoleMass("MSHI0", to_valarray(PHYSICAL(MSHI0))));
   pole_masses.push_back(TPoleMass("MSHIPM", to_valarray(PHYSICAL(MSHIPM))));
   pole_masses.push_back(TPoleMass("MChaI", to_valarray(PHYSICAL(MChaI))));
   pole_masses.push_back(TPoleMass("MChiI", to_valarray(PHYSICAL(MChiI))));
   pole_masses.push_back(TPoleMass("MSHp0", to_valarray(PHYSICAL(MSHp0))));
   pole_masses.push_back(TPoleMass("MSHpp", to_valarray(PHYSICAL(MSHpp))));
   pole_masses.push_back(TPoleMass("MChiP", to_valarray(PHYSICAL(MChiP))));

   if (model.do_calculate_sm_pole_masses()) {
      pole_masses.push_back(TPoleMass("MFd", to_valarray(PHYSICAL(MFd))));
      pole_masses.push_back(TPoleMass("MFe", to_valarray(PHYSICAL(MFe))));
      pole_masses.push_back(TPoleMass("MFu", to_valarray(PHYSICAL(MFu))));
      pole_masses.push_back(TPoleMass("MFv", to_valarray(PHYSICAL(MFv))));
      pole_masses.push_back(TPoleMass("MVG", to_valarray(PHYSICAL(MVG))));
      pole_masses.push_back(TPoleMass("MVP", to_valarray(PHYSICAL(MVP))));
      pole_masses.push_back(TPoleMass("MVWm", to_valarray(PHYSICAL(MVWm))));
      pole_masses.push_back(TPoleMass("MVZ", to_valarray(PHYSICAL(MVZ))));

   }
}

template <class T>
void CNE6SSM_drbar_values_writer::extract_drbar_masses(const CNE6SSM<T>& model)
{
   drbar_masses.clear();
   drbar_masses_scale = model.get_scale();
   drbar_masses_inputs = model.get_input();
   drbar_masses_problems = model.get_problems();

   drbar_masses.push_back(TMass("MVG", to_valarray(MODELPARAMETER(MVG))));
   drbar_masses.push_back(TMass("MGlu", to_valarray(MODELPARAMETER(MGlu))));
   drbar_masses.push_back(TMass("MFv", to_valarray(MODELPARAMETER(MFv))));
   drbar_masses.push_back(TMass("MChaP", to_valarray(MODELPARAMETER(MChaP))));
   drbar_masses.push_back(TMass("MVP", to_valarray(MODELPARAMETER(MVP))));
   drbar_masses.push_back(TMass("MVZ", to_valarray(MODELPARAMETER(MVZ))));
   drbar_masses.push_back(TMass("MVZp", to_valarray(MODELPARAMETER(MVZp))));
   drbar_masses.push_back(TMass("MSd", to_valarray(MODELPARAMETER(MSd))));
   drbar_masses.push_back(TMass("MSv", to_valarray(MODELPARAMETER(MSv))));
   drbar_masses.push_back(TMass("MSu", to_valarray(MODELPARAMETER(MSu))));
   drbar_masses.push_back(TMass("MSe", to_valarray(MODELPARAMETER(MSe))));
   drbar_masses.push_back(TMass("MSDX", to_valarray(MODELPARAMETER(MSDX))));
   drbar_masses.push_back(TMass("Mhh", to_valarray(MODELPARAMETER(Mhh))));
   drbar_masses.push_back(TMass("MAh", to_valarray(MODELPARAMETER(MAh))));
   drbar_masses.push_back(TMass("MHpm", to_valarray(MODELPARAMETER(MHpm))));
   drbar_masses.push_back(TMass("MChi", to_valarray(MODELPARAMETER(MChi))));
   drbar_masses.push_back(TMass("MCha", to_valarray(MODELPARAMETER(MCha))));
   drbar_masses.push_back(TMass("MFe", to_valarray(MODELPARAMETER(MFe))));
   drbar_masses.push_back(TMass("MFd", to_valarray(MODELPARAMETER(MFd))));
   drbar_masses.push_back(TMass("MFu", to_valarray(MODELPARAMETER(MFu))));
   drbar_masses.push_back(TMass("MFDX", to_valarray(MODELPARAMETER(MFDX))));
   drbar_masses.push_back(TMass("MSHI0", to_valarray(MODELPARAMETER(MSHI0))));
   drbar_masses.push_back(TMass("MSHIPM", to_valarray(MODELPARAMETER(MSHIPM))));
   drbar_masses.push_back(TMass("MChaI", to_valarray(MODELPARAMETER(MChaI))));
   drbar_masses.push_back(TMass("MChiI", to_valarray(MODELPARAMETER(MChiI))));
   drbar_masses.push_back(TMass("MSHp0", to_valarray(MODELPARAMETER(MSHp0))));
   drbar_masses.push_back(TMass("MSHpp", to_valarray(MODELPARAMETER(MSHpp))));
   drbar_masses.push_back(TMass("MChiP", to_valarray(MODELPARAMETER(MChiP))));
   drbar_masses.push_back(TMass("MVWm", to_valarray(MODELPARAMETER(MVWm))));
}

template <class T>
void CNE6SSM_drbar_values_writer::extract_drbar_susy_pars(const CNE6SSM<T>& model)
{
   drbar_susy_pars.clear();
   drbar_susy_pars_scale = model.get_scale();
   drbar_susy_pars_inputs = model.get_input();
   drbar_susy_pars_problems = model.get_problems();

   drbar_susy_pars.push_back(TParameter("Yd", 0, 3, 3, to_valarray(MODELPARAMETER(Yd))));
   drbar_susy_pars.push_back(TParameter("hE", 0, 3, 2, to_valarray(MODELPARAMETER(hE))));
   drbar_susy_pars.push_back(TParameter("Ye", 0, 3, 3, to_valarray(MODELPARAMETER(Ye))));
   drbar_susy_pars.push_back(TParameter("SigmaL", 0, 1, 1, to_valarray(MODELPARAMETER(SigmaL))));
   drbar_susy_pars.push_back(TParameter("KappaPr", 0, 1, 1, to_valarray(MODELPARAMETER(KappaPr))));
   drbar_susy_pars.push_back(TParameter("Sigmax", 0, 1, 1, to_valarray(MODELPARAMETER(Sigmax))));
   drbar_susy_pars.push_back(TParameter("gD", 0, 3, 3, to_valarray(MODELPARAMETER(gD))));
   drbar_susy_pars.push_back(TParameter("Kappa", 0, 3, 3, to_valarray(MODELPARAMETER(Kappa))));
   drbar_susy_pars.push_back(TParameter("Lambda12", 0, 2, 2, to_valarray(MODELPARAMETER(Lambda12))));
   drbar_susy_pars.push_back(TParameter("Lambdax", 0, 1, 1, to_valarray(MODELPARAMETER(Lambdax))));
   drbar_susy_pars.push_back(TParameter("fu", 0, 3, 2, to_valarray(MODELPARAMETER(fu))));
   drbar_susy_pars.push_back(TParameter("fd", 0, 3, 2, to_valarray(MODELPARAMETER(fd))));
   drbar_susy_pars.push_back(TParameter("Yu", 0, 3, 3, to_valarray(MODELPARAMETER(Yu))));
   drbar_susy_pars.push_back(TParameter("MuPr", 1, 1, 1, to_valarray(MODELPARAMETER(MuPr))));
   drbar_susy_pars.push_back(TParameter("MuPhi", 1, 1, 1, to_valarray(MODELPARAMETER(MuPhi))));
   drbar_susy_pars.push_back(TParameter("XiF", 2, 1, 1, to_valarray(MODELPARAMETER(XiF))));
   drbar_susy_pars.push_back(TParameter("g1", 0, 1, 1, to_valarray(MODELPARAMETER(g1))));
   drbar_susy_pars.push_back(TParameter("g2", 0, 1, 1, to_valarray(MODELPARAMETER(g2))));
   drbar_susy_pars.push_back(TParameter("g3", 0, 1, 1, to_valarray(MODELPARAMETER(g3))));
   drbar_susy_pars.push_back(TParameter("g1p", 0, 1, 1, to_valarray(MODELPARAMETER(g1p))));
   drbar_susy_pars.push_back(TParameter("vd", 1, 1, 1, to_valarray(MODELPARAMETER(vd))));
   drbar_susy_pars.push_back(TParameter("vu", 1, 1, 1, to_valarray(MODELPARAMETER(vu))));
   drbar_susy_pars.push_back(TParameter("vs", 1, 1, 1, to_valarray(MODELPARAMETER(vs))));
   drbar_susy_pars.push_back(TParameter("vsb", 1, 1, 1, to_valarray(MODELPARAMETER(vsb))));
   drbar_susy_pars.push_back(TParameter("vphi", 1, 1, 1, to_valarray(MODELPARAMETER(vphi))));
}

template <class T>
void CNE6SSM_drbar_values_writer::extract_drbar_soft_pars(const CNE6SSM<T>& model)
{
   drbar_soft_pars.clear();
   drbar_soft_pars_scale = model.get_scale();
   drbar_soft_pars_inputs = model.get_input();
   drbar_soft_pars_problems = model.get_problems();

   drbar_soft_pars.push_back(TParameter("TYd", 1, 3, 3, to_valarray(MODELPARAMETER(TYd))));
   drbar_soft_pars.push_back(TParameter("ThE", 1, 3, 2, to_valarray(MODELPARAMETER(ThE))));
   drbar_soft_pars.push_back(TParameter("TYe", 1, 3, 3, to_valarray(MODELPARAMETER(TYe))));
   drbar_soft_pars.push_back(TParameter("TSigmaL", 1, 1, 1, to_valarray(MODELPARAMETER(TSigmaL))));
   drbar_soft_pars.push_back(TParameter("TKappaPr", 1, 1, 1, to_valarray(MODELPARAMETER(TKappaPr))));
   drbar_soft_pars.push_back(TParameter("TSigmax", 1, 1, 1, to_valarray(MODELPARAMETER(TSigmax))));
   drbar_soft_pars.push_back(TParameter("TgD", 1, 3, 3, to_valarray(MODELPARAMETER(TgD))));
   drbar_soft_pars.push_back(TParameter("TKappa", 1, 3, 3, to_valarray(MODELPARAMETER(TKappa))));
   drbar_soft_pars.push_back(TParameter("TLambda12", 1, 2, 2, to_valarray(MODELPARAMETER(TLambda12))));
   drbar_soft_pars.push_back(TParameter("TLambdax", 1, 1, 1, to_valarray(MODELPARAMETER(TLambdax))));
   drbar_soft_pars.push_back(TParameter("Tfu", 1, 3, 2, to_valarray(MODELPARAMETER(Tfu))));
   drbar_soft_pars.push_back(TParameter("Tfd", 1, 3, 2, to_valarray(MODELPARAMETER(Tfd))));
   drbar_soft_pars.push_back(TParameter("TYu", 1, 3, 3, to_valarray(MODELPARAMETER(TYu))));
   drbar_soft_pars.push_back(TParameter("BMuPr", 2, 1, 1, to_valarray(MODELPARAMETER(BMuPr))));
   drbar_soft_pars.push_back(TParameter("BMuPhi", 2, 1, 1, to_valarray(MODELPARAMETER(BMuPhi))));
   drbar_soft_pars.push_back(TParameter("LXiF", 3, 1, 1, to_valarray(MODELPARAMETER(LXiF))));
   drbar_soft_pars.push_back(TParameter("mq2", 2, 3, 3, to_valarray(MODELPARAMETER(mq2))));
   drbar_soft_pars.push_back(TParameter("ml2", 2, 3, 3, to_valarray(MODELPARAMETER(ml2))));
   drbar_soft_pars.push_back(TParameter("mHd2", 2, 1, 1, to_valarray(MODELPARAMETER(mHd2))));
   drbar_soft_pars.push_back(TParameter("mHu2", 2, 1, 1, to_valarray(MODELPARAMETER(mHu2))));
   drbar_soft_pars.push_back(TParameter("md2", 2, 3, 3, to_valarray(MODELPARAMETER(md2))));
   drbar_soft_pars.push_back(TParameter("mu2", 2, 3, 3, to_valarray(MODELPARAMETER(mu2))));
   drbar_soft_pars.push_back(TParameter("me2", 2, 3, 3, to_valarray(MODELPARAMETER(me2))));
   drbar_soft_pars.push_back(TParameter("ms2", 2, 1, 1, to_valarray(MODELPARAMETER(ms2))));
   drbar_soft_pars.push_back(TParameter("msbar2", 2, 1, 1, to_valarray(MODELPARAMETER(msbar2))));
   drbar_soft_pars.push_back(TParameter("mH1I2", 2, 2, 2, to_valarray(MODELPARAMETER(mH1I2))));
   drbar_soft_pars.push_back(TParameter("mH2I2", 2, 2, 2, to_valarray(MODELPARAMETER(mH2I2))));
   drbar_soft_pars.push_back(TParameter("mSI2", 2, 3, 3, to_valarray(MODELPARAMETER(mSI2))));
   drbar_soft_pars.push_back(TParameter("mDx2", 2, 3, 3, to_valarray(MODELPARAMETER(mDx2))));
   drbar_soft_pars.push_back(TParameter("mDxbar2", 2, 3, 3, to_valarray(MODELPARAMETER(mDxbar2))));
   drbar_soft_pars.push_back(TParameter("mHp2", 2, 1, 1, to_valarray(MODELPARAMETER(mHp2))));
   drbar_soft_pars.push_back(TParameter("mHpbar2", 2, 1, 1, to_valarray(MODELPARAMETER(mHpbar2))));
   drbar_soft_pars.push_back(TParameter("mphi2", 2, 1, 1, to_valarray(MODELPARAMETER(mphi2))));
   drbar_soft_pars.push_back(TParameter("MassB", 1, 1, 1, to_valarray(MODELPARAMETER(MassB))));
   drbar_soft_pars.push_back(TParameter("MassWB", 1, 1, 1, to_valarray(MODELPARAMETER(MassWB))));
   drbar_soft_pars.push_back(TParameter("MassG", 1, 1, 1, to_valarray(MODELPARAMETER(MassG))));
   drbar_soft_pars.push_back(TParameter("MassBp", 1, 1, 1, to_valarray(MODELPARAMETER(MassBp))));
}

template <class T>
void CNE6SSM_drbar_values_writer::extract_drbar_mixings(const CNE6SSM<T>& model)
{
   drbar_mixings.clear();
   drbar_mixings_scale = model.get_scale();
   drbar_mixings_inputs = model.get_input();
   drbar_mixings_problems = model.get_problems();

   drbar_mixings.push_back(TMixing("ZD", 6, true, to_valarray(MODELPARAMETER(ZD))));
   drbar_mixings.push_back(TMixing("ZV", 3, true, to_valarray(MODELPARAMETER(ZV))));
   drbar_mixings.push_back(TMixing("ZU", 6, true, to_valarray(MODELPARAMETER(ZU))));
   drbar_mixings.push_back(TMixing("ZE", 6, true, to_valarray(MODELPARAMETER(ZE))));
   drbar_mixings.push_back(TMixing("ZDX", 6, true, to_valarray(MODELPARAMETER(ZDX))));
   drbar_mixings.push_back(TMixing("ZH", 5, true, to_valarray(MODELPARAMETER(ZH))));
   drbar_mixings.push_back(TMixing("ZA", 5, true, to_valarray(MODELPARAMETER(ZA))));
   drbar_mixings.push_back(TMixing("ZP", 2, true, to_valarray(MODELPARAMETER(ZP))));
   drbar_mixings.push_back(TMixing("ZN", 8, false, to_valarray(MODELPARAMETER(ZN))));
   drbar_mixings.push_back(TMixing("UM", 2, false, to_valarray(MODELPARAMETER(UM))));
   drbar_mixings.push_back(TMixing("UP", 2, false, to_valarray(MODELPARAMETER(UP))));
   drbar_mixings.push_back(TMixing("ZEL", 3, false, to_valarray(MODELPARAMETER(ZEL))));
   drbar_mixings.push_back(TMixing("ZER", 3, false, to_valarray(MODELPARAMETER(ZER))));
   drbar_mixings.push_back(TMixing("ZDL", 3, false, to_valarray(MODELPARAMETER(ZDL))));
   drbar_mixings.push_back(TMixing("ZDR", 3, false, to_valarray(MODELPARAMETER(ZDR))));
   drbar_mixings.push_back(TMixing("ZUL", 3, false, to_valarray(MODELPARAMETER(ZUL))));
   drbar_mixings.push_back(TMixing("ZUR", 3, false, to_valarray(MODELPARAMETER(ZUR))));
   drbar_mixings.push_back(TMixing("ZDXL", 3, false, to_valarray(MODELPARAMETER(ZDXL))));
   drbar_mixings.push_back(TMixing("ZDXR", 3, false, to_valarray(MODELPARAMETER(ZDXR))));
   drbar_mixings.push_back(TMixing("UHI0", 7, true, to_valarray(MODELPARAMETER(UHI0))));
   drbar_mixings.push_back(TMixing("UHIPM", 4, true, to_valarray(MODELPARAMETER(UHIPM))));
   drbar_mixings.push_back(TMixing("ZMI", 2, false, to_valarray(MODELPARAMETER(ZMI))));
   drbar_mixings.push_back(TMixing("ZPI", 2, false, to_valarray(MODELPARAMETER(ZPI))));
   drbar_mixings.push_back(TMixing("ZNI", 7, false, to_valarray(MODELPARAMETER(ZNI))));
   drbar_mixings.push_back(TMixing("UHp0", 2, true, to_valarray(MODELPARAMETER(UHp0))));
   drbar_mixings.push_back(TMixing("UHpp", 2, true, to_valarray(MODELPARAMETER(UHpp))));
   drbar_mixings.push_back(TMixing("ZNp", 2, false, to_valarray(MODELPARAMETER(ZNp))));
}

} // namespace flexiblesusy

#endif
