// ====================================================================
// Implementation of convergence tester for CNE6SSM semianalytic solver
// ====================================================================

#include "CNE6SSM_semianalytic_convergence_tester.hpp"
#include <cmath>
#include <algorithm>
#include "wrappers.hpp"

namespace flexiblesusy {

#define OLD(p) ol.get_##p()
#define NEW(p) ne.get_##p()

#define OLD1(p,i) ol.get_##p()(i)
#define NEW1(p,i) ne.get_##p()(i)

#define OLD2(p,i,j) ol.get_##p(i,j)
#define NEW2(p,i,j) ne.get_##p(i,j)

#define OLD3(p,i,j,k) ol.get_##p(i,j,k)
#define NEW3(p,i,j,k) ne.get_##p(i,j,k)

#define OLD4(p,i,j,k,l) ol.get_##p(i,j,k,l)
#define NEW4(p,i,j,k,l) ne.get_##p(i,j,k,l)

CNE6SSM_convergence_tester<Semianalytic>::CNE6SSM_convergence_tester(CNE6SSM<Semianalytic>* model, double accuracy_goal)
   : Convergence_tester_DRbar<CNE6SSM<Semianalytic> >(model, accuracy_goal)
{
}

CNE6SSM_convergence_tester<Semianalytic>::~CNE6SSM_convergence_tester()
{
}

double CNE6SSM_convergence_tester<Semianalytic>::max_rel_diff_inner() const
{
   const CNE6SSM<Semianalytic>& ol = get_last_iteration_model();
   const CNE6SSM<Semianalytic>& ne = get_model();

   double diff[84] = { 0. };

   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i] = MaxRelDiff(OLD2(Yd,i,j),NEW2(Yd,i,j));
      }
   }
   for (unsigned j = 0; j < 2; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 9] = MaxRelDiff(OLD2(hE,i,j),NEW2(hE,i,j));
      }
   }
   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 15] = MaxRelDiff(OLD2(Ye,i,j),NEW2(Ye,i,j));
      }
   }
   diff[24] = MaxRelDiff(OLD(SigmaL),NEW(SigmaL));
   diff[25] = MaxRelDiff(OLD(KappaPr),NEW(KappaPr));
   diff[26] = MaxRelDiff(OLD(Sigmax),NEW(Sigmax));
   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 27] = MaxRelDiff(OLD2(gD,i,j),NEW2(gD,i,j));
      }
   }
   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 36] = MaxRelDiff(OLD2(Kappa,i,j),NEW2(Kappa,i,j));
      }
   }
   for (unsigned j = 0; j < 2; ++j) {
      for (unsigned i = 0; i < 2; ++i) {
         diff[2 * j + i + 45] = MaxRelDiff(OLD2(Lambda12,i,j),NEW2(Lambda12,i,j));
      }
   }
   diff[49] = MaxRelDiff(OLD(Lambdax),NEW(Lambdax));
   for (unsigned j = 0; j < 2; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 50] = MaxRelDiff(OLD2(fu,i,j),NEW2(fu,i,j));
      }
   }
   for (unsigned j = 0; j < 2; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 56] = MaxRelDiff(OLD2(fd,i,j),NEW2(fd,i,j));
      }
   }
   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 62] = MaxRelDiff(OLD2(Yu,i,j),NEW2(Yu,i,j));
      }
   }
   diff[71] = MaxRelDiff(OLD(MuPr),NEW(MuPr));
   diff[72] = MaxRelDiff(OLD(MuPhi),NEW(MuPhi));
   diff[73] = MaxRelDiff(OLD(XiF),NEW(XiF));
   diff[74] = MaxRelDiff(OLD(g1),NEW(g1));
   diff[75] = MaxRelDiff(OLD(g2),NEW(g2));
   diff[76] = MaxRelDiff(OLD(g3),NEW(g3));
   diff[77] = MaxRelDiff(OLD(g1p),NEW(g1p));
   diff[78] = MaxRelDiff(OLD(vd),NEW(vd));
   diff[79] = MaxRelDiff(OLD(vu),NEW(vu));
   diff[80] = MaxRelDiff(OLD(vs),NEW(vs));
   diff[81] = MaxRelDiff(OLD(vsb),NEW(vsb));
   diff[82] = MaxRelDiff(OLD(vphi),NEW(vphi));
   diff[83] = MaxRelDiff(OLD(QS),NEW(QS));

   return *std::max_element(diff, diff + 84);

}

double CNE6SSM_convergence_tester<Semianalytic>::max_rel_diff_outer() const
{
   const CNE6SSM<Semianalytic>& ol = get_last_iteration_model();
   const CNE6SSM<Semianalytic>& ne = get_model();

   double diff[81] = { 0 };

   diff[0] = MaxRelDiff(OLD(MGlu),NEW(MGlu));
   diff[1] = MaxRelDiff(OLD(MChaP),NEW(MChaP));
   diff[2] = MaxRelDiff(OLD(MVZp),NEW(MVZp));
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 3] = MaxRelDiff(OLD1(MSd,i),NEW1(MSd,i));
   }
   for (unsigned i = 0; i < 3; i++) {
      diff[i + 9] = MaxRelDiff(OLD1(MSv,i),NEW1(MSv,i));
   }
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 12] = MaxRelDiff(OLD1(MSu,i),NEW1(MSu,i));
   }
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 18] = MaxRelDiff(OLD1(MSe,i),NEW1(MSe,i));
   }
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 24] = MaxRelDiff(OLD1(MSDX,i),NEW1(MSDX,i));
   }
   for (unsigned i = 0; i < 5; i++) {
      diff[i + 30] = MaxRelDiff(OLD1(Mhh,i),NEW1(Mhh,i));
   }
   for (unsigned i = 2; i < 5; i++) {
      diff[i + 35] = MaxRelDiff(OLD1(MAh,i),NEW1(MAh,i));
   }
   for (unsigned i = 1; i < 2; i++) {
      diff[i + 40] = MaxRelDiff(OLD1(MHpm,i),NEW1(MHpm,i));
   }
   for (unsigned i = 0; i < 8; i++) {
      diff[i + 42] = MaxRelDiff(OLD1(MChi,i),NEW1(MChi,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 50] = MaxRelDiff(OLD1(MCha,i),NEW1(MCha,i));
   }
   for (unsigned i = 0; i < 3; i++) {
      diff[i + 52] = MaxRelDiff(OLD1(MFDX,i),NEW1(MFDX,i));
   }
   for (unsigned i = 0; i < 7; i++) {
      diff[i + 55] = MaxRelDiff(OLD1(MSHI0,i),NEW1(MSHI0,i));
   }
   for (unsigned i = 0; i < 4; i++) {
      diff[i + 62] = MaxRelDiff(OLD1(MSHIPM,i),NEW1(MSHIPM,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 66] = MaxRelDiff(OLD1(MChaI,i),NEW1(MChaI,i));
   }
   for (unsigned i = 0; i < 7; i++) {
      diff[i + 68] = MaxRelDiff(OLD1(MChiI,i),NEW1(MChiI,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 75] = MaxRelDiff(OLD1(MSHp0,i),NEW1(MSHp0,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 77] = MaxRelDiff(OLD1(MSHpp,i),NEW1(MSHpp,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 79] = MaxRelDiff(OLD1(MChiP,i),NEW1(MChiP,i));
   }

   return *std::max_element(diff, diff + 81);

}

} // namespace flexiblesusy
