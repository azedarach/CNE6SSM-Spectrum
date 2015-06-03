
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_derivs_root_finder

#include <boost/test/unit_test.hpp>

#include "derivs_root_finder.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_parabola_2dim )
{
   struct Parabola {
      static int func(const gsl_vector* x, void*, gsl_vector* f) {
         const double y = gsl_vector_get(x, 0);
         const double z = gsl_vector_get(x, 1);
         gsl_vector_set(f, 0, y*(y - 5.0));
         gsl_vector_set(f, 1, z*(z - 1.0));
         return GSL_SUCCESS;
      }

      static int dfunc(const gsl_vector* x, void*, gsl_matrix* J) {
         const double y = gsl_vector_get(x, 0);
         const double z = gsl_vector_get(x, 1);
         gsl_matrix_set(J, 0, 0, 2*y - 5.0);
         gsl_matrix_set(J, 0, 1, 0.0);
         gsl_matrix_set(J, 1, 0, 0.0);
         gsl_matrix_set(J, 1, 1, 2*z - 1.0);
         return GSL_SUCCESS;
      }

      static int fdfunc(const gsl_vector* x, void*, gsl_vector* f, gsl_matrix* J) {
         const double y = gsl_vector_get(x, 0);
         const double z = gsl_vector_get(x, 1);
         gsl_vector_set(f, 0, y*(y - 5.0));
         gsl_vector_set(f, 1, z*(z - 1.0));
         gsl_matrix_set(J, 0, 0, 2*y - 5.0);
         gsl_matrix_set(J, 0, 1, 0.0);
         gsl_matrix_set(J, 1, 0, 0.0);
         gsl_matrix_set(J, 1, 1, 2*z - 1.0);
         return GSL_SUCCESS;
      }
   };

   Derivs_root_finder<2> root_finder(Parabola::func, Parabola::dfunc, 
                                     Parabola::fdfunc, NULL, 100, 1.0e-5);
   const double start[2] = { 10, 10 };
   const int status = root_finder.find_root(start);

   BOOST_CHECK_EQUAL(status, GSL_SUCCESS);
   BOOST_CHECK_CLOSE_FRACTION(root_finder.get_root(0), 5.0, 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(root_finder.get_root(1), 1.0, 1.0e-5);
}

struct Parabola {
   static unsigned number_of_function_calls;
   static unsigned number_of_jacobian_calls;
   static unsigned number_of_combined_calls;

   static int func(const gsl_vector* x, void*, gsl_vector* f) {
      const double y = gsl_vector_get(x, 0);
      const double z = gsl_vector_get(x, 1);
      gsl_vector_set(f, 0, y*(y - 5.0));
      gsl_vector_set(f, 1, z*(z - 1.0));
      number_of_function_calls++;
      return GSL_SUCCESS;
   }

   static int dfunc(const gsl_vector* x, void*, gsl_matrix* J) {
      const double y = gsl_vector_get(x, 0);
      const double z = gsl_vector_get(x, 1);
      gsl_matrix_set(J, 0, 0, 2*y - 5.0);
      gsl_matrix_set(J, 0, 1, 0.0);
      gsl_matrix_set(J, 1, 0, 0.0);
      gsl_matrix_set(J, 1, 1, 2*z - 1.0);
      number_of_jacobian_calls++;
      return GSL_SUCCESS;
   }

   static int fdfunc(const gsl_vector* x, void*, gsl_vector* f, gsl_matrix* J) {
      const double y = gsl_vector_get(x, 0);
      const double z = gsl_vector_get(x, 1);
      gsl_vector_set(f, 0, y*(y - 5.0));
      gsl_vector_set(f, 1, z*(z - 1.0));
      gsl_matrix_set(J, 0, 0, 2*y - 5.0);
      gsl_matrix_set(J, 0, 1, 0.0);
      gsl_matrix_set(J, 1, 0, 0.0);
      gsl_matrix_set(J, 1, 1, 2*z - 1.0);
      number_of_combined_calls++;
      return GSL_SUCCESS;
   };
};

unsigned Parabola::number_of_function_calls = 0;
unsigned Parabola::number_of_jacobian_calls = 0;
unsigned Parabola::number_of_combined_calls = 0;

BOOST_AUTO_TEST_CASE( test_number_of_calls )
{
   const double precision = 1.0e-5;
   Derivs_root_finder<2> root_finder(Parabola::func, Parabola::dfunc,
                                     Parabola::fdfunc, NULL, 100, precision);
   const double start[2] = { 10, 10 };
   int status = GSL_SUCCESS;

   const gsl_multiroot_fdfsolver_type* solvers[] =
      {gsl_multiroot_fdfsolver_hybridj, gsl_multiroot_fdfsolver_hybridsj,
       gsl_multiroot_fdfsolver_gnewton, gsl_multiroot_fdfsolver_newton};

   for (std::size_t i = 0; i < sizeof(solvers)/sizeof(*solvers); ++i) {
      Parabola::number_of_function_calls = 0;
      Parabola::number_of_jacobian_calls = 0;
      Parabola::number_of_combined_calls = 0;

      root_finder.set_solver_type(solvers[i]);
      status = root_finder.find_root(start);

      BOOST_REQUIRE(status == GSL_SUCCESS);
      BOOST_CHECK_CLOSE_FRACTION(root_finder.get_root(0), 5.0, precision);
      BOOST_CHECK_CLOSE_FRACTION(root_finder.get_root(1), 1.0, precision);
      BOOST_MESSAGE("solver type " << i << " used " << Parabola::number_of_function_calls << " function calls");
      BOOST_MESSAGE("solver type " << i << " used " << Parabola::number_of_jacobian_calls << " Jacobian calls");
      BOOST_MESSAGE("solver type " << i << " used " << Parabola::number_of_combined_calls << " combined calls");
   }
}
