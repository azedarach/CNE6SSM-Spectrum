// ====================================================================
// Test suite for implementation of grid scanner class
// ====================================================================

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_fixed_point_iterator

#include <iostream>
#include <vector>
#include <boost/test/unit_test.hpp>

#include "grid_scanner.hpp"

using namespace flexiblesusy;

template <typename T>
bool check_equal_vectors(const std::vector<T>& v, 
                         const std::vector<T>& u)
{
   if (v.size() != u.size()) {
      return false;
   } else {
      for (std::size_t i = 0, n = v.size(); i < n; ++i) {
         if (v[i] != u[i]) return false;
      }
   }

   return true;
}

template <typename T>
void print_vec(std::ostream & out, const std::vector<T>& v)
{
   out << "(";
   for (std::size_t i = 0; i < v.size(); ++i) {
      out << v[i];
      if (i == v.size() - 1) {
         out << ")\n";
      } else {
         out << ", ";
      }
   }
}

BOOST_AUTO_TEST_CASE( test_step_forward )
{
   std::vector<std::size_t> dims = {3, 3, 3};

   Grid_scanner scan(dims);

   std::vector<std::size_t> posn = {1, 2, 0};
   std::vector<std::size_t> expect(posn);
   ++expect.at(2);

   scan.set_position(posn);

   scan.step_forward();

   bool result = check_equal_vectors(scan.get_position(), expect);

   BOOST_CHECK(result == true);

   posn.at(1) = 0;
   posn.at(2) = 2;
   expect.at(1) = 1;
   expect.at(2) = 0;

   scan.set_position(posn);

   scan.step_forward();

   result = check_equal_vectors(scan.get_position(), expect);

   BOOST_CHECK(result == true);

   posn.at(0) = 2;
   posn.at(1) = 2;
   posn.at(2) = 2;

   expect.at(0) = 0;
   expect.at(1) = 0;
   expect.at(2) = 0;

   scan.set_position(posn);

   scan.step_forward();

   result = check_equal_vectors(scan.get_position(), expect);

   BOOST_CHECK(result == true);
   BOOST_CHECK(scan.has_finished() == true);
}

BOOST_AUTO_TEST_CASE( test_step_backward )
{
  std::vector<std::size_t> dims = {3, 3, 3};

   Grid_scanner scan(dims);

   std::vector<std::size_t> posn = {1, 2, 0};
   std::vector<std::size_t> expect(posn);
   expect.at(1) = 1;
   expect.at(2) = 2;

   scan.set_position(posn);

   scan.step_backward();

   bool result = check_equal_vectors(scan.get_position(), expect);

   BOOST_CHECK(result == true);

   posn.at(2) = 1;
   expect = posn;
   --expect.at(2);

   scan.set_position(posn);

   scan.step_backward();

   result = check_equal_vectors(scan.get_position(), expect);

   BOOST_CHECK(result == true);

   posn.at(0) = 1;
   posn.at(1) = 0;
   posn.at(2) = 0;

   expect.at(0) = 0;
   expect.at(1) = 2;
   expect.at(2) = 2;

   scan.set_position(posn);

   scan.step_backward();

   result = check_equal_vectors(scan.get_position(), expect);

   BOOST_CHECK(result == true);

   posn.at(0) = 0;
   posn.at(1) = 0;
   posn.at(2) = 0;

   expect.at(0) = 2;
   expect.at(1) = 2;
   expect.at(2) = 2;

   scan.set_position(posn);

   scan.step_backward();

   result = check_equal_vectors(scan.get_position(), expect);

   BOOST_CHECK(result == true);
   BOOST_CHECK(scan.has_finished() == true);
}

BOOST_AUTO_TEST_CASE( test_num_points )
{
   std::vector<std::size_t> dims = {12, 2, 92, 42};

   Grid_scanner scan(dims);

   std::size_t expect = 92736;

   BOOST_CHECK_EQUAL(expect, scan.get_num_points());
}

BOOST_AUTO_TEST_CASE( test_index )
{
   std::vector<std::size_t> dims = {4, 2, 6, 3};

   Grid_scanner scan(dims);

   std::vector<std::size_t> posn = {3, 1, 1, 2};

   std::size_t expect = 131;

   scan.set_position(posn);

   BOOST_CHECK_EQUAL(expect, scan.get_current_index());
}

BOOST_AUTO_TEST_CASE( test_1d_index )
{
   std::vector<std::size_t> dims = {3};

   Grid_scanner scan(dims);

   std::vector<std::size_t> posn = {2};

   std::size_t expect = 2;

   scan.set_position(posn);

   BOOST_CHECK_EQUAL(expect, scan.get_current_index());
}

BOOST_AUTO_TEST_CASE( test_number_of_steps )
{
   int count = 0;

   std::vector<std::size_t> dims = {30, 30, 5, 23};

   Grid_scanner scan(dims);

   while (!scan.has_finished()) {
      ++count;
      scan.step_forward();
   }

   int expected = 103500;

   BOOST_CHECK_EQUAL(count, expected);
}

BOOST_AUTO_TEST_CASE( test_scan )
{
   struct MyExample {
      static const std::vector<int> func(std::vector<std::size_t> grid_pt)
         {
            const int x_incr = 1;
            const int y_incr = -2;
            const int z_incr = 3;

            std::vector<int> temp(grid_pt.size(), 0);
            
            temp.at(0) = 2 + x_incr * grid_pt.at(0);
            temp.at(1) = 1 + y_incr * grid_pt.at(1);
            temp.at(2) = z_incr * grid_pt.at(2);
            
            return temp; 
         }
   };
   
   // 3 points in each direction of a 3D grid
   std::vector<std::size_t> dims = {3, 2, 5};
   Grid_scanner scan(dims);
   std::vector<int> posn;
   while (!scan.has_finished()) {
      if (scan.get_current_index() == 23) {
         posn = MyExample::func(scan.get_position());
      }
      scan.step_forward();
   }

   std::vector<int> expect = {4, 1, 9};

   bool result = check_equal_vectors(posn, expect);

   BOOST_CHECK(result == true);

}
