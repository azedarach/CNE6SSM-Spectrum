#include "lowNE6SSM_info.hpp"

const char * get_parameter_name_from_index(int i)
{
   return flexiblesusy::lowNE6SSM_info::parameter_names[i];
}

#include <boost/python.hpp>

BOOST_PYTHON_MODULE(_lowne6ssm)
{
   using namespace boost::python;
   def("get_parameter_name_from_index",
       get_parameter_name_from_index);
}
