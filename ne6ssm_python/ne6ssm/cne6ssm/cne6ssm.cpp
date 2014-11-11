#include "CNE6SSM_info.hpp"
#include <cstring>

#include <boost/python.hpp>

using namespace boost::python;

static list get_all_parameter_names()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARAMETERS;
        ++i) {
      t.append(flexiblesusy::CNE6SSM_info::parameter_names[i]);
   }
   return t;
}

static list get_all_parameter_latex_names()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARAMETERS;
        ++i) {
      t.append(flexiblesusy::CNE6SSM_info::parameter_latex_names[i]);
   }
   return t;
}

static list get_all_particle_names()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARTICLES; ++i) {
      t.append(flexiblesusy::CNE6SSM_info::particle_names[i]);
   }
   return t;
}

static list get_all_particle_latex_names()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARTICLES; ++i) {
      t.append(flexiblesusy::CNE6SSM_info::particle_latex_names[i]);
   }
   return t;
}

static list get_all_particle_multiplicities()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARTICLES; ++i) {
      t.append(flexiblesusy::CNE6SSM_info::particle_multiplicities[i]);
   }
   return t;
}

static list get_all_input_names()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_INPUTS; ++i) {
      t.append(flexiblesusy::CNE6SSM_info::input_names[i]);
   }
   return t;
}

static list get_all_input_latex_names()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_INPUTS; ++i) {
      t.append(flexiblesusy::CNE6SSM_info::input_latex_names[i]);
   }
   return t;
}

const char* convert_to_input_latex_name(const char* name)
{
   std::size_t index = 0;
   while (index < flexiblesusy::CNE6SSM_info::NUMBER_OF_INPUTS) {
      if (!strcmp(name, flexiblesusy::CNE6SSM_info::input_names[index])) {
         return flexiblesusy::CNE6SSM_info::input_latex_names[index];
      }
      ++index;
   }
   return "";
}

BOOST_PYTHON_MODULE(_cne6ssm)
{
   def("get_all_parameter_names",
       get_all_parameter_names);
   def("get_all_parameter_latex_names",
       get_all_parameter_latex_names);
   def("get_all_particle_names",
       get_all_particle_names);
   def("get_all_particle_latex_names",
       get_all_particle_latex_names);
   def("get_all_particle_multiplicities",
       get_all_particle_multiplicities);
   def("get_all_input_names",
       get_all_input_names);
   def("get_all_input_latex_names",
       get_all_input_latex_names);
   def("convert_to_input_latex_name",
       convert_to_input_latex_name);
}
