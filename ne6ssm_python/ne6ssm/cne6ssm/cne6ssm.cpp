#include "CNE6SSM_info.hpp"

#include <boost/python.hpp>

using namespace boost::python;

static list get_parameter_names()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARAMETERS;
        ++i) {
      t.append(flexiblesusy::CNE6SSM_info::parameter_names[i]);
   }
   return t;
}

static list get_parameter_latex_names()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARAMETERS;
        ++i) {
      t.append(flexiblesusy::CNE6SSM_info::parameter_latex_names[i]);
   }
   return t;
}

static list get_particle_names()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARTICLES; ++i) {
      t.append(flexiblesusy::CNE6SSM_info::particle_names[i]);
   }
   return t;
}

static list get_particle_latex_names()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARTICLES; ++i) {
      t.append(flexiblesusy::CNE6SSM_info::particle_latex_names[i]);
   }
   return t;
}

static list get_particle_multiplicities()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARTICLES; ++i) {
      t.append(flexiblesusy::CNE6SSM_info::particle_multiplicities[i]);
   }
   return t;
}

static list get_input_names()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_INPUTS; ++i) {
      t.append(flexiblesusy::CNE6SSM_info::input_names[i]);
   }
   return t;
}

static list get_input_latex_names()
{
   list t;
   for (std::size_t i = 0; i < flexiblesusy::CNE6SSM_info::NUMBER_OF_INPUTS; ++i) {
      t.append(flexiblesusy::CNE6SSM_info::input_latex_names[i]);
   }
   return t;
}

BOOST_PYTHON_MODULE(_cne6ssm)
{
   def("get_parameter_names",
       get_parameter_names);
   def("get_parameter_latex_names",
       get_parameter_latex_names);
   def("get_particle_names",
       get_particle_names);
   def("get_particle_latex_names",
       get_particle_latex_names);
   def("get_particle_multiplicities",
       get_particle_multiplicities);
   def("get_input_names",
       get_input_names);
   def("get_input_latex_names",
       get_input_latex_names);
}
