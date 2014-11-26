#include "CNE6SSM_info.hpp"
#include <string>
#include <stdlib.h>
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

bool is_particle_name(const char* name)
{
   std::string input_name(name); //< much more convenient
   std::string particle_name;
   // see if index is attached
   std::size_t open_bracket_posn = input_name.find_first_of('(');
   if (open_bracket_posn != std::string::npos) {
      particle_name = input_name.substr(0, open_bracket_posn); 
   } else {
      particle_name = input_name;
   }

   // remove leading 'M', if present
   std::size_t mass_prefix_posn = particle_name.find_first_of('M');
   if (mass_prefix_posn != std::string::npos) {
      particle_name = particle_name.substr(1);
   }

   std::size_t loc = 0;
   while (loc < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARTICLES) {
      if (!strcmp(particle_name.c_str(), flexiblesusy::CNE6SSM_info::particle_names[loc]))
         return true;
      ++loc;
   }

   return false;
}

const char * convert_to_particle_latex_name(const char* name)
{
   std::string input_name(name); //< much more convenient
   std::string particle_name;
   // see if index is attached
   std::size_t open_bracket_posn = input_name.find_first_of('(');
   if (open_bracket_posn != std::string::npos) {
      particle_name = input_name.substr(0, open_bracket_posn); 
   } else {
      particle_name = input_name;
   }

   // remove leading 'M', if present
   bool is_mass = false;
   std::size_t mass_prefix_posn = particle_name.find_first_of('M');
   if (mass_prefix_posn != std::string::npos) {
      particle_name = particle_name.substr(1);
      is_mass = true;
   }

   std::size_t loc = 0;
   while (loc < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARTICLES) {
      if (!strcmp(particle_name.c_str(), flexiblesusy::CNE6SSM_info::particle_names[loc]))
         break;
      ++loc;
   }

   if (loc == flexiblesusy::CNE6SSM_info::NUMBER_OF_PARTICLES)
      return "";

   std::string latex_name(flexiblesusy::CNE6SSM_info::particle_latex_names[loc]);

   // if index is provided, locate closing bracket
   // and check that index is valid
   if (open_bracket_posn != std::string::npos) {
      std::size_t close_bracket_posn = input_name.find_first_of(')', open_bracket_posn);
      if (close_bracket_posn == std::string::npos) {
         return "";
      }
      std::string gen = input_name.substr(open_bracket_posn + 1, close_bracket_posn - open_bracket_posn - 1);
      std::size_t gen_num = atoi(gen.c_str());
      if (gen_num <= 0 || gen_num > flexiblesusy::CNE6SSM_info::particle_multiplicities[loc]) {
         return "";
      } else {
         latex_name += "_{" + gen + "}";
      }
   }

   if (is_mass) {
      latex_name = "m_{" + latex_name + "}";
   }

   return latex_name.c_str();
}

bool is_parameter_name(const char* name)
{
   std::size_t loc = 0;
   while (loc < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARAMETERS) {
      if (!strcmp(name, flexiblesusy::CNE6SSM_info::parameter_names[loc]))
         return true;
      ++loc;
   }

   return false;
}

const char* convert_to_parameter_latex_name(const char* name)
{
   std::size_t index = 0;
   while (index < flexiblesusy::CNE6SSM_info::NUMBER_OF_PARAMETERS) {
      if (!strcmp(name, flexiblesusy::CNE6SSM_info::parameter_names[index])) {
         return flexiblesusy::CNE6SSM_info::parameter_latex_names[index];
      }
      ++index;
   }
   return "";
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

bool is_input_name(const char* name)
{
   std::size_t loc = 0;
   while (loc < flexiblesusy::CNE6SSM_info::NUMBER_OF_INPUTS) {
      if (!strcmp(name, flexiblesusy::CNE6SSM_info::input_names[loc]))
         return true;
      ++loc;
   }

   return false;
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

const char* convert_to_latex_name(const char* name)
{
   if (is_particle_name(name)) {
      return convert_to_particle_latex_name(name);
   } else if (is_parameter_name(name)) {
      return convert_to_parameter_latex_name(name);
   } else if (is_input_name(name)) {
      return convert_to_input_latex_name(name);
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
   def("convert_to_particle_latex_name",
       convert_to_particle_latex_name);
   def("get_all_input_names",
       get_all_input_names);
   def("get_all_input_latex_names",
       get_all_input_latex_names);
   def("convert_to_input_latex_name",
       convert_to_input_latex_name);
   def("convert_to_latex_name",
       convert_to_latex_name);
}
