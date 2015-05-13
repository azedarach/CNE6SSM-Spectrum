// ====================================================================
// Test implementation of an interface class allowing for reading
// and writing of input parameters from SLHA file.
// ====================================================================

#ifndef CNE6SSM_INPUT_PARAMETERS_SLHA_H
#define CNE6SSM_INPUT_PARAMETERS_SLHA_H

#include <Eigen/Core>

#include <string>

namespace flexiblesusy {

struct SLHA_input_entry {
   std::string block;
   std::string name;
   unsigned key;

   SLHA_input_entry(const std::string& block_, const std::string& name_,
                    unsigned key_)
      : block(block_), name(name_), key(key_)
      {}
};

struct SLHA_input_block {
   std::string block;
   std::string name;

   SLHA_input_block(const std::string& block_, const std::string& name_)
      : block(block_), name(name_)
      {}
};

struct CNE6SSM_input_parameters_slha {
   // special blocks: MINPAR and EXTPAR
   std::vector<SLHA_input_entry> minpar;
   std::vector<SLHA_input_entry> extpar;

   // additional input entries
   std::vector<SLHA_input_entry> entries;

   // additional input blocks
   std::vector<SLHA_input_block> blocks;
};

} // namespace flexiblesusy

#endif
