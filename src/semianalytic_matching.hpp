// ====================================================================
// Matching class to be used with semianalytic solver
// ====================================================================

#ifndef SEMIANALYTIC_MATCHING_H
#define SEMIANALYTIC_MATCHING_H

#include "matching.hpp"

namespace flexiblesusy {

class Semianalytic;
class Semianalytic_model;

template<>
class Matching<Semianalytic> {
public:
   virtual ~Matching() {}
   virtual void match_low_to_high_scale_model() = 0;
   virtual void match_high_to_low_scale_model() = 0;
   virtual double get_scale() const = 0;
   virtual void set_models(Semianalytic_model* low, Semianalytic_model* high) = 0;
};

}

#endif
