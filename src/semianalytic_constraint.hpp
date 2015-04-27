// ====================================================================
// Constraint class to be used with semianalytic solver
// ====================================================================

#ifndef SEMIANALYTIC_CONSTRAINT_H
#define SEMIANALYTIC_CONSTRAINT_H

#include "constraint.hpp"
#include "logger.hpp"

namespace flexiblesusy {

class Semianalytic;
class Semianalytic_model;

template<>
class Constraint<Semianalytic> {
public:
   virtual ~Constraint() {}
   virtual void apply() = 0;                        ///< apply constraint
   virtual double get_scale() const = 0;            ///< get scale where to apply
   virtual void set_model(Semianalytic_model*) = 0; ///< set model where to apply the constraint
};

} // namespace flexiblesusy

#endif
