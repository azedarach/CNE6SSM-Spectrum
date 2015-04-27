// ====================================================================
// Model class to be used with semianalytic solver
// ====================================================================

#ifndef SEMIANALYTIC_MODEL_H
#define SEMIANALYTIC_MODEL_H

#include "model.hpp"

#include <string>
#include <ostream>

namespace flexiblesusy {

class Semianalytic_model : public Model {
public:
   virtual ~Semianalytic_model() {}
   virtual void clear_problems() {}
   virtual std::string name() const { return "unnamed"; }
   virtual void print(std::ostream& out) const { out << "Model: " << name(); }
   friend std::ostream& operator<<(std::ostream& out, const Semianalytic_model& model) {
      model.print(out);
      return out;
   }
   virtual void set_precision(double) = 0;
};

}

#endif
