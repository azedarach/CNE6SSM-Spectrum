// ====================================================================
// Initial guesser class to be used with semianalytic solver
// ====================================================================

#ifndef SEMIANALYTIC_INITIAL_GUESSER_H
#define SEMIANALYTIC_INITIAL_GUESSER_H

#include "initial_guesser.hpp"

namespace flexiblesusy {

class Semianalytic;

template<>
class Initial_guesser<Semianalytic> {
public:
   virtual ~Initial_guesser() {}
   virtual void guess() = 0;
};

}

#endif
