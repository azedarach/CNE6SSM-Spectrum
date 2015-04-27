// ====================================================================
// Running precision calculator class to be used with semianalytic solver
// ====================================================================

#ifndef SEMIANALYTIC_RUNNING_PRECISION_H
#define SEMIANALYTIC_RUNNING_PRECISION_H

namespace flexiblesusy {

class Semianalytic_running_precision {
public:
   virtual ~Semianalytic_running_precision() {}
   virtual double get_precision(unsigned) = 0;
};

class Semianalytic_constant_precision : public Semianalytic_running_precision {
public:
   explicit Semianalytic_constant_precision(double);
   virtual ~Semianalytic_constant_precision();
   virtual double get_precision(unsigned);
private:
   double precision;
};

class Semianalytic_increasing_precision : public Semianalytic_running_precision {
public:
   Semianalytic_increasing_precision(double, double);
   virtual ~Semianalytic_increasing_precision();
   virtual double get_precision(unsigned);
private:
   double decreasing_factor;
   double minimum_precision;
};

} // namespace flexiblesusy

#endif
