// ====================================================================
// Test implementation of SUSY scale constraint for semianalytic
// algorithm
// ====================================================================

#ifndef CNE6SSM_SEMIANALYTIC_SUSY_SCALE_CONSTRAINT_H
#define CNE6SSM_SEMIANALYTIC_SUSY_SCALE_CONSTRAINT_H

#include "CNE6SSM_susy_scale_constraint.hpp"
#include "CNE6SSM_semianalytic_input_parameters.hpp"
#include "CNE6SSM_semianalytic_high_scale_constraint.hpp"
#include "semianalytic_constraint.hpp"

namespace flexiblesusy {

template <class T>
class CNE6SSM;

class Semianalytic;

template<>
class CNE6SSM_susy_scale_constraint<Semianalytic> : public Constraint<Semianalytic> {
public:
   CNE6SSM_susy_scale_constraint();
   CNE6SSM_susy_scale_constraint(CNE6SSM<Semianalytic>*);
   virtual ~CNE6SSM_susy_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Semianalytic_model*);

   void clear();
   double get_initial_scale_guess() const;
   const CNE6SSM_input_parameters<Semianalytic>& get_input_parameters() const;
   CNE6SSM<Semianalytic>* get_model() const;
   void initialize();

   CNE6SSM_high_scale_constraint<Semianalytic>* get_input_scale_constraint() const;
   void set_input_scale_constraint(CNE6SSM_high_scale_constraint<Semianalytic>*);

protected:
   void update_scale();

private:
   double scale;
   double initial_scale_guess;
   CNE6SSM<Semianalytic>* model;
   CNE6SSM_high_scale_constraint<Semianalytic>* high_constraint;
};

} // namespace flexiblesusy

#endif
