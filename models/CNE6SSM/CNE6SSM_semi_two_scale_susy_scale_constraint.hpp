// ====================================================================
// Test implementation of SUSY scale constraint to be used in
// semianalytic version of the two-scale algorithm
// ====================================================================

#ifndef CNE6SSM_SEMI_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H
#define CNE6SSM_SEMI_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H

#include "CNE6SSM_semi_susy_scale_constraint.hpp"
#include "CNE6SSM_semi_two_scale_input_parameters.hpp"
#include "CNE6SSM_semi_two_scale_high_scale_constraint.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class CNE6SSM_semianalytic;

class Two_scale;

template<>
class CNE6SSM_semianalytic_susy_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   CNE6SSM_semianalytic_susy_scale_constraint();
   CNE6SSM_semianalytic_susy_scale_constraint(CNE6SSM_semianalytic<Two_scale>*);
   virtual ~CNE6SSM_semianalytic_susy_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   const CNE6SSM_semianalytic_input_parameters<Two_scale>& get_input_parameters() const;
   CNE6SSM_semianalytic<Two_scale>* get_model() const;
   void initialize();

   CNE6SSM_semianalytic_high_scale_constraint<Two_scale>* get_input_scale_constraint() const;
   void set_input_scale_constraint(CNE6SSM_semianalytic_high_scale_constraint<Two_scale>*);

protected:
   void update_scale();

private:
   double scale;
   double initial_scale_guess;
   CNE6SSM_semianalytic<Two_scale>* model;
   CNE6SSM_semianalytic_high_scale_constraint<Two_scale>* high_constraint;
};

} // namespace flexiblesusy

#endif
