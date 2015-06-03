// ====================================================================
// Test implementation of high scale constraint to be used in
// semianalytic version of the two-scale algorithm
// ====================================================================

#ifndef CNE6SSM_SEMI_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H
#define CNE6SSM_SEMI_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H

#include "CNE6SSM_semi_high_scale_constraint.hpp"
#include "CNE6SSM_semi_two_scale_input_parameters.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class CNE6SSM_semianalytic;

class Two_scale;

template<>
class CNE6SSM_semianalytic_high_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   CNE6SSM_semianalytic_high_scale_constraint();
   CNE6SSM_semianalytic_high_scale_constraint(CNE6SSM_semianalytic<Two_scale>*);
   virtual ~CNE6SSM_semianalytic_high_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   const CNE6SSM_semianalytic_input_parameters<Two_scale>& get_input_parameters() const;
   CNE6SSM_semianalytic<Two_scale>* get_model() const;
   void initialize();
   void set_initial_scale_guess(double); ///< set initial scale guess before iteration
   void set_scale(double); ///< fix unification scale (0 = unfixed)

protected:
   void update_scale();

private:
   double scale;
   double initial_scale_guess;
   CNE6SSM_semianalytic<Two_scale>* model;
};

} // namespace flexiblesusy

#endif
