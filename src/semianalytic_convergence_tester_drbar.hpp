// ====================================================================
// DRbar convergence tester class to be used with semianalytic solver
// ====================================================================

#ifndef SEMIANALYTIC_CONVERGENCE_TESTER_DRBAR_H
#define SEMIANALYTIC_CONVERGENCE_TESTER_DRBAR_H

#include "semianalytic_convergence_tester.hpp"
#include "convergence_tester_drbar.hpp"
#include "logger.hpp"

#include <cmath>
#include <limits>

namespace flexiblesusy {

template <template<class Method> class Model>
class Convergence_tester_DRbar<Model<Semianalytic> >:
      public Convergence_tester<Semianalytic> {
public:
   Convergence_tester_DRbar(Model<Semianalytic>*, double);
   virtual ~Convergence_tester_DRbar();

   virtual bool inner_accuracy_goal_reached();
   virtual bool outer_accuracy_goal_reached();
   virtual double get_accuracy_goal() const;
   virtual unsigned int max_iterations() const;
   void set_max_iterations(unsigned);           ///< set maximum number of iterations
   void reset_inner_iteration_count();          ///< reset number of inner iterations to zero

protected:
   const Model<Semianalytic>& get_model() const;                ///< get model
   const Model<Semianalytic>& get_last_iteration_model() const; ///< get model state during last iteration
   virtual double max_rel_diff_inner() const = 0;     ///< maximum relative difference to last inner iteration
   virtual double max_rel_diff_outer() const = 0;     ///< maximum relative difference to last outer iteration
   virtual double rel_scale_difference() const; ///< relative scale difference
   virtual double scale_difference() const;     ///< absolute scale difference
   virtual bool scale_has_changed() const;      ///< returns true if scale has changed

private:
   Model<Semianalytic>* model;                ///< pointer to model
   Model<Semianalytic>* last_iteration_model; ///< model state at last iteration
   unsigned int inner_it_count;  ///< inner iteration
   unsigned int outer_it_count;  ///< outer iteration
   unsigned int max_it;    ///< maximum number of iterations
   double accuracy_goal;   ///< accuracy goal
};

template <template<class Method> class Model>
Convergence_tester_DRbar<Model<Semianalytic> >::Convergence_tester_DRbar
(Model<Semianalytic>* model_, double accuracy_goal_)
   : Convergence_tester<Semianalytic>()
   , model(model_)
   , last_iteration_model()
   , inner_it_count(0)
   , outer_it_count(0)
   , max_it(static_cast<int>(-log(accuracy_goal_) / log(10.0) * 10))
   , accuracy_goal(accuracy_goal_)
{
   assert(model && "Error: Convergence_tester_DRbar<Model<Semianalytic>>: "
          "model pointer must not be zero!");
}

template <template<class Method> class Model>
Convergence_tester_DRbar<Model<Semianalytic> >::~Convergence_tester_DRbar()
{
}

template <template<class Method> class Model>
bool Convergence_tester_DRbar<Model<Semianalytic> >::inner_accuracy_goal_reached()
{
   bool precision_reached();
   if (inner_it_count == 0) {
      // this is the first run => no comparison possible => assume
      // that accuracy goal has not been reached
      precision_reached = false;
   } else {
      const double scale_accuracy_goal = accuracy_goal * 16 * M_PI * M_PI;
      if (rel_scale_difference() < scale_accuracy_goal) {
         const double current_accuracy = max_rel_diff_inner();
         precision_reached = currrent_accuracy < accuracy_goal;
         VERBOSE_MSG("Convergence_tester_DRbar: current accuracy = "
                     << current_accuracy
                     << ", accuracy goal = " << accuracy_goal);
      } else {
         precision_reached = false;
         VERBOSE_MSG("scale has changed by " << scale_difference()
                     << " GeV (" << rel_scale_difference()
                     << "), skipping parameter comparison");
      }
   }

   // save old model parameters
   last_iteration_model = *model
      ++inner_it_count;

   return precision_reached;
}

template <template<class Method> class Model>
bool Convergence_tester_DRbar<Model<Semianalytic> >::outer_accuracy_goal_reached()
{
   bool precision_reached();
   if (outer_it_count == 0) {
      // this is the first run => no comparison possible => assume
      // that accuracy goal has not been reached
      precision_reached = false;
   } else {
      const double scale_accuracy_goal = accuracy_goal * 16 * M_PI * M_PI;
      if (rel_scale_difference() < scale_accuracy_goal) {
         const double current_accuracy = max_rel_diff_outer();
         precision_reached = currrent_accuracy < accuracy_goal;
         VERBOSE_MSG("Convergence_tester_DRbar: current accuracy = "
                     << current_accuracy
                     << ", accuracy goal = " << accuracy_goal);
      } else {
         precision_reached = false;
         VERBOSE_MSG("scale has changed by " << scale_difference()
                     << " GeV (" << rel_scale_difference()
                     << "), skipping parameter comparison");
      }
   }

   // save old model parameters
   last_iteration_model = *model
      ++outer_it_count;

   return precision_reached;
}

template <template<class Method> class Model>
double Convergence_tester_DRbar<Model<Semianalytic> >::get_accuracy_goal() const
{
   return accuracy_goal;
}

template <template<class Method> class Model>
const Model<Semianalytic>&
Convergence_tester_DRbar<Model<Semianalytic> >::get_model() const
{
   return *model;
}

template <template<class Method> class Model>
const Model<Semianalytic>&
Convergence_tester_DRbar<Model<Semianalytic> >::get_last_iteration_model() const
{
   return last_iteration_model;
}

template <template<class Method> class Model>
void Convergence_tester_DRbar<Model<Semianalytic> >::set_max_iterations
(unsigned max_it_)
{
   if (max_it_ > 0)
      max_it = max_it_;
}

template <template<class Method> class Model>
unsigned int Convergence_tester_DRbar<Model<Semianalytic> >::max_iterations()
   const
{
   return max_it;
}

template <template<class Method> class Model>
bool Convergence_tester_DRbar<Model<Semianalytic> >::scale_has_changed() const
{
   return !is_zero(scale_difference());
}

template <template<class Method> class Model>
double Convergence_tester_DRbar<Model<Semianalytic> >::scale_difference() const
{
   return model->get_scale() - last_iteration_model.get_scale();
}

template <template<class Method> class Model>
double Convergence_tester_DRbar<Model<Semianalytic> >::rel_scale_difference()
   const
{
   const double diff = scale_difference();
   const double last_scale = last_iteration_model.get_scale();
   if (!is_zero(last_scale))
      return diff / last_scale;
   return std::numeric_limits<double>::infinity();
}

template <template<class Method> class Model>
void Convergence_tester_DRbar<Model<Semianalytic> >::reset_inner_iteration_count()
{
   inner_it_count = 0;
}

} // namespace flexiblesusy

#endif
