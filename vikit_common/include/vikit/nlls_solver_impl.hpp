/*
 * Abstract Nonlinear Least-Squares Solver Class
 *
 * nlls_solver.h
 *
 *  Created on: Nov 5, 2012
 *      Author: cforster
 */

#ifndef LM_SOLVER_IMPL_HPP_
#define LM_SOLVER_IMPL_HPP_

#include <stdexcept>

template <int D, typename T>
void vk::NLLSSolver<D, T>::
optimize(ModelType& model)
{
  if(method_ == GaussNewton)
    optimizeGaussNewton(model);
  else if(method_ == LevenbergMarquardt)
    optimizeLevenbergMarquardt(model);
}

template <int D, typename T>
void vk::NLLSSolver<D, T>::
optimizeGaussNewton(ModelType& model)
{
  // Compute weight scale
  if(use_weights_)
    computeResiduals(model, false, true);

  // Save the old model to rollback in case of unsuccessful update
  ModelType old_model(model);

  // perform iterative estimation
  for (iter_ = 0; iter_<n_iter_; ++iter_)
  {
    rho_ = 0;
    startIteration();


    H_.setZero();
    Jres_.setZero();

    // compute initial error
    n_meas_ = 0;
    double new_chi2 = computeResiduals(model, true);

    // solve the linear system
    if(!solve())
    {
      // matrix was singular and could not be computed
      cout << "Matrix is close to singular! Stop Optimizing." << endl;
      cout << "H = " << H_ << endl;
      cout << "Jres = " << Jres_ << endl;
      stop_ = true;
    }

    // check if error increased since last optimization
    if((iter_ > 0 && new_chi2 > chi2_) || stop_)
    {
      if(verbose_)
      {
        cout << "It. " << iter_
             << "\t Failure"
             << "\t new_chi2 = " << new_chi2
             << "\t Error increased. Stop optimizing."
             << endl;
      }
      model = old_model; // rollback
      break;
    }

    // update the model
    ModelType new_model;
    update(model, new_model);
    old_model = model;
    model = new_model;

    chi2_ = new_chi2;

    if(verbose_)
    {
      cout << "It. " << iter_
           << "\t Success"
           << "\t new_chi2 = " << new_chi2
           << "\t n_meas = " << n_meas_
           << endl;
    }

    finishIteration();

    // stop when converged, i.e. update step too small
    if(vk::norm_max(x_)<=eps_)
      break;
  }

  // reset
  is_initial_chi2_provided_ = false;
}

template <int D, typename T>
void vk::NLLSSolver<D, T>::
optimizeLevenbergMarquardt(ModelType& model)
{
  // Compute weight scale
  if(use_weights_)
    computeResiduals(model, false, true);

  // compute the initial error
  if(!is_initial_chi2_provided_)
    chi2_ = computeResiduals(model, true);

  if(verbose_)
    cout << "init chi2 = " << chi2_
         << "\t n_meas = " << n_meas_
         << endl;

  // TODO: compute initial lambda
  // Hartley and Zisserman: "A typical init value of lambda is 10^-3 times the
  // average of the diagonal elements of J'J"

  // Compute Initial Lambda
  if(mu_ < 0)
  {
    double H_max_diag = 0;
    double tau = 1e-4;
    for(size_t j=0; j<D; ++j)
      H_max_diag = max(H_max_diag, fabs(H_(j,j)));
    mu_ = tau*H_max_diag;
  }

  // perform iterative estimation
  for (iter_ = 0; iter_<n_iter_; ++iter_)
  {
    rho_ = 0;
    startIteration();

    // try to compute and update, if it fails, try with increased mu
    n_trials_ = 0;
    do
    {
      // init variables
      ModelType new_model;
      double new_chi2 = -1;
      H_.setZero();
      //H_ = mu_ * Matrix<double,D,D>::Identity(D,D);
      Jres_.setZero();

      // compute initial error
      n_meas_ = 0;
      computeResiduals(model, true);

      // add damping term:
      H_ += (H_.diagonal()*mu_).asDiagonal();

      // solve the linear system
      if(solve())
      {
        // update the model
        update(model, new_model);

        // compute error with new model and compare to old error
        n_meas_ = 0;
        new_chi2 = computeResiduals(new_model, false);
        rho_ = chi2_-new_chi2;
      }
      else
      {
        // matrix was singular and could not be computed
        cout << "Matrix is close to singular!" << endl;
        cout << "H = " << H_ << endl;
        cout << "Jres = " << Jres_ << endl;
        rho_ = -1;
      }

      if(rho_>0)
      {
        // update decrased the error -> success
        model = new_model;
        chi2_ = new_chi2;
        stop_ = vk::norm_max(x_)<=eps_;
        mu_ *= max(1./3., min(1.-pow(2*rho_-1,3), 2./3.));
        nu_ = 2.;
        if(verbose_)
        {
          cout << "It. " << iter_
               << "\t Trial " << n_trials_
               << "\t Success"
               << "\t n_meas = " << n_meas_
               << "\t new_chi2 = " << new_chi2
               << "\t mu = " << mu_
               << "\t nu = " << nu_
               << endl;
        }
      }
      else
      {
        // update increased the error -> fail
        mu_ *= nu_;
        nu_ *= 2.;
        ++n_trials_;
        if (n_trials_ >= n_trials_max_)
          stop_ = true;

        if(verbose_)
        {
          cout << "It. " << iter_
               << "\t Trial " << n_trials_
               << "\t Failure"
               << "\t n_meas = " << n_meas_
               << "\t new_chi2 = " << new_chi2
               << "\t mu = " << mu_
               << "\t nu = " << nu_
               << endl;
        }
      }

      finishTrial();

    } while(!(rho_>0 || stop_));
    if (stop_)
      break;

    finishIteration();
  }

  // reset
  is_initial_chi2_provided_ = false;
}


template <int D, typename T>
void vk::NLLSSolver<D, T>::
setRobustCostFunction(ScaleEstimatorType scale_estimator_t,
                      WeightFunctionType weight_function_t)
{
  switch(scale_estimator_t)
  {
    case UnitScale:
      printf("Using Unit Scale Estimator\n");
      scale_estimator_ = new robust_cost::UnitScaleEstimator();
      break;
    case TDistScale:
      printf("Using TDistribution Scale Estimator\n");
      scale_estimator_ = new robust_cost::TDistributionScaleEstimator();
      break;
    case MADScale:
      printf("Using MAD Scale Estimator\n");
      scale_estimator_ = new robust_cost::MADScaleEstimator();
    break;
    case NormalScale:
      printf("Using Normal Scale Estimator\n");
      scale_estimator_ = new robust_cost::NormalDistributionScaleEstimator();
      break;
    default: std::runtime_error("Scale Estimator unknown");
  }

  switch(weight_function_t)
  {
    case UnitWeight:
      printf("Using Unit Weight Function\n");
      weight_function_ = new robust_cost::UnitWeightFunction();
      break;
    case TDistWeight:
      printf("Using TDistribution Weight Function\n");
      weight_function_ = new robust_cost::TDistributionWeightFunction();
      break;
    case TukeyWeight:
      printf("Using Tukey Weight Function\n");
      weight_function_ = new robust_cost::TukeyWeightFunction();
      break;
    case HuberWeight:
      printf("Using Huber Weight Function\n");
      weight_function_ = new robust_cost::HuberWeightFunction();
      break;
    default: std::runtime_error("Weight function unknown");
  }
  use_weights_ = true;
}

template <int D, typename T>
void vk::NLLSSolver<D, T>::
setInitialChi2(double c)
{
  if(!is_initial_chi2_provided_)
  {
    chi2_ = c;
    is_initial_chi2_provided_ = true;
  }
}

template <int D, typename T>
void vk::NLLSSolver<D, T>::
reset()
{
  mu_ = mu_init_;
  nu_ = nu_init_;
  n_meas_ = 0;
  n_iter_ = n_iter_init_;
  iter_ = 0;
  stop_ = false;
  is_initial_chi2_provided_ = false;
}

template <int D, typename T>
inline const double& vk::NLLSSolver<D, T>::
getChi2() const
{
  return chi2_;
}

#endif /* LM_SOLVER_IMPL_HPP_ */
