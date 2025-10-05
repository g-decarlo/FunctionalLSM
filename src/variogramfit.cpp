// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#include "variogramfit.hpp"
#include "variogramfunctions.hpp" // Ensure header is included for declarations
#include "LBFGS/LBFGSB.h"
#include <iostream>

namespace LocallyStationaryModels {
using namespace cd;
using namespace LBFGSpp;

TargetFunction::TargetFunction(const cd::matrixptr& empiricvariogram, const cd::matrixptr& squaredweights,
                               const cd::vectorptr& mean_x, const cd::vectorptr& mean_y, const size_t& x0, const std::vector<std::string>& id, const size_t& dim)
  : m_empiricvariogram(empiricvariogram)
  , m_squaredweights(squaredweights)
  , m_mean_x(mean_x)
  , m_mean_y(mean_y)
  , m_x0(x0)
  , m_gammaisoptr(make_variogramiso(id[0]))
  , m_crosscov(id, dim, id.size())
  , m_r(id.size())
  , m_dim(dim) {};

double TargetFunction::operator()(const cd::vector& params)
{
  if (m_dim != 1){
    return this->operator()(params, m_r);
  }
  
  VariogramFunction& gammaiso = *(m_gammaisoptr);
  vector w = m_squaredweights->row(m_x0);
  vector truegamma(m_empiricvariogram->rows());
  
  double sill;
  bool has_nugget = (params.size() == 5); // Models with nugget have 5 params (lambda1, lambda2, phi, sigma, tau2)
  
  if (has_nugget) {
    double sigma2 = params[3] * params[3];
    double tau2 = params[4];
    sill = sigma2 + tau2;
  } else {
    double sigma2 = params[3] * params[3];
    sill = sigma2;
  }
  
  for (size_t h = 0; h < truegamma.size(); ++h) {
    double correlation = gammaiso.correlation(params, m_mean_x->operator[](h), m_mean_y->operator[](h));
    truegamma[h] = sill * (1.0 - correlation);
  }
  
  vector empiricgamma = m_empiricvariogram->col(m_x0);
  return w.dot((truegamma - empiricgamma).cwiseProduct(truegamma - empiricgamma));
}

double TargetFunction::operator()(const cd::vector& params, const size_t& r)
{
  vector w = m_squaredweights->row(m_x0);
  size_t rows = m_empiricvariogram->rows();
  matrix truegamma(rows, m_dim);
  matrix empiricgamma = m_empiricvariogram->block(0,m_x0*m_dim,rows,m_dim);
  std::vector<vector> paramvec;
  unsigned n_params = 3; //TODO Fix variable number of parameters as a member of variogram function
  for (size_t i = 0, totparams = 0; i < r; i++ ){
    paramvec.push_back(params.segment(totparams, n_params + m_dim*(m_dim+1)/2));
    totparams += n_params + m_dim*(m_dim+1)/2;
  }
  for (size_t h = 0; h < rows/m_dim ; ++h) {
    truegamma.block(h*m_dim,0,m_dim,m_dim) = sqrt(w[h])*m_crosscov(paramvec, m_mean_x->operator[](h), m_mean_y->operator[](h));
    empiricgamma.block(h*m_dim,0,m_dim,m_dim) *= sqrt(w[h]);
  }
  return (truegamma - empiricgamma).squaredNorm();
}

double TargetFunction::operator()(const cd::vector& params, vector& grad)
{
  if (m_dim != 1){
    return this->operator()(params, grad, m_r);
  }
  
  for (size_t i = 0; i < params.size(); ++i) {
    vector paramsdeltaplus(params);
    vector paramsdeltaminus(params);
    
    double increment = 1e-6 * std::max(std::abs(params[i]), 1.0);
    if (increment < 1e-9) increment = 1e-9; 
    
    paramsdeltaplus[i] += increment;
    paramsdeltaminus[i] -= increment;
    
    double f_plus = TargetFunction::operator()(paramsdeltaplus);
    double f_minus = TargetFunction::operator()(paramsdeltaminus);
    
    if (!std::isfinite(f_plus) || !std::isfinite(f_minus)) {
      grad[i] = 0.0;
    } else {
      grad[i] = (f_plus - f_minus) / (2.0 * increment);
    }
  }
  return TargetFunction::operator()(params);
}

double TargetFunction::operator()(const cd::vector& params, vector& grad, const size_t& r)
{
  for (size_t i = 0; i < params.size(); ++i) {
    vector paramsdeltaplus(params);
    vector paramsdeltaminus(params);
    
    double increment = 1e-6 * std::max(std::abs(params[i]), 1.0);
    if (increment < 1e-9) increment = 1e-9;
    
    paramsdeltaplus[i] += increment;
    paramsdeltaminus[i] -= increment;
    
    double f_plus = TargetFunction::operator()(paramsdeltaplus, r);
    double f_minus = TargetFunction::operator()(paramsdeltaminus, r);
    
    if (!std::isfinite(f_plus) || !std::isfinite(f_minus)) {
      grad[i] = 0.0;
    } else {
      grad[i] = (f_plus - f_minus) / (2.0 * increment);
    }
  }
  return TargetFunction::operator()(params, r);
}

Opt::Opt(const cd::matrixptr& empiricvariogram, const cd::matrixptr& squaredweights, const cd::vectorptr& mean_x,
         const cd::vectorptr& mean_y, const size_t& dim, const std::vector<std::string>& id, const cd::vector& initialparameters,
         const cd::vector& lowerbound, const cd::vector& upperbound)
  : m_empiricvariogram(empiricvariogram)
  , m_squaredweights(squaredweights)
  , m_mean_x(mean_x)
  , m_mean_y(mean_y)
  , m_id(id)
  , m_initialparameters(initialparameters)
  , m_lowerbound(lowerbound)
  , m_upperbound(upperbound)
  , m_dim(dim)
{
  m_solutions = std::make_shared<matrix>(matrix::Zero(m_empiricvariogram->cols()/m_dim, m_initialparameters.size()));
};

vector Opt::findonesolution(const size_t& pos) const
{
  TargetFunction fun(m_empiricvariogram, m_squaredweights, m_mean_x, m_mean_y, pos, m_id, m_dim);
  
  LBFGSBParam<double> param;
  param.epsilon = 1e-6;
  param.max_iterations = 100;
  LBFGSBSolver<double> solver(param);
  
  Eigen::VectorXd lb(m_lowerbound);
  Eigen::VectorXd ub(m_upperbound);
  
  cd::vector x(m_initialparameters);
  double fx;
  
  try {
    solver.minimize(fun, x, fx, lb, ub);
  } catch (const std::exception& e) {
    std::cerr << "L-BFGS-B optimization failed at anchor point " << pos << ": " << e.what() << std::endl;
    x = m_initialparameters; // Return initial parameters on failure
  }
  return x;
}

void Opt::findallsolutions()
{
#pragma omp parallel for
  for (size_t i = 0; i < m_empiricvariogram->cols()/m_dim; ++i) {
    vector sol = findonesolution(i);
    m_solutions->row(i) = sol;
  }
}

cd::matrixptr Opt::get_solutions() const { return m_solutions; }

} // namespace LocallyStationaryModels

