// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#include <chrono>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "variogramfit.hpp"
#include "LBFGS/LBFGS.h"

namespace LocallyStationaryModels {

using namespace CppAD;
using namespace std::chrono;
using cd = CppAD::vector<double>;

// class VariogramFit
VariogramFit::VariogramFit(const Eigen::MatrixXd& empvar,
                           const Eigen::MatrixXd& sqweights,
                           const size_t& dim,
                           const std::vector<std::string>& id,
                           const Eigen::VectorXd& x,
                           const Eigen::VectorXd& y)
  : m_empvar(empvar)
  , m_sqweights(sqweights)
  , m_dim(dim)
  , m_crosscov(id, dim, id.size())
  , m_x(x)
  , m_y(y) {};

template <class Type>
Type VariogramFit::operator()(const CppAD::vector<Type>& parameters)
{
  Type cost = 0;
  matrix_t<Type> B(m_dim, m_dim);
  vector_t<Type> params(m_dim * (m_dim + 1) / 2);
  vector_t<Type> p(m_r * 4);
  matrix_t<Type> cross(m_dim, m_dim);
  
  for (size_t u = 0; u < m_r; u++) {
    p[4 * u] = parameters[m_r * m_dim * (m_dim + 1) / 2 + 4 * u];
    p[4 * u + 1] = parameters[m_r * m_dim * (m_dim + 1) / 2 + 4 * u + 1];
    p[4 * u + 2] = parameters[m_r * m_dim * (m_dim + 1) / 2 + 4 * u + 2];
    p[4 * u + 3] = parameters[m_r * m_dim * (m_dim + 1) / 2 + 4 * u + 3];
  }
  for (size_t i = 0; i < m_empvar.cols(); i++) {
    cross.setZero();
    for (size_t u = 0; u < m_r; u++) {
      size_t count = 0;
      for (size_t k = 0; k < m_dim; k++) {
        for (size_t l = k; l < m_dim; l++) {
          params[count] = parameters[u * m_dim * (m_dim + 1) / 2 + count];
          count++;
        }
      }
      B = m_crosscov.reconstruct_matrix(params, u);
      
      vector_t<Type> temp_p(4);
      temp_p[0] = p[4 * u];
      temp_p[1] = p[4 * u + 1];
      temp_p[2] = p[4 * u + 2];
      temp_p[3] = p[4 * u + 3];
      cross += B * (1 - m_crosscov.m_v[u]->correlation(temp_p, m_x[i], m_y[i]));
    }
    
    matrix_t<Type> diff(m_dim, m_dim);
    for (size_t k = 0; k < m_dim; k++) {
      for (size_t l = 0; l < m_dim; l++) {
        diff(k, l) = (cross(k, l) - m_empvar(i, k * m_dim + l));
      }
    }
    cost += (diff.array() * diff.array() * m_sqweights(i, 0)).sum();
  }
  return cost;
}

// class VariogramFitScalar
VariogramFitScalar::VariogramFitScalar(const Eigen::VectorXd& empvar,
                                       const Eigen::VectorXd& weights,
                                       const std::vector<std::string>& id,
                                       const Eigen::VectorXd& x,
                                       const Eigen::VectorXd& y)
  : m_empvar(empvar)
  , m_weights(weights)
  , m_x(x)
  , m_y(y)
{
  m_v = make_variogramiso(id[0]);
}

template <class Type>
Type VariogramFitScalar::operator()(const CppAD::vector<Type>& parameters)
{
  Type cost = 0;
  for (size_t i = 0; i < m_empvar.size(); i++) {
    cost += m_weights[i] * CppAD::pow(m_v->operator()(parameters, m_x[i], m_y[i]) - m_empvar[i], 2);
  }
  return cost;
}

// class VariogramFitLBFGS
VariogramFitLBFGS::VariogramFitLBFGS(const Eigen::VectorXd& empvar,
                                     const Eigen::VectorXd& weights,
                                     const std::vector<std::string>& id,
                                     const Eigen::VectorXd& x,
                                     const Eigen::VectorXd& y,
                                     const size_t& dim,
                                     const size_t& r)
  : m_empvar(empvar)
  , m_weights(weights)
  , m_x(x)
  , m_y(y)
  , m_crosscov(id, dim, r)
  , m_dim(dim)
  , m_r(r) {};

double VariogramFitLBFGS::operator()(const Eigen::VectorXd& parameters, Eigen::VectorXd& grad)
{
  double cost = 0;
  cd::vector ind(parameters.size());
  for (size_t i = 0; i < parameters.size(); i++) {
    ind[i] = parameters[i];
  }
  auto fun = [this](cd::vector p) {
    double cost = 0;
    Eigen::MatrixXd B(m_dim, m_dim);
    Eigen::VectorXd params(m_dim * (m_dim + 1) / 2);
    Eigen::VectorXd par(m_r * 4);
    Eigen::MatrixXd cross(m_dim, m_dim);
    
    for (size_t u = 0; u < m_r; u++) {
      par[4 * u] = p[m_r * m_dim * (m_dim + 1) / 2 + 4 * u];
      par[4 * u + 1] = p[m_r * m_dim * (m_dim + 1) / 2 + 4 * u + 1];
      par[4 * u + 2] = p[m_r * m_dim * (m_dim + 1) / 2 + 4 * u + 2];
      par[4 * u + 3] = p[m_r * m_dim * (m_dim + 1) / 2 + 4 * u + 3];
    }
    
    for (size_t i = 0; i < m_empvar.size(); i++) {
      cross.setZero();
      for (size_t u = 0; u < m_r; u++) {
        size_t count = 0;
        for (size_t k = 0; k < m_dim; k++) {
          for (size_t l = k; l < m_dim; l++) {
            params[count] = p[u * m_dim * (m_dim + 1) / 2 + count];
            count++;
          }
        }
        B = m_crosscov.reconstruct_matrix(params, u);
        Eigen::VectorXd temp_p(4);
        temp_p << par[4 * u], par[4 * u + 1], par[4 * u + 2], par[4 * u + 3];
        cross += B * (1 - m_crosscov.m_v[u]->correlation(temp_p, m_x[i], m_y[i]));
      }
      
      Eigen::MatrixXd diff(m_dim, m_dim);
      for (size_t k = 0; k < m_dim; k++) {
        for (size_t l = 0; l < m_dim; l++) {
          diff(k, l) = (cross(k, l) - m_empvar(i, k * m_dim + l));
        }
      }
      cost += (diff.array() * diff.array() * m_weights(i, 0)).sum();
    }
    return cost;
  };
  cd::vector g = CppAD::gradient(fun, ind);
  for (size_t i = 0; i < parameters.size(); i++) {
    grad[i] = g[i];
  }
  cost = fun(ind);
  return cost;
}

std::shared_ptr<VariogramFunction> make_variogramiso(const std::string& id)
{
  if (id == "exponential" || id == "esponenziale" || id == "exp") {
    return std::make_shared<Exponential>();
  }
  if (id == "matern" || id == "Matern") {
    return std::make_shared<Matern>();
  }
  if (id == "gaussian" || id == "Gaussian") {
    return std::make_shared<Gaussian>();
  }
  if (id == "nugget" || id == "Nugget") {
    return std::make_shared<Nugget>();
  }
  if (id == "exponentialnugget" || id == "ExponentialNugget") {
    return std::make_shared<ExponentialNugget>();
  }
  // using the following method we can set directly from R passing a string a constant value for nu
  if (id.substr(0, 13) == "maternNuFixed") {
    try {
      double NU = std::stod(id.substr(14));
      return std::make_shared<MaternNuFixed>(NU);
    } catch (std::exception& e) {
      return std::make_shared<Exponential>();
    }
  }
  
  if (id.substr(0, 14) == "maternNuNugget") {
    try {
      double NU = std::stod(id.substr(15));
      return std::make_shared<MaternNuNugget>(NU);
    } catch (std::exception& e) {
      return std::make_shared<Exponential>();
    }
  }
  return std::make_shared<Exponential>();
}
} // namespace LocallyStationaryModels
