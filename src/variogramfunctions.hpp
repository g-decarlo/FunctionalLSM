// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#ifndef LOCALLY_STATIONARY_MODES_VARIOGRAMFUNCTIONS
#define LOCALLY_STATIONARY_MODES_VARIOGRAMFUNCTIONS

#include "traits.hpp"
#include <string>
#include <memory>

namespace LocallyStationaryModels {
class VariogramFunction {
protected:
  double compute_anisotropic_h(
      const double& lambda1, const double& lambda2, const double& phi, const double& x, const double& y);
  
public:
  VariogramFunction() = default;
  virtual ~VariogramFunction() = default; // Virtual destructor for base class
  
  // Base virtual functions that will be overridden by derived classes
  virtual double operator()(const cd::vector& params1, const cd::vector& params2, const double& x, const double& y);
  virtual double operator()(const cd::vector& params, const double& x, const double& y);
  virtual double correlation(const cd::vector& params1, const cd::vector& params2, const double& x, const double& y);
  virtual double correlation(const cd::vector& params, const double& x, const double& y) = 0; // Pure virtual, must be implemented
};

class Exponential : public VariogramFunction {
public:
  Exponential() = default;
  double correlation(const cd::vector& params, const double& x, const double& y) override;
};

class ExponentialNugget : public VariogramFunction {
public:
  ExponentialNugget() = default;
  double correlation(const cd::vector& params, const double& x, const double& y) override;
  double operator()(const cd::vector& params, const double& x, const double& y) override;
};

class Matern : public VariogramFunction {
public:
  Matern() = default;
  double correlation(const cd::vector& params, const double& x, const double& y) override;
};

class MaternNuFixed : public VariogramFunction {
private:
  double m_nu = 0.5;
public:
  MaternNuFixed(const double& nu) : m_nu(nu) {};
  double correlation(const cd::vector& params, const double& x, const double& y) override;
};

class Gaussian : public VariogramFunction {
public:
  Gaussian() = default;
  double correlation(const cd::vector& params, const double& x, const double& y) override;
};

class Nugget : public VariogramFunction {
public:
  Nugget() = default;
  double correlation(const cd::vector& params, const double& x, const double& y) override;
};

class MaternNuNugget : public VariogramFunction {
private:
  double m_nu = 0.5;
public:
  MaternNuNugget(const double& nu) : m_nu(nu) {};
  double correlation(const cd::vector& params, const double& x, const double& y) override;
  double operator()(const cd::vector& params, const double& x, const double& y) override;
};

// Factory function to create the correct model object
std::shared_ptr<VariogramFunction> make_variogramiso(const std::string& id);

} // namespace LocallyStationaryModels

#endif // LOCALLY_STATIONARY_MODES_VARIOGRAMFUNCTIONS
