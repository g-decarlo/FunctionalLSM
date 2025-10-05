// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#ifndef LOCALLY_STATIONARY_MODES_VARIOGRAMFUNCTIONS
#define LOCALLY_STATIONARY_MODES_VARIOGRAMFUNCTIONS

#include "traits.hpp"
#include <Rcpp.h>
#include <string>
#include <memory>
#include <algorithm>

namespace LocallyStationaryModels {

// Forward declaration for cd::vector if traits.hpp is minimal
namespace cd {
using vector = Eigen::VectorXd;
using matrix = Eigen::MatrixXd;
}

class VariogramFunction {
protected:
  double compute_anisotropic_h(
      const double& lambda1, const double& lambda2, const double& phi, const double& x, const double& y);
  
public:
  VariogramFunction() = default;
  virtual ~VariogramFunction() = default;
  
  virtual double operator()(const cd::vector& params1, const cd::vector& params2, const double& x, const double& y);
  virtual double operator()(const cd::vector& params, const double& x, const double& y);
  virtual double correlation(const cd::vector& params1, const cd::vector& params2, const double& x, const double& y);
  virtual double correlation(const cd::vector& params, const double& x, const double& y) = 0;
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



inline std::shared_ptr<VariogramFunction> make_variogramiso(const std::string& id) {
  std::string lower_id = id;
  std::transform(lower_id.begin(), lower_id.end(), lower_id.begin(), ::tolower);
  
  if (lower_id == "exponential" || lower_id == "esponenziale" || lower_id == "exp") {
    return std::make_shared<Exponential>();
  } 
  else if (lower_id == "matern") {
    return std::make_shared<Matern>();
  } 
  else if (lower_id == "gaussian") {
    return std::make_shared<Gaussian>();
  }
  else if (lower_id == "nugget") {
    return std::make_shared<Nugget>();
  }
  else if (lower_id == "exponentialnugget") {
    return std::make_shared<ExponentialNugget>();
  }
  else {
    std::string base_nugget_long = "maternununugget";
    std::string base_nugget_short = "maternunugget";
    std::string base_fixed_long = "maternunufixed";
    std::string base_fixed_short = "maternufixed";
    
    if (lower_id.rfind(base_nugget_long, 0) == 0 || lower_id.rfind(base_nugget_short, 0) == 0) {
      try {
        std::string nu_str;
        if (lower_id.rfind(base_nugget_long, 0) == 0) {
          nu_str = lower_id.substr(base_nugget_long.length());
        } else {
          nu_str = lower_id.substr(base_nugget_short.length());
        }
        
        // More robustly find the start of the number
        size_t first_digit = nu_str.find_first_of("0123456789.");
        if (first_digit == std::string::npos) {
          Rcpp::stop("No numeric value for nu found in '" + id + "'.");
        }
        nu_str = nu_str.substr(first_digit);
        
        double NU = std::stod(nu_str);
        return std::make_shared<MaternNuNugget>(NU);
      } catch (const std::exception& e) {
        Rcpp::stop("Failed to parse nu value from '" + id + "'. Details: " + e.what());
      }
    }
    else if (lower_id.rfind(base_fixed_long, 0) == 0 || lower_id.rfind(base_fixed_short, 0) == 0) {
      try {
        std::string nu_str;
        if (lower_id.rfind(base_fixed_long, 0) == 0) {
          nu_str = lower_id.substr(base_fixed_long.length());
        } else {
          nu_str = lower_id.substr(base_fixed_short.length());
        }
        
        size_t first_digit = nu_str.find_first_of("0123456789.");
        if (first_digit == std::string::npos) {
          Rcpp::stop("No numeric value for nu found in '" + id + "'.");
        }
        nu_str = nu_str.substr(first_digit);
        
        double NU = std::stod(nu_str);
        return std::make_shared<MaternNuFixed>(NU);
      } catch (const std::exception& e) {
        Rcpp::stop("Failed to parse nu value from '" + id + "'. Details: " + e.what());
      }
    }
    else {
      Rcpp::stop("Invalid variogram model id: '" + id + "'");
    }
  }
  // This part of the code is now unreachable, ensuring no undefined behavior.
  return nullptr;
}

} // namespace LocallyStationaryModels

#endif // LOCALLY_STATIONARY_MODES_VARIOGRAMFUNCTIONS
