// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#include <Rmath.h> 
#include "variogramfunctions.hpp"

namespace LocallyStationaryModels {
using namespace cd;

double VariogramFunction::compute_anisotropic_h(
    const double& lambda1, const double& lambda2, const double& phi, const double& x, const double& y)
{
  double xx = x * x;
  double yy = y * y;
  double xy = x * y;
  double denominator = lambda1 * lambda1 * lambda2 * lambda2;
  if (denominator < 1e-12) { return 1e12; }
  double argument = (lambda2 * lambda2 * xx * cos(phi) * cos(phi) + lambda1 * lambda1 * yy * cos(phi) * cos(phi)
                       + lambda1 * lambda1 * xx * sin(phi) * sin(phi) + lambda2 * lambda2 * yy * sin(phi) * sin(phi)
                       - lambda1 * lambda1 * xy * sin(2 * phi) + lambda2 * lambda2 * xy * sin(2 * phi))
                       / denominator;
                       return sqrt(std::max(0.0, argument));
}

double VariogramFunction::correlation(
    const cd::vector& params1, const cd::vector& params2, const double& x, const double& y)
{
  if (params1.size() < 3 || params2.size() < 3) {
    Rcpp::stop("Parameter vectors for non-stationary correlation must have at least 3 elements (lambda1, lambda2, phi).");
  }
  double lambda1_1 = params1[0]; double lambda2_1 = params1[1]; double phi_1 = params1[2];
  double lambda1_2 = params2[0]; double lambda2_2 = params2[1]; double phi_2 = params2[2];
  cd::matrix rot1(2,2), rot2(2,2), eig1(2,2), eig2(2,2);
  rot1 << cos(phi_1), -sin(phi_1)  , +sin(phi_1), cos(phi_1);
  rot2 << cos(phi_2), -sin(phi_2)  , +sin(phi_2), cos(phi_2);
  eig1 << 1.0/lambda1_1, 0, 0, 1.0/lambda2_1;
  eig2 << 1.0/lambda1_2, 0, 0, 1.0/lambda2_2;
  cd::matrix anis1(rot1*eig1*eig1*rot1.transpose());
  cd::matrix anis2(rot2*eig2*eig2*rot2.transpose());
  Eigen::Matrix<double, 2, 2> anistot((anis1+anis2)/2);
  
  double det = (2*anistot).determinant();
  if (det <= 1e-9) {

    return 0.0;
  }
  
  cd::vector s(2); s << x, y;
  double h_squared = s.dot(anistot.inverse() * s);
  double h = sqrt(std::max(0.0, h_squared));
  return 2.0 * sqrt( (lambda1_1 * lambda1_2 * lambda2_1 * lambda2_2) / det ) * exp(-h);
}

double VariogramFunction::operator()(
    const cd::vector& params1, const cd::vector& params2, const double& x, const double& y)
{
  if (params1.size() < 4 || params2.size() < 4) {
    Rcpp::stop("Parameter vectors for non-stationary variogram must have at least 4 elements (lambda1, lambda2, phi, sigma).");
  }
  double sigma_1 = params1[3]; double sigma_2 = params2[3];
  return sigma_1*sigma_2*(1.0 - this->correlation(params1, params2, x, y));
}

double VariogramFunction::operator()(const cd::vector& params, const double& x, const double& y)
{
  if (params.size() < 4) {
    Rcpp::stop("Parameter vector must have at least 4 elements (lambda1, lambda2, phi, sigma).");
  }
  return params[3]*params[3]*(1-this->correlation(params, x, y));
}

double Exponential::correlation(const cd::vector& params, const double& x, const double& y)
{
  if (params.size() < 4) {
    Rcpp::stop("Exponential model requires 4 parameters (lambda1, lambda2, phi, sigma).");
  }
  return exp(-compute_anisotropic_h(params[0], params[1], params[2], x, y));
}

double ExponentialNugget::correlation(const cd::vector& params, const double& x, const double& y)
{
  if (params.size() < 5) {
    Rcpp::stop("ExponentialNugget model requires 5 parameters (lambda1, lambda2, phi, sigma, tau).");
  }
  double h = compute_anisotropic_h(params[0], params[1], params[2], x, y);
  if (h < 1e-9) { return 1.0; }
  
  double sigma2 = params[3] * params[3];
  double tau2 = params[4];
  double total_sill = sigma2 + tau2;
  if (total_sill < 1e-9) { return 0.0; }
  
  return (1.0 - tau2 / total_sill) * exp(-h);
}

double ExponentialNugget::operator()(const cd::vector& params, const double& x, const double& y)
{
  if (params.size() < 5) {
    Rcpp::stop("ExponentialNugget model requires 5 parameters (lambda1, lambda2, phi, sigma, tau).");
  }
  return (params[3]*params[3]+params[4])*(1-this->correlation(params, x, y));
}

double Matern::correlation(const cd::vector& params, const double& x, const double& y)
{
  if (params.size() < 5) {
    Rcpp::stop("Matern model requires 5 parameters (lambda1, lambda2, phi, sigma, nu).");
  }
  double h = compute_anisotropic_h(params[0], params[1], params[2], x, y);
  if (h < 1e-9) { return 1.0; }
  double nu = params[4];
  if (nu <= 0) {
    Rcpp::stop("Matern smoothness parameter 'nu' (params[4]) must be positive.");
  }
  double h_nu = std::sqrt(2 * nu) * h;
  double result = (std::pow(h_nu, nu) * R::bessel_k(h_nu, nu, 1.0) * exp(-h_nu)) / (std::tgamma(nu) * std::pow(2.0, nu - 1.0));
  if (!std::isfinite(result)) { return 0.0; }
  return result;
}

double Nugget::correlation(const cd::vector& params, const double& x, const double& y)
{
  if (params.size() < 4) {
    Rcpp::stop("Nugget model requires 4 parameters (lambda1, lambda2, phi, sigma).");
  }
  double h = compute_anisotropic_h(params[0], params[1], params[2], x, y);
  if (h < 1e-9) { return 1.0; }
  return 0.0;
}

double MaternNuFixed::correlation(const cd::vector& params, const double& x, const double& y)
{
  if (params.size() < 4) {
    Rcpp::stop("MaternNuFixed model requires 4 parameters (lambda1, lambda2, phi, sigma).");
  }
  double h = compute_anisotropic_h(params[0], params[1], params[2], x, y);
  if (h < 1e-9) { return 1.0; }
  double nu = m_nu;
  double h_nu = std::sqrt(2 * nu) * h;
  double result = (std::pow(h_nu, nu) * R::bessel_k(h_nu, nu, 1.0) * exp(-h_nu)) / (std::tgamma(nu) * std::pow(2.0, nu - 1.0));
  if (!std::isfinite(result)) { return 0.0; }
  return result;
}

double MaternNuNugget::correlation(const cd::vector& params, const double& x, const double& y)
{
  if (params.size() < 5) {
    Rcpp::stop("MaternNuNugget model requires 5 parameters (lambda1, lambda2, phi, sigma, tau).");
  }
  double h = compute_anisotropic_h(params[0], params[1], params[2], x, y);
  if (h < 1e-9) { return 1.0; }
  
  double sigma2 = params[3] * params[3];
  double tau2 = params[4];
  double total_sill = sigma2 + tau2;
  if (total_sill < 1e-9) { return 0.0; }
  
  double nu = m_nu;
  double h_nu = std::sqrt(2 * nu) * h;
  double matern_corr = (std::pow(h_nu, nu) * R::bessel_k(h_nu, nu, 1.0) * exp(-h_nu)) / (std::tgamma(nu) * std::pow(2.0, nu - 1.0));
  if (!std::isfinite(matern_corr)) { matern_corr = 0.0; }
  
  return (1.0 - tau2 / total_sill) * matern_corr;
}

double MaternNuNugget::operator()(const cd::vector& params, const double& x, const double& y)
{
  if (params.size() < 5) {
    Rcpp::stop("MaternNuNugget model requires 5 parameters (lambda1, lambda2, phi, sigma, tau).");
  }
  return (params[3]*params[3]+params[4])*(1-this->correlation(params, x, y));
}

double Gaussian::correlation(const cd::vector& params, const double& x, const double& y)
{
  if (params.size() < 4) {
    Rcpp::stop("Gaussian model requires 4 parameters (lambda1, lambda2, phi, sigma).");
  }
  double h = compute_anisotropic_h(params[0], params[1], params[2], x, y);
  return exp(-h * h);
}


} // namespace LocallyStationaryModels

