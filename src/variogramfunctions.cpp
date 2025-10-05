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
  cd::vector s(2); s << x, y;
  double h_squared = s.dot(anistot.inverse() * s);
  double h = sqrt(std::max(0.0, h_squared));
  return 2.0 * sqrt( (lambda1_1 * lambda1_2 * lambda2_1 * lambda2_2) / ( (2*anistot).determinant() ) ) * exp(-h);
}
double VariogramFunction::operator()(
    const cd::vector& params1, const cd::vector& params2, const double& x, const double& y)
{
  double sigma_1 = params1[3]; double sigma_2 = params2[3];
  return sigma_1*sigma_2*(1.0 - this->correlation(params1, params2, x, y));
}
double VariogramFunction::operator()(const cd::vector& params, const double& x, const double& y)
{
  return params[3]*params[3]*(1-this->correlation(params, x, y));
}
double Exponential::correlation(const cd::vector& params, const double& x, const double& y)
{
  return exp(-compute_anisotropic_h(params[0], params[1], params[2], x, y));
}
double ExponentialNugget::correlation(const cd::vector& params, const double& x, const double& y)
{
  double h = compute_anisotropic_h(params[0], params[1], params[2], x, y);
  if (h < 1e-9) { return 1.0; }
  return (1.0 - params[4]/(params[3]*params[3]+params[4]))*exp(-h);
}
double ExponentialNugget::operator()(const cd::vector& params, const double& x, const double& y)
{
  return (params[3]*params[3]+params[4])*(1-this->correlation(params, x, y));
}
double Matern::correlation(const cd::vector& params, const double& x, const double& y)
{
  double h = compute_anisotropic_h(params[0], params[1], params[2], x, y);
  if (h < 1e-9) { return 1.0; }
  double nu = params[4];
  double h_nu = std::sqrt(2 * nu) * h;
  double result = (std::pow(h_nu, nu) * bessel_k(h_nu, nu, 1.0) * exp(-h_nu)) / (std::tgamma(nu) * std::pow(2.0, nu - 1.0));
  if (!std::isfinite(result)) { return 0.0; }
  return result;
}
double Nugget::correlation(const cd::vector& params, const double& x, const double& y)
{
  double h = compute_anisotropic_h(params[0], params[1], params[2], x, y);
  if (h < 1e-9) { return 1.0; }
  return 0.0;
}
double MaternNuFixed::correlation(const cd::vector& params, const double& x, const double& y)
{
  double h = compute_anisotropic_h(params[0], params[1], params[2], x, y);
  if (h < 1e-9) { return 1.0; }
  double nu = m_nu;
  double h_nu = std::sqrt(2 * nu) * h;
  double result = (std::pow(h_nu, nu) * bessel_k(h_nu, nu, 1.0) * exp(-h_nu)) / (std::tgamma(nu) * std::pow(2.0, nu - 1.0));
  if (!std::isfinite(result)) { return 0.0; }
  return result;
}
double MaternNuNugget::correlation(const cd::vector& params, const double& x, const double& y)
{
  double h = compute_anisotropic_h(params[0], params[1], params[2], x, y);
  if (h < 1e-9) { return 1.0; }
  double sigma = params[3]; double tau2 = params[4]; double nu = m_nu;
  double h_nu = std::sqrt(2 * nu) * h;
  double matern_corr = (std::pow(h_nu, nu) * bessel_k(h_nu, nu, 1.0) * exp(-h_nu)) / (std::tgamma(nu) * std::pow(2.0, nu - 1.0));
  if (!std::isfinite(matern_corr)) { matern_corr = 0.0; }
  return (1.0 - tau2/(sigma*sigma+tau2)) * matern_corr;
}
double MaternNuNugget::operator()(const cd::vector& params, const double& x, const double& y)
{
  return (params[3]*params[3]+params[4])*(1-this->correlation(params, x, y));
}
double Gaussian::correlation(const cd::vector& params, const double& x, const double& y)
{
  double h = compute_anisotropic_h(params[0], params[1], params[2], x, y);
  return exp(-h * h);
}

// The implementation of make_variogramiso has been moved to the header file.

} // namespace LocallyStationaryModels

