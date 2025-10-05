// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
<<<<<<< HEAD
#include <Rmath.h>
=======

>>>>>>> main
#include "variogramfunctions.hpp"

namespace LocallyStationaryModels {
using namespace cd;

double VariogramFunction::compute_anisotropic_h(
    const double& lambda1, const double& lambda2, const double& phi, const double& x, const double& y)
{
<<<<<<< HEAD
  double xx = x * x;
  double yy = y * y;
  double xy = x * y;
  
  return sqrt((lambda2 * lambda2 * xx * cos(phi) * cos(phi) + lambda1 * lambda1 * yy * cos(phi) * cos(phi)
                 + lambda1 * lambda1 * xx * sin(phi) * sin(phi) + lambda2 * lambda2 * yy * sin(phi) * sin(phi)
                 - lambda1 * lambda1 * xy * sin(2 * phi) + lambda2 * lambda2 * xy * sin(2 * phi))
                 / (lambda1 * lambda1 * lambda2 * lambda2));
=======
    double xx = x * x;
    double yy = y * y;
    double xy = x * y;

    return sqrt((lambda2 * lambda2 * xx * cos(phi) * cos(phi) + lambda1 * lambda1 * yy * cos(phi) * cos(phi)
                    + lambda1 * lambda1 * xx * sin(phi) * sin(phi) + lambda2 * lambda2 * yy * sin(phi) * sin(phi)
                    - lambda1 * lambda1 * xy * sin(2 * phi) + lambda2 * lambda2 * xy * sin(2 * phi))
        / (lambda1 * lambda1 * lambda2 * lambda2));
>>>>>>> main
}

double VariogramFunction::correlation(
    const cd::vector& params1, const cd::vector& params2, const double& x, const double& y)
<<<<<<< HEAD
{
  double lambda1_1 = params1[0];
  double lambda2_1 = params1[1];
  double phi_1 = params1[2];
  double lambda1_2 = params2[0];
  double lambda2_2 = params2[1];
  double phi_2 = params2[2];
  double sigma_1 = params1[3];
  double sigma_2 = params2[3];
  
  cd::matrix rot1(2,2), rot2(2,2), eig1(2,2), eig2(2,2);
  
  rot1 << cos(phi_1), -sin(phi_1)  , +sin(phi_1), cos(phi_1);
  rot2 << cos(phi_2), -sin(phi_2)  , +sin(phi_2), cos(phi_2);
  eig1 << 1.0/lambda1_1, 0, 0, 1.0/lambda2_1;
  eig2 << 1.0/lambda1_2, 0, 0, 1.0/lambda2_2;
  
  cd::matrix anis1(rot1*eig1*eig1*rot1.transpose());
  cd::matrix anis2(rot2*eig2*eig2*rot2.transpose());
  Eigen::Matrix<double, 2, 2> anistot((anis1+anis2)/2);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2,2> > eigen(anistot);
  
  cd::vector s(2);
  s << x, y;
  
  return 2 * sqrt( (lambda1_1 * lambda1_2 * lambda2_1 * lambda2_2) / ( (2*anistot).determinant() ) ) * exp( -sqrt( s.dot(anistot.inverse() * s) ) );
=======
{   
    double lambda1_1 = params1[0];
    double lambda2_1 = params1[1];
    double phi_1 = params1[2];
    double lambda1_2 = params2[0];
    double lambda2_2 = params2[1];
    double phi_2 = params2[2];
    double sigma_1 = params1[3];
    double sigma_2 = params2[3];

    cd::matrix rot1(2,2), rot2(2,2), eig1(2,2), eig2(2,2);

    rot1 << cos(phi_1), -sin(phi_1)  , +sin(phi_1), cos(phi_1);
    rot2 << cos(phi_2), -sin(phi_2)  , +sin(phi_2), cos(phi_2);
    eig1 << lambda1_1, 0, 0, lambda2_1;
    eig2 << lambda1_2, 0, 0, lambda2_2;

    cd::matrix anis1(rot1*eig1*eig1*rot1.transpose());
    cd::matrix anis2(rot2*eig2*eig2*rot2.transpose());
    Eigen::Matrix<double, 2, 2> anistot(anis1/2+anis2/2);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2,2> > eigen(anistot);
    Eigen::RowVectorXd paramaniso(4);
    double phi = acos(abs(eigen.eigenvectors()(1,0)));
    double lambda1 = sqrt(eigen.eigenvalues()(0));
    double lambda2 = sqrt(eigen.eigenvalues()(1));

    paramaniso << lambda1, lambda2, phi , sqrt(sigma_1*sigma_2);
    cd::vector s(2);
    s << x, y;
    return 2*(sqrt(lambda1_1*lambda1_2*lambda2_1*lambda2_2/((2*anistot).determinant())))*(exp(-sqrt(s.dot(anistot.inverse()*s))));
    //return 2*(sqrt(lambda1_1*lambda1_2*lambda2_1*lambda2_2/((2*anistot).determinant())))*(this->correlation(paramaniso,x,y)); TODO Temporary fix
>>>>>>> main
}

double VariogramFunction::operator()(
    const cd::vector& params1, const cd::vector& params2, const double& x, const double& y)
<<<<<<< HEAD
{
  double sigma_1 = params1[3];
  double sigma_2 = params2[3];
  return sigma_1*sigma_2*(1.0 - this->correlation(params1, params2, x, y));
}

double VariogramFunction::operator()(const cd::vector& params, const double& x, const double& y)
{
  return params[3]*params[3]*(1-this->correlation(params, x, y));
=======
{   
    double sigma_1 = params1[3];
    double sigma_2 = params2[3];
    return sigma_1*sigma_2*(this->correlation(params1, params2, x, y));
}

double VariogramFunction::operator()(const cd::vector& params, const double& x, const double& y)
{   
    return params[3]*params[3]*(1-this->correlation(params, x, y));
>>>>>>> main
}


double Exponential::correlation(const cd::vector& params, const double& x, const double& y)
<<<<<<< HEAD
{
  double lambda1 = params[0];
  double lambda2 = params[1];
  double phi = params[2];
  double h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
  return exp(-h);
}

double ExponentialNugget::correlation(const cd::vector& params, const double& x, const double& y)
{
  double lambda1 = params[0];
  double lambda2 = params[1];
  double phi = params[2];
  double sigma = params[3];
  double tau2 = params[4];
  double h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
  if (std::abs(x) < 1e-9 && std::abs(y) < 1e-9) {
    return 1.0;
  }
  
  return (1.0 - tau2/(sigma*sigma+tau2))*exp(-h);
  
=======
{   
    double lambda1 = params[0];
    double lambda2 = params[1];
    double phi = params[2];
    double sigma = params[3];
    double h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
    return exp(-h);
}

double ExponentialNugget::correlation(const cd::vector& params, const double& x, const double& y)
{   
    double lambda1 = params[0];
    double lambda2 = params[1];
    double phi = params[2];
    double sigma = params[3];
    double tau2 = params[4];
    double h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
    if (std::abs(x) < Tolerances::min_norm && std::abs(y) < Tolerances::min_norm) {
        return 1;
    }

    return (1 - tau2/(sigma*sigma+tau2))*exp(-h);
        
>>>>>>> main
}

double ExponentialNugget::operator()(const cd::vector& params, const double& x, const double& y)
{
<<<<<<< HEAD
  return (params[3]*params[3]+params[4])*(1-this->correlation(params, x, y));
=======
    return (params[3]*params[3]+params[4])*(1-this->correlation(params, x, y));
>>>>>>> main
}

double Matern::correlation(const cd::vector& params, const double& x, const double& y)
{
<<<<<<< HEAD
  double lambda1 = params[0];
  double lambda2 = params[1];
  double phi = params[2];
  double nu = params[4];
  
  if (std::abs(x) < 1e-9 && std::abs(y) < 1e-9) {
    return 1.0;
  }
  
  double h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
  double h_nu = std::sqrt(2 * nu) * h;
  return (std::pow(h_nu, nu) * bessel_k(h_nu, nu, 1.0)) / (std::tgamma(nu) * std::pow(2.0, nu - 1.0));
=======
    double lambda1 = params[0];
    double lambda2 = params[1];
    double phi = params[2];
    double sigma = params[3];
    double nu = params[4];

    if (std::abs(x) < Tolerances::min_norm && std::abs(y) < Tolerances::min_norm) {
        return Tolerances::infinity;
    }

    double h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
    return 
            - std::pow(std::sqrt(2 * nu) * h, nu) * std::cyl_bessel_k(nu, std::sqrt(2 * nu) * h)
                / (std::tgamma(nu) * std::pow(2, nu - 1));
>>>>>>> main
}

double Nugget::correlation(const cd::vector& params, const double& x, const double& y)
{
<<<<<<< HEAD
  if (std::abs(x) < 1e-9 && std::abs(y) < 1e-9) {
    return 1;
  }
  return 0;
=======
    if (std::abs(x) < Tolerances::min_norm && std::abs(y) < Tolerances::min_norm) {
        return 1;
    }
    return 0;
>>>>>>> main
}

double MaternNuFixed::correlation(const cd::vector& params, const double& x, const double& y)
{
<<<<<<< HEAD
  double lambda1 = params[0];
  double lambda2 = params[1];
  double phi = params[2];
  double nu = m_nu;
  
  if (std::abs(x) < 1e-9 && std::abs(y) < 1e-9) {
    return 1.0;
  }
  
  double h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
  double h_nu = std::sqrt(2 * nu) * h;
  return (std::pow(h_nu, nu) * bessel_k(h_nu, nu, 1.0)) / (std::tgamma(nu) * std::pow(2.0, nu - 1.0));
=======
    double lambda1 = params[0];
    double lambda2 = params[1];
    double phi = params[2];
    double sigma = params[3];
    double nu = m_nu;

    if (std::abs(x) < Tolerances::min_norm && std::abs(y) < Tolerances::min_norm) {
        return Tolerances::infinity;
    }

    double h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
    return
             std::pow(std::sqrt(2 * nu) * h, nu) * std::cyl_bessel_k(nu, std::sqrt(2 * nu) * h)
                / (std::tgamma(nu) * std::pow(2, nu - 1));
>>>>>>> main
}

double MaternNuNugget::correlation(const cd::vector& params, const double& x, const double& y)
{
<<<<<<< HEAD
  double lambda1 = params[0];
  double lambda2 = params[1];
  double phi = params[2];
  double sigma = params[3];
  double tau2 = params[4];
  double nu = m_nu;
  
  if (std::abs(x) < 1e-9 && std::abs(y) < 1e-9) {
    return 1.0;
  }
  
  double h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
  double h_nu = std::sqrt(2 * nu) * h;
  
  return (1.0 - tau2/(sigma*sigma+tau2)) * (std::pow(h_nu, nu) * bessel_k(h_nu, nu, 1.0)) / (std::tgamma(nu) * std::pow(2.0, nu - 1.0));
  
  
=======
    double lambda1 = params[0];
    double lambda2 = params[1];
    double phi = params[2];
    double sigma = params[3];
    double tau2 = params[4];
    double nu = m_nu;

    if (std::abs(x) < Tolerances::min_norm && std::abs(y) < Tolerances::min_norm) {
        return 1;
    }

    double h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
    
    return
             (1-tau2/(sigma*sigma+tau2))*std::pow(std::sqrt(2 * nu) * h, nu) * std::cyl_bessel_k(nu, std::sqrt(2 * nu) * h)
                / (std::tgamma(nu) * std::pow(2, nu - 1));
    

>>>>>>> main
}

double MaternNuNugget::operator()(const cd::vector& params, const double& x, const double& y)
{
<<<<<<< HEAD
  return (params[3]*params[3]+params[4])*(1-this->correlation(params, x, y));
=======
    return (params[3]*params[3]+params[4])*(1-this->correlation(params, x, y));
>>>>>>> main
}

double Gaussian::correlation(const cd::vector& params, const double& x, const double& y)
{
<<<<<<< HEAD
  double lambda1 = params[0];
  double lambda2 = params[1];
  double phi = params[2];
  
  double h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
  return exp(-h * h);
=======
    double lambda1 = params[0];
    double lambda2 = params[1];
    double phi = params[2];
    double sigma = params[3];

    double h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
    return exp(-h * h);
>>>>>>> main
}

std::shared_ptr<VariogramFunction> make_variogramiso(const std::string& id)
{
<<<<<<< HEAD
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
=======
    if (id == "exponential" || id == "esponenziale") {
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
>>>>>>> main
