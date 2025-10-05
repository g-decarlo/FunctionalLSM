// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#include "coregionalization.hpp" 
#include <iostream>   

namespace LocallyStationaryModels {
using namespace cd;

cd::matrix CrossCovariance::operator()(const std::vector<cd::vector>& params, const double& x, const double& y) const{
    cd::matrix A_mat(m_dim,m_dim);
    double Rns;
    cd::matrix CrossVario = Eigen::MatrixXd::Zero(m_dim,m_dim);
    for(size_t i = 0; i < m_r; i++){
        A_mat = matricize(params[i].segment(3,m_dim*(m_dim+1)/2));
        Rns = vario_functions[i]->correlation(params[i].segment(0,4), x, y);//corr including 3, to be fixed TODO
        CrossVario += A_mat*A_mat.transpose()*(1-Rns);
    }
    return CrossVario;
    
}

cd::matrix CrossCovariance::operator()(const std::vector<cd::vector>& params1, const std::vector<cd::vector>& params2, const double& x, const double& y) const{
    cd::matrix A_mat1(m_dim,m_dim), A_mat2(m_dim,m_dim);
    double Rns;
    cd::matrix CrossVario = Eigen::MatrixXd::Zero(m_dim,m_dim);
    for(size_t i = 0; i < m_r; i++){
        A_mat1 = matricize(params1[i].segment(3,m_dim*(m_dim+1)/2));
        A_mat2 = matricize(params2[i].segment(3,m_dim*(m_dim+1)/2));
        Rns = vario_functions[i]->correlation(params1[i].segment(0,4), params2[i].segment(0,4), x, y);//corr including 3, to be fixed TODO
        CrossVario += A_mat2*A_mat1.transpose()*Rns;
    }
    return CrossVario;
}

CrossCovariance::CrossCovariance(const std::vector<std::string>& vario_ids, const size_t& dim, const size_t& r): m_r(r), m_dim(dim){
    for(size_t i = 0; i < m_r; i++){
        this->vario_functions.push_back( make_variogramiso(vario_ids[i]) );
    }
}

cd::matrix CrossCovariance::matricize(const cd::vector& vec) const{
    cd::matrix mat = Eigen::MatrixXd::Zero(m_dim,m_dim);
    for(size_t i = 0, k = 0; i < m_dim; i++ ){
        for(size_t j = 0; j <= i; j++ ){
           mat(i,j) = vec(k);
           k++; 
        }
    }
    return mat;
}



} // namespace LocallyStationaryModels
