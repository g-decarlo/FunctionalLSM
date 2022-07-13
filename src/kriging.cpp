// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#include "kriging.hpp"
#include <iostream>

namespace LocallyStationaryModels {
using namespace cd;

vectorind Predictor::build_neighbourhood(const cd::vector& pos) const
{
    vectorind n;
    for (size_t i = 0; i < m_data->rows(); ++i) {
        const vector& datapos = m_data->row(i);
        // if datapos is in a neighbourhood of radius m_b
        if ((pos - datapos).norm() < m_b) {
            n.push_back(i);
        }
    }
    return n;
}

vectorind Predictor::build_neighbourhood(const size_t& pos) const
{
    vectorind n;
    const vector& pospos = m_a->row(pos);
    for (size_t i = 0; i < m_data->rows(); ++i) {
        const vector& posi = m_data->row(i);
        // if pos is in a neighbourhood of radius m_b
        if ((pospos - posi).norm() < m_b) {
            n.push_back(i);
        }
    }
    return n;
}

cd::vector Predictor::build_eta(cd::vector& params, vectorind& neighbourhood) const
{
    
    size_t n = neighbourhood.size();
    matrix gamma(n, n);
    VariogramFunction& gammaiso = *(m_gammaisoptr);
    // compute gamma
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            const vector& posi = m_data->row(neighbourhood[i]);
            const vector& posj = m_data->row(neighbourhood[j]);
            cd::vector s = posi - posj;
            gamma(i, j) = gammaiso(params, s[0], s[1]);
        }
    }

    vector ones = vector::Ones(n);

    // compute eta
    vector gammaones(n);
    gammaones = gamma.fullPivHouseholderQr().solve(ones);
    double denominator = ones.dot(gammaones);
    vector eta = (gammaones) / denominator;
    return eta;
}

cd::matrix Predictor::build_etavec( cd::vector& params, vectorind& neighbourhood) const
{
    size_t n = neighbourhood.size();
    matrix gamma(Eigen::MatrixXd::Zero((n+1)*m_dim, (n+1)*m_dim));
    matrix eta((n+1)*m_dim, m_dim);
    std::vector<vector> paramvec;
    for (size_t i = 0, totparams = 0; i < m_id.size(); i++ ){
            unsigned n_params = 3;
                paramvec.push_back(params.segment(totparams, n_params + m_dim*(m_dim+1)/2));
            totparams += n_params + m_dim*(m_dim+1)/2;
            }
    // compute matrix
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            const vector& posi = m_data->row(neighbourhood[i]);
            const vector& posj = m_data->row(neighbourhood[j]);
            cd::vector s = posi - posj;
            gamma.block(i*m_dim, j*m_dim, m_dim, m_dim) = m_crosscov(paramvec, 0., 0.) - m_crosscov(paramvec, s[0], s[1]);
            gamma.block(n*m_dim, j*m_dim, m_dim, m_dim) = Eigen::MatrixXd::Identity(m_dim,m_dim);
            gamma.block(i*m_dim, n*m_dim, m_dim, m_dim) = Eigen::MatrixXd::Identity(m_dim,m_dim);

        }
    }
    // compute eta
    //matrix eta((n+1)*m_dim, m_dim);
    matrix zeroone(Eigen::MatrixXd::Zero((n+1)*m_dim,m_dim));
    zeroone.block(n*m_dim,0,m_dim,m_dim) = Eigen::MatrixXd::Identity(m_dim,m_dim);
    eta = gamma.fullPivHouseholderQr().solve(zeroone);
    return eta;
}

std::pair<cd::vector, matrix> Predictor::build_etakriging(const cd::vector& params, const cd::vector& pos) const
{
    size_t n = m_data->rows();

    vector etakriging(n);
    vector C0(n);
    matrix correlationmatrix(n, n);
    double sigma2 = params[3] * params[3];
    VariogramFunction& gammaiso = *(m_gammaisoptr);
    // compute the corralation matrix and C0
    for (size_t i = 0; i < n; ++i) {
        const vector& posi = m_data->row(i);
        cd::vector s0 = posi - pos;
        cd::vector paramsi = m_smt.smooth_vector(posi);
        C0(i) = gammaiso(params, paramsi, s0[0], s0[1]);
    }
    // compute etakriging
    etakriging = *m_kriging_matrix_inverse*C0;
    // compute the variance
    matrix krigingvariance(1,1);
    krigingvariance(0,0) = params(3) * params(3) - C0.transpose() * etakriging;
    return std::make_pair(etakriging, krigingvariance);
}

std::pair<cd::matrix, cd::matrix> Predictor::build_etakrigingvec(const cd::vector& params, const cd::vector& pos) const
{
    size_t n = m_data->rows();

    matrix etakriging(n*m_dim, m_dim);
    matrix C0(n*m_dim, m_dim);
    matrix correlationmatrix(n*m_dim, n*m_dim);
    std::vector<vector> paramvec;
    for (size_t i = 0, totparams = 0; i < m_id.size(); i++ ){
            unsigned n_params = 3;
                paramvec.push_back(params.segment(totparams, n_params + m_dim*(m_dim+1)/2));
            totparams += n_params + m_dim*(m_dim+1)/2;
            }
    // compute the corralation matrix and C0
    for (size_t i = 0; i < n; ++i) {
        const vector& posi = m_data->row(i);
        cd::vector s0 = posi - pos;
        cd::vector paramsi = m_smt.smooth_vector(posi);
        std::vector<vector> paramveci;
        for (size_t i = 0, totparams = 0; i < m_id.size(); i++ ){
            unsigned n_params = 3;
                paramveci.push_back(paramsi.segment(totparams, n_params + m_dim*(m_dim+1)/2));
            totparams += n_params + m_dim*(m_dim+1)/2;
            }
        C0.block(i*m_dim,0,m_dim,m_dim) = m_crosscov(paramveci, paramvec, s0[0], s0[1]);
    }
    // compute etakriging
    etakriging = *m_kriging_matrix_inverse*C0;
    // compute the variance
    matrix krigingvariance(m_crosscov(paramvec, paramvec, 0., 0.));
    for (size_t i = 0; i < n; ++i) {
        krigingvariance -= etakriging.block(i*m_dim,0,m_dim,m_dim).transpose()*C0.block(i*m_dim,0,m_dim,m_dim);
    }

    return std::make_pair(etakriging, krigingvariance);
}


template <> vector Predictor::predict_mean<cd::vector, cd::vector>(const cd::vector& pos) const
{
    return m_smt.smooth_vector(pos).tail(m_z->cols());//aggiustare per matern variabile
}

template <> vector Predictor::predict_mean<size_t, vector>(const size_t& pos) const
{
    // find the value of the parameters relative to the anchorpoint in row pos
    cd::vector params = m_smt.smooth_vector(m_a->row(pos));
    // find the anchropoints in its neighbourhood
    vectorind neighbourhood = build_neighbourhood(pos);
    size_t n = neighbourhood.size();
    // build eta
    if(m_dim == 1){
        vector eta(build_eta(params, neighbourhood));

        vector result = Eigen::VectorXd::Zero(m_z->cols());
        // compute the mean of z in position pos
        for (size_t i = 0; i < n; ++i) {
            result += eta(i) * m_z->row(neighbourhood[i]);
        }

        return result;
    }
    matrix eta = build_etavec( params, neighbourhood);

    vector result = Eigen::VectorXd::Zero(m_z->cols());
        // compute the mean of z in position pos
    for (size_t i = 0; i < n; ++i) {
        vector row = m_z->row(neighbourhood[i]);
        result += eta.block(i*m_dim,0,m_dim,m_dim).transpose() * row;
    }

    return result;

    
}

template <> cd::matrix Predictor::predict_mean<cd::matrix, cd::matrix>(const cd::matrix& pos) const
{
    matrix result(pos.rows(),m_z->cols());
    #pragma omp parallel for
    for (size_t i = 0; i < pos.rows(); ++i) {
        Eigen::VectorXd row(predict_mean<cd::vector, vector>(pos.row(i)));
        result.row(i) = row;
    }
    return result;
}

template <>
std::pair<vector, matrix> Predictor::predict_z<cd::vector, std::pair<vector, matrix>>(const cd::vector& pos) const
{
    size_t n = m_data->rows();
    // predict the mean of z in pos
    vector m0 = predict_mean<cd::vector, vector>(pos);
    vector result = m0;
    // find the value of the parameters in pos
    cd::vector params = m_smt.smooth_vector(pos);
    // build etakriging and calculate the variance
    if(m_dim == 1){
    std::pair<vector, matrix> fulletakriging(build_etakriging(params, pos));
    vector& etakriging = fulletakriging.first;
    // predict the value of z(pos)
    for (size_t i = 0; i < n; ++i) {
        result += etakriging(i) * (m_z->row(i) - m_means->row(i));
    }
    // return z(pos) and the kriging variance
    return std::make_pair(result, fulletakriging.second);
    }

    std::pair<matrix, matrix> fulletakriging(build_etakrigingvec(params, pos));
    matrix etakriging = fulletakriging.first;
    // predict the value of z(pos)
    for (size_t i = 0; i < n; ++i) {
        Eigen::RowVectorXd row = m_z->row(i) - m_means->row(i);
        result += etakriging.block(i*m_dim, 0, m_dim, m_dim).transpose() * row.transpose();
    }
    // return z(pos) and the kriging variance
    return std::make_pair(result, fulletakriging.second);
}

template <> cd::matrix Predictor::predict_z<cd::matrix, cd::matrix>(const cd::matrix& pos) const
{
    matrix result(pos.rows(),m_z->cols() + 1);
    #pragma omp parallel for
    for (size_t i = 0; i < pos.rows(); ++i) {
        std::pair<vector, matrix> prediction = predict_z<cd::vector, std::pair<vector, matrix>>(pos.row(i));
        Eigen::RowVectorXd row(prediction.first);
        result.block(i, 0, 1, m_z->cols()) = row;
        result(i, m_z->cols()) = prediction.second(0,0);
    }
    return result;
}

cd::matrix Predictor::compute_kriging_matrix_inverse(){

    VariogramFunction& gammaiso = *(m_gammaisoptr);
    size_t n = m_data->rows();
    matrix kriging_matrix(n,n);
    for (size_t i = 0; i < n; ++i) {
        const vector& posi = m_data->row(i);
        cd::vector paramsi = m_smt.smooth_vector(posi);
        for (size_t j = i; j < n; ++j) {
            const vector& posj = m_data->row(j);
            cd::vector paramsj = m_smt.smooth_vector(posj);
            cd::vector s = posi - posj;
            kriging_matrix(i, j) = gammaiso(paramsi, paramsj, s[0], s[1]);
            kriging_matrix(j, i) = gammaiso(paramsi, paramsj, s[0], s[1]);
        }
    }
    return 	kriging_matrix.llt().solve(Eigen::MatrixXd::Identity(n,n));
}

cd::matrix Predictor::compute_kriging_matrix_inverse_vec(){

    size_t n = m_data->rows();
    matrix kriging_matrix(n*m_dim,n*m_dim);
    for (size_t i = 0; i < n; ++i) {
        const vector& posi = m_data->row(i);
        cd::vector paramsi = m_smt.smooth_vector(posi);
        std::vector<vector> paramveci;
        for (size_t i = 0, totparams = 0; i < m_id.size(); i++ ){
            unsigned n_params = 3;
                paramveci.push_back(paramsi.segment(totparams, n_params + m_dim*(m_dim+1)/2));
            totparams += n_params + m_dim*(m_dim+1)/2;
            }

        for (size_t j = i; j < n; ++j) {
            const vector& posj = m_data->row(j);
            cd::vector paramsj = m_smt.smooth_vector(posj);
            std::vector<vector> paramvecj;
            for (size_t i = 0, totparams = 0; i < m_id.size(); i++ ){
            unsigned n_params = 3;
                paramvecj.push_back(paramsj.segment(totparams, n_params + m_dim*(m_dim+1)/2));
            totparams += n_params + m_dim*(m_dim+1)/2;
            }
            cd::vector s = posi - posj;
            kriging_matrix.block(i*m_dim, j*m_dim, m_dim, m_dim) = m_crosscov(paramveci, paramvecj, s[0], s[1]);
            kriging_matrix.block(j*m_dim, i*m_dim, m_dim, m_dim) = kriging_matrix.block(i*m_dim, j*m_dim, m_dim, m_dim).transpose();
        }
    }
    return 	kriging_matrix.llt().solve(Eigen::MatrixXd::Identity(n*m_dim,n*m_dim));
}

Predictor::Predictor(
    const std::vector<std::string>& id, const cd::matrixptr& z, const size_t& dim,const Smt& mysmt, const double& b, const cd::matrixptr& data, const cd::matrixptr& anchorpoints, const bool& predict_y)
    : m_gammaisoptr(make_variogramiso(id[0]))
    , m_crosscov(id, dim, id.size())
    , m_z(z)
    , m_smt(mysmt)
    , m_b(b)
    , m_data(data)
    , m_dim(dim)
    , m_a(anchorpoints)
    , m_id(id)
{
    m_means = std::make_shared<matrix>(z->rows(),z->cols());
    m_anchor_means = std::make_shared<matrix>(m_a->rows(),z->cols());
    // build a vector with the prediction of the mean of z in every anchorpoint to speed up the next computations
    
    #pragma omp parallel for
    for (size_t i = 0; i < m_a->rows(); ++i) {
        Eigen::VectorXd row(predict_mean<size_t, vector>(i));
        m_anchor_means->row(i) = row;
    }
    m_smt.add_mean(m_anchor_means);
    
    *m_means = predict_mean<cd::matrix, cd::matrix>(*m_data);

    if(predict_y){
        if(m_dim == 1)
            m_kriging_matrix_inverse = std::make_shared<matrix>(compute_kriging_matrix_inverse());
        else
            m_kriging_matrix_inverse = std::make_shared<matrix>(compute_kriging_matrix_inverse_vec());
    }

};

Predictor::Predictor()
    : m_gammaisoptr(make_variogramiso("esponenziale")) {}
} // namespace LocallyStationaryModels
