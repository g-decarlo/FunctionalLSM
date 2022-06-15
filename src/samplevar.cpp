// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#include "samplevar.hpp"

namespace LocallyStationaryModels {
using namespace cd;

void SampleVar::build_samplevar(const cd::matrixptr& data, const cd::matrixptr& anchorpoints, const cd::matrixptr& z)
{
    m_grid.build_grid(data, m_n_angles, m_n_intervals);

    m_kernel.build_kernel(data, anchorpoints);

    // d is the matrix with the coordinates of the initial points
    const matrix& d = *(data);
    // z is the vector with z(d)
    const matrix& zz = *(z);
    // a is the matrix with the coordinates of the anchor points
    const matrix& a = *(anchorpoints);
    const matrixIptr g = m_grid.get_grid();
    const matrix& K = *(m_kernel.get_kernel());

    size_t n = g->rows();

    size_t max_index = g->maxCoeff();

    size_t N = a.rows();

    m_variogram = std::make_shared<matrix>(matrix::Zero((max_index + 1)*m_dim, N*m_dim));
    m_denominators = std::make_shared<matrix>(matrix::Zero(max_index + 1, N));

    // Z is the matrix that contains the squared norm of the difference between each possible pair z_i and z_j
    matrix Z(n*m_dim, n*m_dim);
    if (m_dim == 1){//if trace/univariate
        #pragma omp parallel for
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                Z(i,j) = (zz.row(i) - zz.row(j)).squaredNorm();
                Z(j, i) = Z(i, j);
            }
        }
    }
    else {////if cross/multivariate
        #pragma omp parallel for
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                Z.block(i*m_dim, j*m_dim, m_dim, m_dim) = (zz.row(i) - zz.row(j)).transpose()*(zz.row(i) - zz.row(j));
                Z.block(j*m_dim, i*m_dim, m_dim, m_dim) = Z.block(i*m_dim, j*m_dim, m_dim, m_dim);
            }
        }
    }

    #pragma omp parallel
    {   
        int k = 0;
        // for every location in d
        #pragma omp for
        for (size_t l = 0; l < N; ++l) {
            Eigen::VectorXi counters = Eigen::VectorXi::Zero(max_index + 1);
            // for every couple of locations in d
            for (size_t i = 0; i < n - 1; ++i) {
                for (size_t j = i + 1; j < n; ++j) {
                    // if the vector between i and j belongs to the cell k
                    k = g->operator()(i, j);
                    if (k >= 0) {
                        double prodotto = K(l, i) * K(l, j);
                        m_variogram->block(k*m_dim, l*m_dim, m_dim, m_dim) += prodotto * Z.block(i*m_dim, j*m_dim, m_dim, m_dim);
                        m_denominators->operator()(k, l) += prodotto;
                        counters[k]++;
                    }
                }
            }
            for (size_t u = 0; u < max_index + 1; ++u) {
                if (counters[u] != 0) {
                    m_variogram->block(u*m_dim, l*m_dim, m_dim, m_dim) /= (2 * m_denominators->operator()(u, l));
                }
            }
        }
    }
    build_squaredweights();
}

void SampleVar::build_squaredweights()
{
    const matrixIptr g = m_grid.get_grid();
    const vectorptr normh = m_grid.get_normh();

    size_t N = m_denominators->cols();
    size_t htot = normh->rows();

    m_squaredweights = std::make_shared<matrix>(matrix::Zero(N, htot));

    #pragma omp parallel for
    for (size_t k = 0; k < N; ++k) {
        for (size_t h = 0; h < htot; ++h) {
            if (normh->operator()(h) != 0) {
                m_squaredweights->operator()(k, h) = m_denominators->operator()(h, k) / normh->operator()(h);
            }
        }
    }
}

SampleVar::SampleVar(
    const std::string& kernel_id, const size_t& n_angles, const size_t& n_intervals, const double& epsilon, const size_t& dim)
    : m_kernel(kernel_id, epsilon)
    , m_grid("pizza", epsilon)
    , m_n_angles(n_angles)
    , m_n_intervals(n_intervals)
    , m_dim{dim} {};

SampleVar::SampleVar()
    : m_kernel()
    , m_grid() {};

const matrixptr SampleVar::get_variogram() const { return m_variogram; }

const matrixptr SampleVar::get_denominators() const { return m_denominators; }

const matrixptr SampleVar::get_squaredweights() const { return m_squaredweights; }

const vectorptr SampleVar::get_x() const { return m_grid.get_x(); }

const vectorptr SampleVar::get_y() const { return m_grid.get_y(); }

const matrixptr SampleVar::get_kernel() const { return m_kernel.get_kernel(); }

const matrixIptr SampleVar::get_grid() const { return m_grid.get_grid(); }

const vectorptr SampleVar::get_normh() const { return m_grid.get_normh(); }
} // namespace LocallyStationaryModels
