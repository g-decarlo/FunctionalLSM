// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#include "kernel.hpp"
#include "kernelfunctions.hpp"

namespace LocallyStationaryModels {
using namespace cd;

Kernel::Kernel(const std::string& id, const double& epsilon)
    : m_epsilon(epsilon)
    , m_f(kf::make_kernel(id)) {};

Kernel::Kernel()
    : Kernel("Gaussian", 1.) {};

double Kernel::operator()(const vector& x, const vector& y) const { return m_f(x, y, m_epsilon); }

void Kernel::build_kernel(const matrixptr& data, const matrixptr& anchorpoints)
{
    size_t n = data->rows();

    size_t N = anchorpoints->rows();

    m_k->resize(N, n);

    // fill each component of m_k with the value of the kernel function evaluated between the i-th anchor point and the j-th
    // initial point
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < n; ++j) {
            m_k->operator()(i, j) = this->operator()(anchorpoints->row(i), data->row(j));
        }
    }

    // create a vector with the sum on each row of m_k
    vector sums(N);
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        sums(i) = (m_k->row(i)).sum();
    }

    // divide each element of m_k by the sum of the elements of its row to obtained the normalized version of the kernel
    // matrix K*
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < n; ++j) {
            m_k->operator()(i, j) /= sums(i);
        }
    }
}

void Kernel::build_simple_kernel(const matrixptr& coordinates)
{
    size_t n = coordinates->rows();
    m_k->resize(n, n);
    // fill each component of m_k with the kernel function evaluated between the i-th and the j-th point of d
    #pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            m_k->operator()(i, j) = this->operator()(coordinates->row(i), coordinates->row(j));
            m_k->operator()(j, i) = m_k->operator()(i, j);
        }
    }
}

void Kernel::build_simple_kernel(const matrixptr& coordinates, const double& epsilon)
{
    m_epsilon = epsilon;
    build_simple_kernel(coordinates);
}

const matrixptr Kernel::get_kernel() const { return m_k; }
} // namespace LocallyStationaryModels
