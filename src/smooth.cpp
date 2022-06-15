// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#include "smooth.hpp"

namespace LocallyStationaryModels {
using namespace cd;

double Smt::smooth_value(const size_t& pos, const size_t& n) const
{
    const matrix& K = *(m_kernel.get_kernel());

    double numerator = 0;
    double denominator = 0;

    for (size_t i = 0; i < m_anchorpos->rows(); ++i) {
        numerator += K(pos, i) * m_solutions->operator()(i, n);
        denominator += K(pos, i);
    }
    if (denominator < std::numeric_limits<double>::min()) {
        return 0;
    }
    return numerator / denominator;
}

double Smt::smooth_value(const cd::vector& pos, const size_t& n) const
{
    double numerator = 0;
    double denominator = 0;

    for (size_t i = 0; i < m_anchorpos->rows(); ++i) {
        numerator += m_kernel(pos, m_anchorpos->row(i)) * m_solutions->operator()(i, n);
        denominator += m_kernel(pos, m_anchorpos->row(i));
    }
    if (denominator < std::numeric_limits<double>::min()) {
        return 0;
    }
    return numerator / denominator;
}

Smt::Smt(const cd::matrixptr& solutions, const matrixptr& anchorpos, const double& min_delta, const double& max_delta,
    const std::string& kernel_id)
    : m_anchorpos(anchorpos)
    , m_solutions(solutions)
    , m_kernel(kernel_id, min_delta)
{
    double min_error = std::numeric_limits<double>::infinity();
    m_optimal_delta = (max_delta - min_delta) / 2;
    const size_t n_deltas = Tolerances::n_deltas;
    // find the optimal value of delta via cross-validation
    for (size_t i = 0; i <= n_deltas; i++) {
        double delta = min_delta + i * (max_delta - min_delta) / n_deltas;
        // build a new kernel with bandwidth parameter equal to delta
        m_kernel.build_simple_kernel(anchorpos, delta);
        const matrix& Kk = *(m_kernel.get_kernel());

        double error = 0;
        // find the value of the error function for the current value of delta
        #pragma omp parallel for reduction(+ : error)
        for (size_t j = 0; j < m_anchorpos->rows(); ++j) {
            cd::vector Kkrow = Kk.row(j);
            double predicted_value = smooth_value(j, 3);
            double real_value = m_solutions->operator()(j, 3);
            double weightk2 = (1 - Kk(j, j) / Kkrow.sum()) * (1 - Kk(j, j) / Kkrow.sum());

            error += (real_value - predicted_value) * (real_value - predicted_value) / weightk2;
        }
        if (error < min_error) {
            m_optimal_delta = delta;
            min_error = error;
        }
    }
    // build the final kernel with the optimal value of delta
    m_kernel.build_simple_kernel(m_anchorpos, m_optimal_delta);
}

Smt::Smt(const cd::matrixptr& solutions, const matrixptr& anchorpos, const double delta, const std::string& kernel_id)
    : m_anchorpos(anchorpos)
    , m_solutions(solutions)
    , m_kernel(kernel_id, delta)
    , m_optimal_delta(delta)
{
    m_kernel.build_simple_kernel(m_anchorpos);
}

const cd::matrixptr Smt::get_solutions() const { return m_solutions; }

double Smt::get_optimal_delta() const { return m_optimal_delta; }

const cd::matrixptr Smt::get_anchorpos() const { return m_anchorpos; }

void Smt::add_mean(const cd::matrixptr& means){
    size_t columns = means->cols();
    m_solutions->conservativeResize(m_solutions->rows(), m_solutions->cols()+columns);
    m_solutions->rightCols(columns) = *means;
    return;
}

Smt::Smt()
    : m_kernel() {};
} // namespace LocallyStationaryModels
