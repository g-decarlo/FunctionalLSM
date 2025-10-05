// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#include "grid.hpp"
#include "gridfunctions.hpp"

namespace LocallyStationaryModels {
using namespace cd;

void Grid::build_grid(const matrixptr& data, const size_t& n_angles, const size_t& n_intervals)
{
    m_g = m_f(data, n_angles, n_intervals, m_epsilon);
    build_normh(data);
}

Grid::Grid(const std::string& id, const double& epsilon)
    : m_f(gf::make_grid(id))
    , m_epsilon(epsilon) {};

Grid::Grid()
    : Grid("pizza", 1.) {};

const matrixIptr Grid::get_grid() const { return m_g; }

const vectorptr Grid::get_normh() const { return m_normh; }

const vectorptr Grid::get_x() const { return m_mean_x; }

const vectorptr Grid::get_y() const { return m_mean_y; }

void Grid::build_normh(const matrixptr& data)
{
    const matrix& d = *(data);
    // n is the number of rows of grid which is equal to the number of points in d
    size_t n = m_g->rows();
    // max_index is the maximum index assigned to any pair of points in the grid
    size_t max_index = m_g->maxCoeff() + 1;

    // m_mean_x will be the vector with the x of each cell of the grid (mean of the x of all the pairs inside)
    m_mean_x = std::make_shared<vector>(vector::Zero(max_index));
    // m_mean_y will be the vector with the y of each cell of the grid (mean of the y of all the pairs inside)
    m_mean_y = std::make_shared<vector>(vector::Zero(max_index));
    // m_normh will be the vector with the norm of each cell of the grid (mean of the norm of all the pairs inside)
    m_normh = std::make_shared<vector>(vector::Zero(max_index));
    // nn[k] will count how many times we will update the k-th element of m_mean_x, m_mean_y and m_normh
    // which corresponds to the number of vectors which fell in the k-th cell of the grid
    Eigen::VectorXi nn = Eigen::VectorXi::Zero(max_index);
    // k is a variable used to index elements according to their position in the grid
    // it has to be an integer since we have initialized each component of the grid with -1
    int k = 0;

    // for every couple of index i and j update m_mean_x, m_mean_y and m_normh in poistion k if m_g->operator(i, j) == k
    for (size_t i = 0; i < n - 1; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            k = m_g->operator()(i, j);

            if (k >= 0) {
                m_normh->operator()(k)
                    += sqrt((d(j, 0) - d(i, 0)) * (d(j, 0) - d(i, 0)) + (d(j, 1) - d(i, 1)) * (d(j, 1) - d(i, 1)));
                // because of the way we have constructed the grid we need to add the absolute value of pairs in the
                // first and third quadrant and subtract in the second and fourth ones
                if ((d(j, 0) - d(i, 0)) * (d(j, 1) - d(i, 1)) < 0) {
                    m_mean_x->operator()(k) -= std::abs(d(j, 0) - d(i, 0));
                } else {
                    m_mean_x->operator()(k) += std::abs(d(j, 0) - d(i, 0));
                }
                m_mean_y->operator()(k) += std::abs(d(j, 1) - d(i, 1));
                nn[k]++;
            }
        }
    }
    // now divide element by element m_normh, m_mean_x and m_mean_y by nn to get the sample mean
    for (size_t u = 0; u < max_index; ++u) {
        if (nn[u] != 0) {
            m_normh->operator[](u) /= nn[u];
            m_mean_x->operator[](u) /= nn[u];
            m_mean_y->operator[](u) /= nn[u];
        }
    }
}
} // namespace LocallyStationaryModels
