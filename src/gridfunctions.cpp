// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#include "gridfunctions.hpp"

namespace LocallyStationaryModels {
using namespace cd;

namespace gf {
    matrixIptr pizza(const matrixptr& data, const size_t& n_angles, const size_t& n_intervals, const double& epsilon)
    {
        double pi = Tolerances::pi;
        // create a square matrix of dimension data->rows()^2 and fill it with -1
        matrixIptr grid(std::make_shared<matrixI>(matrixI::Constant(data->rows(), data->rows(), -1)));
        // b is the radius of locally stationary neighbourhood as function of bandwidth parameter epsilon
        double b = 2 * epsilon;
        double cell_length = b / n_intervals;
        double cell_angle = pi / (n_angles);

        // for every couple of points i and j in data compute the position of the vector (j - i) in the grid and fill grid(i, j)
        // accordingly since grid is symmetric we only need to fill the upper triangular part of the matrix
        #pragma omp parallel for
        for (size_t i = 0; i < data->rows() - 1; ++i) {
            for (size_t j = i + 1; j < data->rows(); ++j) {
                double deltax = data->operator()(j, 0) - data->operator()(i, 0);
                double deltay = data->operator()(j, 1) - data->operator()(i, 1);
                double radius = std::sqrt(deltax * deltax + deltay * deltay);

                if (radius >= b) {
                    grid->operator()(i, j) = -1;
                } else if (deltax != 0) {
                    grid->operator()(i, j) = floor(radius / cell_length)
                        + n_intervals * floor((pi / 2 + std::atan(deltay / deltax)) / cell_angle);
                } else {
                    grid->operator()(i, j) = floor(radius / cell_length);
                }
            }
        }
        return grid;
    }

    gridfunction make_grid(const std::string& id) { return pizza; }
} // namespace gf
} // namespace LocallyStationaryModels