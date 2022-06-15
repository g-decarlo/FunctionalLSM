// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#ifndef LOCALLY_STATIONARY_MODELS_TRAITS
#define LOCALLY_STATIONARY_MODELS_TRAITS

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <omp.h>
#include <string>
#include <vector>

#include "Eigen/Dense"
#include "tolerances.hpp"

namespace LocallyStationaryModels {
/**
 * Namespace cd
 * \brief collects all the typedef used inside the code
 */
namespace cd {
    // defining basic types
    using vector = Eigen::VectorXd;
    using matrix = Eigen::MatrixXd;
    using matrixI = Eigen::MatrixXi;
    using vectorptr = std::shared_ptr<vector>;
    using matrixptr = std::shared_ptr<matrix>;
    using matrixIptr = std::shared_ptr<matrixI>;
    using vectorind = std::vector<size_t>;

    // defining function types
    using kernelfunction = std::function<double(const vector&, const vector&, const double&)>;
    using gridfunction = std::function<matrixIptr(const matrixptr&, const size_t&, const size_t&, const double&)>;

} // namespace cd
} // namespace LocallyStationaryModels

#endif // LOCALLY_STATIONARY_MODELS_TRAITS
