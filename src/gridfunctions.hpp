// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#ifndef LOCALLY_STATIONARY_MODELS_GRID_FUNCTIONS
#define LOCALLY_STATIONARY_MODELS_GRID_FUNCTIONS

#include "traits.hpp"

namespace LocallyStationaryModels {
/**
 * Namespace gf
 * \brief collect the grid functions and the corresponding factory function
 */
namespace gf {
    /**
     * \brief this function builds a 2D-grid using a "a fette di pizza" (slices-of-pizza like) algorithm to partition
     * the domain
     * \param data a shared pointer to the matrix of the coordinates
     * \param n_angles number of slices of the pizza
     * \param n_intervals number of the pieces for each slice of the pizza
     * \param epsilon bandwidth parameter epsilon. Same of the kernel
     */
    cd::matrixIptr pizza(
        const cd::matrixptr& data, const size_t& n_angles, const size_t& n_intervals, const double& epsilon);

    /**
     * \brief allow to select between the preferred method to build the grid
     * \param id name of the function of choice
     */
    cd::gridfunction make_grid(const std::string& id);
} // namespace gf
} // namespace LocallyStationaryModels

#endif // LOCALLY_STATIONARY_MODELS_GRID_FUNCTIONS
