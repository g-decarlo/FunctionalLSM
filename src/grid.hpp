// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#ifndef LOCALLY_STATIONARY_MODELS_GRID
#define LOCALLY_STATIONARY_MODELS_GRID

#include "traits.hpp"

namespace LocallyStationaryModels {
/**
 * \brief   class to build the grid
 */
class Grid {
private:
    cd::gridfunction m_f; ///< grid function
    cd::matrixIptr m_g = std::make_shared<cd::matrixI>(0, 0); ///< grid matrix
    cd::vectorptr m_normh
        = nullptr; ///< vector with the norm of each cell of the grid (mean of the norm of all the pairs inside)
    cd::vectorptr m_mean_x
        = nullptr; ///< vector with the x of each cell of the grid (mean of the x of all the pairs inside)
    cd::vectorptr m_mean_y
        = nullptr; ///< vector with the y of each cell of the grid (mean of the y of all the pairs inside)
    double m_epsilon; ///< bandwidth parameter inside

    /**
     * \brief a "helper" function to build the vector containing the position of the centers of the cells of the grid.
     * Each pair of coordinates is assigned to a position of the grid.
     * \param data a shared pointer to the matrix of the coordinates
     */
    void build_normh(const cd::matrixptr& data);

public:
    /**
     * \brief constructor
     * \param id name of the grid function
     * \param epsilon the same epsilon regulating the kernel
     */
    Grid(const std::string& id, const double& epsilon);

    /**
     * \brief constructor. Use the Pizza style by default
     */
    Grid();

    /**
     * \brief build the grid
     * \param data a shared pointer to the matrix of the coordinates
     * \param n_angles number of slices of the pizza
     * \param n_intervals number of the pieces for each slice of the pizza
     */
    void build_grid(const cd::matrixptr& data, const size_t& n_angles, const size_t& n_intervals);

    /**
     * \return a shared pointer to the grid
     */
    const cd::matrixIptr get_grid() const;
    /**
     * \return a shared pointer to m_normh
     */
    const cd::vectorptr get_normh() const;
    /**
     * \return a pointer to the vector containing the xs of the centers of the cells of the grid
     */
    const cd::vectorptr get_x() const;
    /**
     * \return a pointer to the vector containing the ys of the centers of the cells of the grid
     */
    const cd::vectorptr get_y() const;
}; // class Grid
} // namespace LocallyStationaryModels

#endif // LOCALLY_STATIONARY_MODELS_GRID
