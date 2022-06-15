// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#ifndef LOCALLY_STATIONARY_MODELS_KERNEL
#define LOCALLY_STATIONARY_MODELS_KERNEL

#include "traits.hpp"

namespace LocallyStationaryModels {
/**
 * \brief class to compute the kernel matrix
 */
class Kernel {
private:
    double m_epsilon; ///< bandwidth parameter
    cd::kernelfunction m_f; ///< kernel function
    cd::matrixptr m_k = std::make_shared<cd::matrix>(0, 0); ///< kernel matrix

public:
    /**
     * \brief constructor
     * \param id name of the kernel function
     * \param epsilon value of the bandwidth parameter epsilon
     */
    Kernel(const std::string& id, const double& epsilon);

    /**
     * \brief default constuctor with a gaussian kernel and epsilon equal to 1.
     */
    Kernel();

    /**
     * \return m_f(x ,y) where m_f is the kernel function
     */
    double operator()(const cd::vector& x, const cd::vector& y) const;

    /**
     * \brief build the "star" version of the kernel that contains the standardized kernel weights in such
     * a way that each row sums to one
     * \param data a shared pointer to the matrix with the coordinates of the original dataset
     * \param anchorpoints a shared pointer to the matrix with the coordinates of the anchor points
     */
    void build_kernel(const cd::matrixptr& data, const cd::matrixptr& anchorpoints);

    /**
     * \brief build the "standard" version of the kernel needed for smoothing
     * \param coordinates a shared pointer to the matrix with the coordinates
     */
    void build_simple_kernel(const cd::matrixptr& coordinates);

    /**
     * \brief build the "standard" version of the kernel needed for smoothing
     * \param coordinates a shared pointer to the matrix with the coordinates
     * \param epsilon replace the old epsilon with a new value
     */
    void build_simple_kernel(const cd::matrixptr& coordinates, const double& epsilon);

    /**
     * \return a shared pointer to the matrix pointed by m_k
     */
    const cd::matrixptr get_kernel() const;
}; // class Kernel
} // namespace LocallyStationaryModels

#endif // LOCALLY_STATIONARY_MODELS_KERNEL
