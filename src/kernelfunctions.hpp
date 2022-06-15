// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#ifndef LOCALLY_STATIONARY_MODELS_KERNEL_FUNCTIONS
#define LOCALLY_STATIONARY_MODELS_KERNEL_FUNCTIONS

#include "traits.hpp"

namespace LocallyStationaryModels {
/**
 * Namespace kf
 * \brief collect the kernel functions and the corresponding factory function
 */
namespace kf {
    /**
     * \return e^(-norm(x-y)^2/epsilon^2)
     */
    double gaussian(const cd::vector& x, const cd::vector& y, const double& epsilon);

    /**
     * \return 1 only if the norm of the difference between x and y is less the the square of epsilon
     */
    double identity(const cd::vector& x, const cd::vector& y, const double& epsilon);

    /**
     * \brief allow to select between pre-built kernel functions
     * \param id a string with the name of the kernel function you want to use
     */
    cd::kernelfunction make_kernel(const std::string& id);
} // namespace kf
} // namespace LocallyStationaryModels

#endif // LOCALLY_STATIONARY_MODELS_KERNEL_FUNCTION
