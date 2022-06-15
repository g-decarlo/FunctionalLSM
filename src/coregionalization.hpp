// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#ifndef LOCALLY_STATIONARY_MODELS_COREGIONALIZATION
#define LOCALLY_STATIONARY_MODELS_COREGIONALIZATION

#include "variogramfunctions.hpp"

namespace LocallyStationaryModels {
class CrossCovariance {

    size_t m_r;///< Coregionalization dimension
    size_t m_dim;///< Vector dimension
    std::vector<std::shared_ptr<VariogramFunction>> vario_functions;

    cd::matrix matricize(const cd::vector& vec) const;

public:
    CrossCovariance() = default;
    CrossCovariance(const std::vector<std::string>& vario_ids, const size_t& dim, const size_t& r);
    /**
     * \brief return locally stationary cross-variogram operator between 2 locations given their spatial vector difference and
     *  parameters the center of stationarity 
     */
    cd::matrix operator()(const std::vector<cd::vector>& params, const double& x, const double& y) const;
    /**
     * \brief return non-stationary cross-covariance operator between 2 locations given their spatial vector difference and
     *  parameters in such locations
     */
    cd::matrix operator()(const std::vector<cd::vector>& params1, const  std::vector<cd::vector>& params2, const double& x, const double& y) const;
}; // class CrossCovariance

} // namespace LocallyStationaryModels

#endif // LOCALLY_STATIONARY_MODELS_COREGIONALIZATION
