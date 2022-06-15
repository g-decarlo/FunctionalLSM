// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>

#ifndef LOCALLY_STATIONARY_MODELS_KRIGING
#define LOCALLY_STATIONARY_MODELS_KRIGING

#include "smooth.hpp"
#include "traits.hpp"
#include "variogramfit.hpp"

namespace LocallyStationaryModels {
/**
 * \brief class to perform kriging on the data
 */
class Predictor {
private:
    std::shared_ptr<VariogramFunction> m_gammaisoptr; ///< variogram function
    cd::matrixptr m_z = nullptr; ///< z(m_data)
    Smt m_smt; ///< smoother
    double m_b; ///< cutoff-radius of locally stationary neighbourhood
    cd::matrixptr m_means = nullptr; ///< vector with the mean predicted in each anchor point
    cd::matrixptr m_data = nullptr; ///< dataset with the initial points
    cd::matrixptr m_kriging_matrix_inverse = nullptr; ///<kriging covariance matrix in anchor points
    size_t m_dim;
    std::vector<std::string> m_id; ///< names of the chosen variograms
    CrossCovariance m_crosscov;///< cross covariance operator
    cd::matrixptr m_a;

    /**
     * \brief build a vector with the index of the points in the neighbourhood of radius b of the point in position pos
     * \param pos a vector of coordinates or the index of the position of the center of the neighbourhood
     */
    cd::vectorind build_neighbourhood(const cd::vector& pos) const;
    cd::vectorind build_neighbourhood(const size_t& pos) const;

    /**
     * \brief build the vector eta necessary to perform kriging on the mean of Y in a point
     * \param params the params obtained by smoothing in the center of the neighbourhood
     * \param neighbourhood a "neighbourhood" vector build with the previous functions
     */
    cd::vector build_eta(cd::vector& params, cd::vectorind& neighbourhood) const;
    cd::matrix build_etavec( cd::vector& params,  cd::vectorind& neighbourhood) const;
    /**
     * \brief build the vector eta necessary to perform kriging on Y in a point
     * \param params the params obtained by smoothing in the center of the neighbourhood
     * \param neighbourhood a "neighbourhood" vector build with the previous functions
     */
    std::pair<cd::vector, cd::matrix> build_etakriging(const cd::vector& params, const cd::vector& pos) const;
    std::pair<cd::matrix, cd::matrix> build_etakrigingvec(const cd::vector& params, const cd::vector& pos) const;

    cd::matrix compute_kriging_matrix_inverse();
    cd::matrix compute_kriging_matrix_inverse_vec();


public:
    /**
     * \brief constructor
     * \param id name of the variogram function associated with the problem
     * \param z the vector with the value of the function Y in the known points
     * \param mysmt the one used to previously smooth the variogram
     * \param b the radius of the neighbourhood of the point where to perform kriging
     * \param data a shared pointer to the matrix with the coordinates of the original dataset
     */
    Predictor(
        const std::vector<std::string>& id, const cd::matrixptr& z, const size_t& dim, const Smt& mysmt, const double& b, const cd::matrixptr& data, const cd::matrixptr& anchorpoints);
    /**
     * \brief gammaiso set by default to exponential
     */
    Predictor();

    /**
     * \brief predict the mean
     */
    template <typename Input, typename Output> Output predict_mean(const Input& pos) const;

    /**
     * \brief predict Z
     */
    template <typename Input, typename Output> Output predict_z(const Input& pos) const;
}; // class Predictor
} // namespace LocallyStationaryModels

#endif // LOCALLY_STATIONARY_MODELS_KRIGING
