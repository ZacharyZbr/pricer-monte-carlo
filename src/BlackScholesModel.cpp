#include <iostream>
#include <ctime>
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "BlackScholesModel.hpp"

/// \brief Modèle de Black Scholes

// Default Constructor
// BlackScholesModel() {}

// Initialize Black Scholes Model with values
BlackScholesModel::BlackScholesModel(int size, double r, double rho, PnlVect *sigma, PnlVect *spot)
{
    this->sigma_ = sigma;
    this->r_ = r;
    this->rho_ = rho;
    this->size_ = size;
    this->spot_ = spot;
    this->correlationMat_ = pnl_mat_create_from_scalar(size_, size_, rho_);
    this->gaussian_ = pnl_vect_new();
    this->rowChol_ = pnl_vect_create(size_);
    pnl_mat_set_diag(correlationMat_, 1, 0);
    pnl_mat_chol(correlationMat_);
}

void BlackScholesModel::asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng)
{

    pnl_mat_set_col(path, spot_, 0);

    for (int k = 1; k <= nbTimeSteps; k++)
    {
        pnl_vect_rng_normal(gaussian_, size_, rng);

        for (int underlyingAsset = 0; underlyingAsset < size_; underlyingAsset++)
        {

            pnl_mat_get_row(rowChol_, correlationMat_, underlyingAsset);

            double volatility = pnl_vect_get(sigma_, underlyingAsset);

            double scalar_product = pnl_vect_scalar_prod(rowChol_, gaussian_);

            double new_price = pnl_mat_get(path, underlyingAsset, k - 1) * exp((r_ - (volatility * volatility / 2)) * ((T / nbTimeSteps)) + (volatility * sqrt((T / nbTimeSteps)) * scalar_product));

            pnl_mat_set(path, underlyingAsset, k, new_price);
        }
    }
}