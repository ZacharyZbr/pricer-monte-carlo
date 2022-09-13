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
}

void BlackScholesModel::asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng)
{
    printf("Début de la boucle \n");
    for (int underlyingAsset = 0; underlyingAsset < size_ * nbTimeSteps; underlyingAsset += nbTimeSteps)
    {
        printf("asset %d : \n", underlyingAsset);
        for (int k = 0; k < nbTimeSteps; k++)
        {
            double brownian_t = pnl_rng_normal(rng) * pow((k * T / nbTimeSteps), 2);
            double s_0 = spot_->array[underlyingAsset / nbTimeSteps];
            double volatility_0 = sigma_->array[underlyingAsset / nbTimeSteps];
            printf("We put %f value in place %d \n", s_0 * exp((r_ - pow(volatility_0, 2) / 2) + volatility_0 * brownian_t), k + (underlyingAsset));
            path->array[k + (underlyingAsset)] = s_0 * exp((r_ - pow(volatility_0, 2) / 2) + volatility_0 * brownian_t);
        }
    }
}