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
    for (int k = 0; k <= nbTimeSteps; k++)
    {
        printf("k = %f ; \n", k * T / nbTimeSteps);
        // path->array[k] = k * 0.1;
        double brownian_t = pnl_rng_normal(rng) * pow((k * T / nbTimeSteps), 2);
        printf("brownian simulation is : %f \n", brownian_t);
        double s_0 = spot_->array[0];
        double volatility_0 = sigma_->array[0];
        path->array[k] = s_0 * exp((r_ - pow(volatility_0, 2) / 2) + volatility_0 * brownian_t);
    }
}