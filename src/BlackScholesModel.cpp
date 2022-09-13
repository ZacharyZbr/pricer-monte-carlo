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
    for (int k = 0; k++; k <= nbTimeSteps)
    {
        printf("k = %f ; \n", k);
        path->array[k] = k * 0.1;
    }
}