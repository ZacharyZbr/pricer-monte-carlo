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

    PnlMat *correlationMat = pnl_mat_create_from_scalar(size_, size_, rho_);
    pnl_mat_set_diag(correlationMat, 1, 0);
    pnl_mat_chol(correlationMat);
    
    for (int underlyingAsset = 0; underlyingAsset < size_; underlyingAsset++)
    {

        
        
        double volatility = pnl_vect_get(sigma_, underlyingAsset);
        
        
        pnl_mat_set(path, underlyingAsset, 0, pnl_vect_get(spot_, 0));
        
        PnlVect *row_chol = pnl_vect_create(size_);
        pnl_mat_get_row(row_chol, correlationMat, underlyingAsset);
        
        for (int k = 1; k < nbTimeSteps; k++)
        {
            
            PnlVect *G = pnl_vect_new();
            pnl_vect_rng_normal(G, size_, rng);
            double scalar_product = pnl_vect_scalar_prod(row_chol, G);

            double new_price = pnl_mat_get(path, underlyingAsset, k - 1) * exp((r_ - pow(volatility, 2) / 2) * ((T / nbTimeSteps)) + (volatility * sqrt((T / nbTimeSteps)) * scalar_product));

            pnl_mat_set(path, underlyingAsset, k, new_price);
        }
    }
}