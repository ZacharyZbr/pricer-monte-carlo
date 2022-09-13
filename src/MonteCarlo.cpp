#pragma once

#include "Option.hpp"
#include "BlackScholesModel.hpp"
#include "pnl/pnl_random.h"
#include "MonteCarlo.hpp"

MonteCarlo::MonteCarlo(BlackScholesModel *mod, Option *opt, PnlRng *rng, double fdStep, long nbSamples)
{
    this->mod_ = mod;
    this->opt_ = opt;
    this->rng_ = rng;
    this->fdStep_ = fdStep;
    this->nbSamples_ = nbSamples;
}

/**
 * Calcule le prix de l'option à la date 0
 *
 * @param[out] prix valeur de l'estimateur Monte Carlo
 * @param[out] ic écart type de l'estimateur
 */
void MonteCarlo::price(double &prix, double &std_dev)
{
    double meanPayoff = 0;
    int nb_assets = opt_->size_;
    int steps = opt_->nbTimeSteps_;
    PnlMat *pMatrix = pnl_mat_create_from_zero(nb_assets, steps);
    for (long sample = 0; sample <= nbSamples_; sample++)
    {
        mod_->asset(pMatrix, opt_->T_, steps, rng_);
        meanPayoff += opt_->payoff(pMatrix);
    }
    prix = meanPayoff / nbSamples_;
}

/**
 * Calcule le prix de l'option à la date t
 *
 * @param[in]  past contient la trajectoire du sous-jacent
 * jusqu'à l'instant t
 * @param[in] t date à laquelle le calcul est fait
 * @param[out] prix contient le prix
 * @param[out] std_dev contient l'écart type de l'estimateur
 */
void price(const PnlMat *past, double t, double &prix, double &std_dev)
{
    // TODO
}

/**
 * Calcule le delta de l'option à la date t
 *
 * @param[in] past contient la trajectoire du sous-jacent
 * jusqu'à l'instant t
 * @param[in] t date à laquelle le calcul est fait
 * @param[out] delta contient le vecteur de delta
 * @param[out] std_dev contient l'écart type de l'estimateur
 */
void delta(const PnlMat *past, double t, PnlVect *delta, PnlVect *std_dev)
{
    // TODO
}

/**
 * Calcule le delta de l'option à la date 0
 *
 * @param[in] t date à laquelle le calcul est fait
 * @param[out] delta contient le vecteur de delta
 * @param[out] std_dev contient l'écart type de l'estimateur
 */
void delta(PnlVect *delta, PnlVect *std_dev)
{
    // TODO
}