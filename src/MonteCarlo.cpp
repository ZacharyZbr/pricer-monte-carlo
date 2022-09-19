#include "Option.hpp"
#include "BlackScholesModel.hpp"
#include "pnl/pnl_random.h"
#include "MonteCarlo.hpp"
#include "assert.h"

MonteCarlo::MonteCarlo(BlackScholesModel *mod, Option *opt, PnlRng *rng, double fdStep, long nbSamples)
{
    this->mod_ = mod;
    this->opt_ = opt;
    this->rng_ = rng;
    this->fdStep_ = fdStep;
    this->nbSamples_ = nbSamples;
}

MonteCarlo::~MonteCarlo() {}

/**
 * Calcule le prix de l'option à la date 0
 *
 * @param[out] prix valeur de l'estimateur Monte Carlo
 * @param[out] ic écart type de l'estimateur
 */
void MonteCarlo::price(double &prix, double &std_dev)
{
    double meanPayoff = 0;
    double meanPayoffSquared = 0;
    int nb_assets = opt_->size_;
    int steps = opt_->nbTimeSteps_;

    for (long sample = 0; sample < nbSamples_; sample++)
    {

        PnlMat *pMatrix = pnl_mat_create_from_zero(nb_assets, steps + 1);
        mod_->asset(pMatrix, opt_->T_, steps, rng_);

        meanPayoff += opt_->payoff(pMatrix);
        meanPayoffSquared += opt_->payoff(pMatrix) * opt_->payoff(pMatrix);
        pnl_mat_free(&pMatrix);
    }
    prix = exp(-mod_->r_ * opt_->T_) * meanPayoff / nbSamples_;
    double ksiSquared = exp(-2 * mod_->r_ * opt_->T_) * (meanPayoffSquared / nbSamples_ - ((meanPayoff / nbSamples_) * (meanPayoff / nbSamples_)));
    std_dev = sqrt(ksiSquared / nbSamples_);
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
void MonteCarlo::price(const PnlMat *past, double t, double &prix, double &std_dev)
{
    double meanPayoff = 0;
    double meanPayoffSquared = 0;
    int nb_assets = opt_->size_;
    int steps = opt_->nbTimeSteps_;

    for (long sample = 0; sample < nbSamples_; sample++)
    {

        PnlMat *pMatrix = pnl_mat_create_from_zero(nb_assets, steps + 1);
        mod_->asset(pMatrix, t, opt_->T_, steps, rng_, past);

        meanPayoff += opt_->payoff(pMatrix);
        meanPayoffSquared += opt_->payoff(pMatrix) * opt_->payoff(pMatrix);
        pnl_mat_free(&pMatrix); 
    }
    prix = exp(-mod_->r_ * opt_->T_) * meanPayoff / nbSamples_;
    double ksiSquared = exp(-2 * mod_->r_ * opt_->T_) * (meanPayoffSquared / nbSamples_ - ((meanPayoff / nbSamples_) * (meanPayoff / nbSamples_)));
    std_dev = sqrt(ksiSquared / nbSamples_);
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
void MonteCarlo::delta(const PnlMat *past, double t, PnlVect *delta, PnlVect *std_dev)
{
    int nb_assets = opt_->size_;
    int steps = opt_->nbTimeSteps_;

    for (long sample = 0; sample < nbSamples_; sample++)
    {

        PnlMat *pMatrix = pnl_mat_create_from_zero(nb_assets, steps + 1);
        mod_->asset(pMatrix, t, opt_->T_, steps, rng_, past);

        for (int d = 0; d < opt_->size_; d++){
            PnlMat *shiftedMatrixPlus = pnl_mat_create_from_zero(nb_assets, steps + 1);
            mod_->shiftAsset(shiftedMatrixPlus, pMatrix, d, fdStep_, t, opt_->T_ / opt_->nbTimeSteps_);
            PnlMat *shiftedMatrixMinus = pnl_mat_create_from_zero(nb_assets, steps + 1);
            mod_->shiftAsset(shiftedMatrixMinus, pMatrix, d, -fdStep_, t, opt_->T_ / opt_->nbTimeSteps_);
            delta->array[d] += (opt_->payoff(shiftedMatrixPlus) - opt_->payoff(shiftedMatrixMinus)) / pnl_mat_get(past, d, past->n - 1);
            pnl_mat_free(&shiftedMatrixMinus);
            pnl_mat_free(&shiftedMatrixPlus);
        }

        pnl_mat_free(&pMatrix);
    }
    double facteur_mult = exp(-mod_->r_ * (opt_->T_ - t)) / (nbSamples_ * 2 * 0.1);
    pnl_vect_mult_scalar(delta, facteur_mult);
}

/**
 * Calcule le delta de l'option à la date 0
 *
 * @param[in] t date à laquelle le calcul est fait
 * @param[out] delta contient le vecteur de delta
 * @param[out] std_dev contient l'écart type de l'estimateur
 */
void MonteCarlo::delta(PnlVect *delta, PnlVect *std_dev)
{
    int nb_assets = opt_->size_;
    int steps = opt_->nbTimeSteps_;

    for (long sample = 0; sample < nbSamples_; sample++)
    {
        PnlMat *pMatrix = pnl_mat_create_from_zero(nb_assets, steps + 1);
        mod_->asset(pMatrix, opt_->T_, steps, rng_);

        for (int d = 0; d < opt_->size_; d++)
        {
            PnlMat *shiftedMatrixPlus = pnl_mat_create_from_zero(nb_assets, steps + 1);
            mod_->shiftAsset(shiftedMatrixPlus, pMatrix, d, 0.1, 0, opt_->T_ / opt_->nbTimeSteps_);
            PnlMat *shiftedMatrixMinus = pnl_mat_create_from_zero(nb_assets, steps + 1);
            mod_->shiftAsset(shiftedMatrixMinus, pMatrix, d, -0.1, 0, opt_->T_ / opt_->nbTimeSteps_);
            delta->array[d] += (opt_->payoff(shiftedMatrixPlus) - opt_->payoff(shiftedMatrixMinus)) / mod_->spot_->array[d];
            pnl_mat_free(&shiftedMatrixMinus);
            pnl_mat_free(&shiftedMatrixPlus);
        }
        pnl_mat_free(&pMatrix);
    }
    double facteur_mult = exp(-mod_->r_ * opt_->T_) / (nbSamples_ * 2 * 0.1);
    pnl_vect_mult_scalar(delta, facteur_mult);
}