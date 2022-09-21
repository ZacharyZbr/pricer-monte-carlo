#include "Option.hpp"
#include "BlackScholesModel.hpp"
#include "pnl/pnl_random.h"
#include "MonteCarlo.hpp"
#include "assert.h"
#include <omp.h>

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
    PnlMat *pMatrix = pnl_mat_create_from_zero(nb_assets, steps + 1);
    for (long sample = 0; sample < nbSamples_; sample++)
    {
        mod_->asset(pMatrix, opt_->T_, steps, rng_);
        double payoff = opt_->payoff(pMatrix);
        meanPayoff += payoff;
        meanPayoffSquared += payoff * payoff;
    }
    pnl_mat_free(&pMatrix);
    prix = exp(-mod_->r_ * opt_->T_) * meanPayoff / nbSamples_;
    double ksiSquared = exp(-2 * mod_->r_ * opt_->T_) * (meanPayoffSquared / nbSamples_ - ((meanPayoff / nbSamples_) * (meanPayoff / nbSamples_)));
    std_dev = sqrt(ksiSquared / nbSamples_);
}

/**
 * Calcule le prix de l'option à la date 0 avec une boucle parallèle
 *
 * @param[out] prix valeur de l'estimateur Monte Carlo
 * @param[out] ic écart type de l'estimateur
 */
void MonteCarlo::parallelprice(double &prix, double &std_dev)
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
    PnlVect *meanPayoffSquared = pnl_vect_create(std_dev->size);

    for (long sample = 0; sample < nbSamples_; sample++)
    {

        PnlMat *pMatrix = pnl_mat_create_from_zero(nb_assets, steps + 1);
        mod_->asset(pMatrix, t, opt_->T_, steps, rng_, past);
        PnlMat *shiftedMatrixPlus = pnl_mat_create_from_zero(nb_assets, steps + 1);
        PnlMat *shiftedMatrixMinus = pnl_mat_create_from_zero(nb_assets, steps + 1);

        for (int d = 0; d < opt_->size_; d++)
        {

            mod_->shiftAsset(shiftedMatrixPlus, pMatrix, d, fdStep_, t, opt_->T_ / opt_->nbTimeSteps_);

            mod_->shiftAsset(shiftedMatrixMinus, pMatrix, d, -fdStep_, t, opt_->T_ / opt_->nbTimeSteps_);
            delta->array[d] += (opt_->payoff(shiftedMatrixPlus) - opt_->payoff(shiftedMatrixMinus)) / pnl_mat_get(past, d, past->n - 1);
            meanPayoffSquared->array[d] += (opt_->payoff(shiftedMatrixPlus) - opt_->payoff(shiftedMatrixMinus)) * (opt_->payoff(shiftedMatrixPlus) - opt_->payoff(shiftedMatrixMinus));
        }
        pnl_mat_free(&shiftedMatrixMinus);
        pnl_mat_free(&shiftedMatrixPlus);
        pnl_mat_free(&pMatrix);
    }
    for (int d = 0; d < opt_->size_; d++)
    {
        double ksi_carre = (exp(-2 * mod_->r_ * (opt_->T_ - t)) / (0.1 * 0.1 * 2 * 2 * pnl_mat_get(past, d, past->n - 1) * pnl_mat_get(past, d, past->n - 1))) * ((pnl_vect_get(meanPayoffSquared, d) / nbSamples_) - ((delta->array[d] * pnl_mat_get(past, d, past->n - 1)) / nbSamples_) * ((delta->array[d] * pnl_mat_get(past, d, past->n - 1)) / nbSamples_));
        pnl_vect_set(std_dev, d, sqrt(ksi_carre / nbSamples_));
    }

    double facteur_mult = exp(-mod_->r_ * (opt_->T_ - t)) / (nbSamples_ * 2 * 0.1);
    pnl_vect_mult_scalar(delta, facteur_mult);
}
/*
 * @param[in] t date à laquelle le calcul est fait
 * @param[out] delta contient le vecteur de delta
 * @param[out] std_dev contient l'écart type de l'estimateur
 */
void MonteCarlo::delta(PnlVect *delta, PnlVect *std_dev)
{
    int nb_assets = opt_->size_;
    int steps = opt_->nbTimeSteps_;

    PnlVect *meanPayoffSquared = pnl_vect_create_from_zero(std_dev->size);
    PnlMat *shiftedMatrixMinus = pnl_mat_create_from_zero(nb_assets, steps + 1);
    PnlMat *shiftedMatrixPlus = pnl_mat_create_from_zero(nb_assets, steps + 1);
    PnlMat *pMatrix = pnl_mat_create_from_zero(nb_assets, steps + 1);
    for (long sample = 0; sample < nbSamples_; sample++)
    {

        mod_->asset(pMatrix, opt_->T_, steps, rng_);
        for (int d = 0; d < opt_->size_; d++)
        {

            mod_->shiftAsset(shiftedMatrixPlus, pMatrix, d, 0.1, 0, opt_->T_ / opt_->nbTimeSteps_);

            mod_->shiftAsset(shiftedMatrixMinus, pMatrix, d, -0.1, 0, opt_->T_ / opt_->nbTimeSteps_);
            delta->array[d] += (opt_->payoff(shiftedMatrixPlus) - opt_->payoff(shiftedMatrixMinus)) / mod_->spot_->array[d];
            double delta_d = pnl_vect_get(delta, d);
            pnl_vect_set(meanPayoffSquared, d, pnl_vect_get(meanPayoffSquared, d) + delta_d * delta_d);
        }
    }
    pnl_mat_free(&pMatrix);
    pnl_mat_free(&shiftedMatrixMinus);
    pnl_mat_free(&shiftedMatrixPlus);
    for (int d = 0; d < opt_->size_; d++)
    {
        double delta_d = pnl_vect_get(delta, d);
        double ksi_carre = (exp(-2 * mod_->r_ * opt_->T_) / (0.1 * 0.1 * 2 * 2 * mod_->spot_->array[d] * mod_->spot_->array[d])) * ((pnl_vect_get(meanPayoffSquared, d) / nbSamples_) - ((delta->array[d] * mod_->spot_->array[d]) / nbSamples_) * ((delta->array[d] * mod_->spot_->array[d]) / nbSamples_));
        pnl_vect_set(std_dev, d, sqrt(ksi_carre / nbSamples_));
    }

    double facteur_mult = exp(-mod_->r_ * opt_->T_) / (nbSamples_ * 2 * 0.1);
    pnl_vect_mult_scalar(delta, facteur_mult);
    pnl_vect_free(&meanPayoffSquared);
}

void MonteCarlo::paralleldelta(PnlVect *delta, PnlVect *std_dev)
{
    int nb_assets = opt_->size_;
    int steps = opt_->nbTimeSteps_;

    PnlVect *meanPayoffSquared = pnl_vect_create(std_dev->size);
    PnlMat *shiftedMatrixMinus = pnl_mat_create(nb_assets, steps + 1);
    PnlMat *shiftedMatrixPlus = pnl_mat_create(nb_assets, steps + 1);
    PnlMat *pMatrix = pnl_mat_create(nb_assets, steps + 1);

    for (long sample = 0; sample < nbSamples_; sample++)
    {

        mod_->assetP(pMatrix, opt_->T_, steps, rng_);
        int d;

        //#pragma omp parallel shared(pMatrix) private(d)
        //{
        //#pragma omp parallel for
        for (int d = 0; d < opt_->size_; d++)
        {

            mod_->shiftAsset(shiftedMatrixPlus, pMatrix, d, 0.1, 0, opt_->T_ / opt_->nbTimeSteps_);
            mod_->shiftAsset(shiftedMatrixMinus, pMatrix, d, -0.1, 0, opt_->T_ / opt_->nbTimeSteps_);
            delta->array[d] += (opt_->payoff(shiftedMatrixPlus) - opt_->payoff(shiftedMatrixMinus)) / pnl_vect_get(mod_->spot_, d);
            double delta_d = pnl_vect_get(delta, d);
            pnl_vect_set(meanPayoffSquared, d, pnl_vect_get(meanPayoffSquared, d) + delta_d * delta_d);
        }
        //}
    }
    pnl_mat_free(&shiftedMatrixMinus);
    pnl_mat_free(&shiftedMatrixPlus);
    pnl_mat_free(&pMatrix);

    for (int d = 0; d < opt_->size_; d++)
    {
        double delta_d = pnl_vect_get(delta, d);
        double ksi_carre = (1 / (0.1 * 0.1 * 2 * 2)) * (exp(-2 * mod_->r_ * opt_->T_) * 1 / nbSamples_ * pnl_vect_get(meanPayoffSquared, d) - (delta_d / nbSamples_) * (delta_d / nbSamples_) * exp(-2 * mod_->r_ * opt_->T_));
        pnl_vect_set(std_dev, d, sqrt(ksi_carre / nbSamples_));
    }

    double facteur_mult = exp(-mod_->r_ * opt_->T_) / (nbSamples_ * 2 * 0.1);
    pnl_vect_mult_scalar(delta, facteur_mult);
}

void MonteCarlo::PL(const PnlMat *matriceTot, double &PL)
{
    double pas = opt_->T_ / matriceTot->m;
    PnlVect *delta1 = pnl_vect_create(opt_->size_);
    PnlVect *vect_stdDev = pnl_vect_create(opt_->size_);
    double price1;
    double std_dev;
    delta(delta1, vect_stdDev);
    price(price1, std_dev);
    PL = price1 - pnl_vect_scalar_prod(delta1, mod_->spot_);
    for (double k = pas; k < opt_->T_; k += pas)
    {
        PnlVect *deltaMoins = pnl_vect_create(opt_->size_);
        PnlVect *deltaPlus = pnl_vect_create(opt_->size_);
        PnlMat *past = pnl_mat_create_from_zero(opt_->size_, floor(opt_->nbTimeSteps_ + 1));
        pnl_mat_clone(past, matriceTot);
        pnl_mat_resize(past, opt_->size_, floor(k / pas) + 1);
        delta(past, k - pas, deltaMoins, vect_stdDev);
        delta(past, k, deltaPlus, vect_stdDev);
        pnl_vect_minus_vect(deltaPlus, deltaMoins);
        PnlVect *col = pnl_vect_create(opt_->size_);
        pnl_mat_get_col(col, past, past->n - 1);
        PL = PL * exp((mod_->r_ * opt_->T_) / matriceTot->m) - pnl_vect_scalar_prod(deltaPlus, col);
        printf(" PL boucle : %f\n", PL);
    }
}