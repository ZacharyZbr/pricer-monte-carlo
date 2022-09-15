#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "AsianOption.hpp"

/**
 * Calcule la valeur du payoff sur la trajectoire
 *
 * @param[in] path est une matrice de taille (N+1) x d
 * contenant une trajectoire du modèle telle que créée
 * par la fonction asset.
 * @return phi(trajectoire)
 */

AsianOption::AsianOption(double T, int nbTimeSteps, int size, float strike, PnlVect *coefficients)
{
    this->strike_ = strike;
    this->T_ = T;
    this->nbTimeSteps_ = nbTimeSteps;
    this->coefficients_ = coefficients;
    this->size_ = size;
}

double AsianOption::payoff(const PnlMat *path)
{

    double payoff = 0;

    for (int d = 0; d < size_; d++)
    {
        double lambda = pnl_vect_get(coefficients_, d);
        lambda /= (nbTimeSteps_ + 1);
        double sum = 0;
        for (int i = 0; i <= nbTimeSteps_; i++)
        {

            // TODO
            sum += pnl_mat_get(path, d, i);
            // sum += path->array[i + d * nbTimeSteps_]; // path->array[d][i];
        }
        payoff += lambda * sum;
    }
    if (payoff > strike_)
    {
        return payoff - strike_;
    }
    else
    {
        double zero = 0;
        return zero;
    }
};