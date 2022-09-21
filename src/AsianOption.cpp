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

AsianOption::~AsianOption() {}

double AsianOption::payoff(const PnlMat *path)
{

    double payoff = 0;
    char c = 'c';
    PnlMat *path_clone = pnl_mat_create(path->n - 1, path->m -1);
    pnl_mat_clone(path_clone, path);
    pnl_mat_cumsum(path_clone, c);
    PnlVect *col = pnl_vect_create(size_);
    pnl_mat_get_col(col, path_clone, path_clone->n - 1);
    double coeff = (nbTimeSteps_ + 1);
    payoff = pnl_vect_scalar_prod(col , coefficients_) / coeff;
    pnl_mat_free(&path_clone);
    pnl_vect_free(&col);
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