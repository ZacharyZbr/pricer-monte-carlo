#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "BasketOption.hpp"

/// \brief Classe Option abstraite

BasketOption::BasketOption(double T, int nbTimeSteps, int size, PnlVect *lambda, double strike)
{
  this->T_ = T;
  this->nbTimeSteps_ = nbTimeSteps;
  this->size_ = size;
  this->lambda_ = lambda;
  this->strike_ = strike;
}

/**
 * Calcule la valeur du payoff sur la trajectoire
 *
 * @param[in] path est une matrice de taille (N+1) x d
 * contenant une trajectoire du modèle telle que créée
 * par la fonction asset.
 * @return phi(trajectoire)
 */
double BasketOption::payoff(const PnlMat *path)
{
  double somme = 0;
  for (int k = 0; k <= nbTimeSteps_; k++)
  {
    somme += pnl_mat_get(path, size_ - 1, k) * lambda_->array[k];
  }
  double payoff = somme - strike_;
  if (payoff > 0)
  {
    return payoff;
  }
  else
  {
    return 0;
  }
}
