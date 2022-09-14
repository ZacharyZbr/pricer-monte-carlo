#pragma once

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"

/// \brief Classe Option abstraite
class AsianOption : public Option
{
public:
  float strike_; /// strike
  PnlVect *coefficients_;
  /**
   * Calcule la valeur du payoff sur la trajectoire
   *
   * @param[in] path est une matrice de taille (N+1) x d
   * contenant une trajectoire du modèle telle que créée
   * par la fonction asset.
   * @return phi(trajectoire)
   */
  AsianOption(double T, int nbTimeSteps, int size, float strike, PnlVect *coefficients);

  double payoff(const PnlMat *path) = 0;
};