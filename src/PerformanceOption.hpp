#pragma once

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"

/// \brief Classe Option abstraite
class PerformanceOption : public Option
{
public:
  PnlVect *coefficients_;
  /**
   * Calcule la valeur du payoff sur la trajectoire
   *
   * @param[in] path est une matrice de taille (N+1) x d
   * contenant une trajectoire du modèle telle que créée
   * par la fonction asset.
   * @return phi(trajectoire)
   */
  PerformanceOption(double T, int nbTimeSteps, int size, PnlVect *coefficients);
  ~PerformanceOption();
  double payoff(const PnlMat *path) override;
};