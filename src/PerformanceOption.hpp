#pragma once

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

/// \brief Classe Option abstraite
class PerformanceOption :: Option
{
  public:
    PnlVect* coefficients;
    /**
     * Calcule la valeur du payoff sur la trajectoire
     *
     * @param[in] path est une matrice de taille (N+1) x d
     * contenant une trajectoire du modèle telle que créée
     * par la fonction asset.
     * @return phi(trajectoire)
     */
    PerformanceOption(double T, int nbTimeSteps, int size, PnlVect coefficients);
    double payoff(const PnlMat* path) = 0;
};