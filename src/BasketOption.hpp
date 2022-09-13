#pragma once

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"

/// \brief Classe Option abstraite
class BasketOption : public Option 
{
  public:
    //double T_;        /// maturité
    //int nbTimeSteps_; /// nombre de pas de temps de discrétisation
    //int size_;        /// dimension du modèle, redondant avec BlackScholesModel::size_
    PnlVect *lambda_;

    BasketOption(double T_, int nbTimeSteps_, int size_, PnlVect *lambda_);
    /**
     * Calcule la valeur du payoff sur la trajectoire
     *
     * @param[in] path est une matrice de taille (N+1) x d
     * contenant une trajectoire du modèle telle que créée
     * par la fonction asset.
     * @return phi(trajectoire)
     */
    virtual double payoff(const PnlMat* path, double K);
};