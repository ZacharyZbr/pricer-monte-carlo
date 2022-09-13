#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

/// \brief Classe AsianOption
class AsianOption :: Option
{
    /**
     * Calcule la valeur du payoff sur la trajectoire
     *
     * @param[in] path est une matrice de taille (N+1) x d
     * contenant une trajectoire du modèle telle que créée
     * par la fonction asset.
     * @return phi(trajectoire)
     */

    AsianOption::AsianOption(double T, int nbTimeSteps, int size, float strike, PnlVect coefficients){
        this->strike_= strike;
        this->T_ = T;
        this->nbTimeSteps_ = nbTimeSteps;
        this->coefficients_ = coefficients;
        this->size_ = size;
    }

    double AsianOption::payoff(const PnlMat* path){

        double payoff = 0;
        
        for(int d = 0; d < size; size++){
            double lambda = coefficients->array[d];
            lambda /= (nbTimeSteps + 1);
            double sum = 0;
            for(int i = 0; i <= nbTimeSteps; i++){
                sum += path->array[d][i];
            }
            payoff += lambda * sum;
        }
        if (payoff > strike){
            return payoff - strike;
        }
        else{
            double zero = 0;
            return zero;
        }

    };
};