#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

/// \brief Classe Option abstraite
class PerformanceOption :: Option
{
  public:
    /**
     * Calcule la valeur du payoff sur la trajectoire
     *
     * @param[in] path est une matrice de taille (N+1) x d
     * contenant une trajectoire du modèle telle que créée
     * par la fonction asset.
     * @return phi(trajectoire)
     */
    PerformanceOption::PerformanceOption(double T, int nbTimeSteps, int size, PnlVect coefficients){
        this->T_ = T;
        this->nbTimeSteps_ = nbTimeSteps;
        this->coefficients_ = coefficients;
        this->size_ = size;
    }

    double PerformanceOption::payoff(const PnlMat* path){

        double payoff = 1;
        for(int i = 0; i < nbTimeSteps; i++){
            double numerator = 0;
            double denominator = 0;
            for(int d = 0; d < size; d++){
                double lambda = coefficients->array[d];
                numerator += lambda * path->array[d][i+1];
                denominator += lambda * path->array[d][i];
            }
            if (numerator > denominator){
                payoff += numerator/denominator - 1;
            }
        }

        return payoff;
        

    };

};