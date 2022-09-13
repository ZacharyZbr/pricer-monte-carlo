#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "PerformanceOption.hpp"

PerformanceOption::PerformanceOption(double T, int nbTimeSteps, int size, PnlVect *coefficients)
{
    this->T_ = T;
    this->nbTimeSteps_ = nbTimeSteps;
    this->coefficients_ = coefficients;
    this->size_ = size;
}

double PerformanceOption::payoff(const PnlMat *path)
{

    double payoff = 1;
    for (int i = 0; i < nbTimeSteps_; i++)
    {
        double numerator = 0;
        double denominator = 0;
        for (int d = 0; d < size_; d++)
        {
            double lambda = coefficients_->array[d];

            numerator += lambda * pnl_mat_get(path, d, i + 1); // path->array[i + 1 + d * nbTimeSteps_];
            denominator += lambda * pnl_mat_get(path, d, i);   // path->array[i + d * nbTimeSteps_];
        }
        if (numerator > denominator)
        {
            payoff += numerator / denominator - 1;
        }
    }

    return payoff;
};