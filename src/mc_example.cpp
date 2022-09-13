#include <iostream>
#include <ctime>
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "BlackScholesModel.hpp"
#include "pnl/pnl_matrix.h"

using namespace std;

int main()
{
    PnlVect *G = pnl_vect_new();
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    long M = 1E5;
    int dim = 2;
    pnl_rng_sseed(rng, time(NULL));
    int var_1 = 1;
    double var_2 = 2.0;
    BlackScholesModel blackScholesModel1(var_1, 6, var_2, G, G);
    PnlMat *pMatrix = pnl_mat_create_from_zero(1, 3);
    blackScholesModel1.asset(pMatrix, 6, 3, rng);

    // BlackScholesModel * pBlackScholesModel1 = new BlackScholesModel(var_1, var_2, var_2, G, G);

    double acc = 0.,
           var = 0;

    for (int i = 0; i < M; i++)
    {
        pnl_vect_rng_normal(G, dim, rng);
        double tmp = pnl_vect_norm_two(G);
        acc += tmp;
        var += tmp * tmp;
    }

    acc /= M;
    var = var / M - acc * acc;

    // cout << "E[||G||_2] = " << acc << endl;
    // cout << "Var(||G||_2) = " << var << endl;

    pnl_vect_free(&G);
    pnl_rng_free(&rng);
    // cout << "Black Scholes Model is : " << blackScholesModel1.r_ << " " << blackScholesModel1.rho_ << " " << blackScholesModel1.sigma_ << endl;
    // pnl_mat_print(pMatrix);
    return 0;
}
