#include <iostream>
#include <ctime>
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "BasketOption.hpp"
#include "BlackScholesModel.hpp"
#include "MonteCarlo.hpp"

#include "pnl/pnl_matrix.h"
#include "Option.hpp"

using namespace std;

int main()
{
    float maturity = 3;
    int nbTimeStep = 3;
    int size_option = 40;
    double strike = 100;
    PnlVect *Lambda = pnl_vect_create_from_scalar(40, 0.025);
    BasketOption *pBasketOption1 = new BasketOption(maturity, nbTimeStep, size_option, Lambda, strike);

    PnlVect *G = pnl_vect_new();
    PnlVect *Sigma = pnl_vect_create_from_scalar(40, 0.2);
    PnlVect *Spot = pnl_vect_create_from_scalar(40, 100);
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    long M = 1E5;
    int dim = 40;
    pnl_rng_sseed(rng, time(NULL));
    int size = 40;
    double interest_rate = 0.04879;
    double rho = 0.0;
    BlackScholesModel *blackScholesModel1 = new BlackScholesModel(size, interest_rate, rho, Sigma, Spot);

    MonteCarlo *monteCarlo1 = new MonteCarlo(blackScholesModel1, pBasketOption1, rng, 0.0, 50000);
    double price = 0.0;
    double stdev = 0.0;
    monteCarlo1->price(price, stdev);
    printf("Le prix de l'option est : %f \n", price);
    // PnlMat *pMatrix = pnl_mat_create_from_zero(2, 3);
    // blackScholesModel1.asset(pMatrix, 6, 3, rng);

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

    cout << "E[||G||_2] = " << acc << endl;
    cout << "Var(||G||_2) = " << var << endl;

    pnl_vect_free(&G);
    pnl_rng_free(&rng);
    // cout << "Black Scholes Model is : " << blackScholesModel1.r_ << " " << blackScholesModel1.rho_ << " " << blackScholesModel1.sigma_ << endl;
    // pnl_mat_print(pMatrix);

    // test de la basket option
    //  float maturity = 3;
    //  int nbTimeStep = 1;
    //  int size_option = 40;
    //  PnlVect *Lambda = pnl_vect_create_from_scalar(1, 0.025);
    //  BasketOption BasketOption(maturity, nbTimeStep, size_option, Lambda);
    return 0;
}
