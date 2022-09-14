#include <iostream>
#include <ctime>
//#include "gtest/gtest.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "BasketOption.hpp"
#include "BlackScholesModel.hpp"
#include "MonteCarlo.hpp"

#include "pnl/pnl_matrix.h"
#include "Option.hpp"

int main(){

    float maturity = 3;
    int nbTimeStep = 2;
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

    return 0;
}