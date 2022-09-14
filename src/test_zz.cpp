#include <iostream>
#include <ctime>
//#include "gtest/gtest.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "BasketOption.hpp"
#include "BlackScholesModel.hpp"
#include "MonteCarlo.hpp"
#include "jlparser/parser.hpp"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"

int main(int argc, char** argv){

    double T, r, strike, correlation;
    PnlVect *spot, *sigma, *divid, *payoff_coefficients;
    int size, nbTimeStep;
    size_t n_samples;

    char* infile = argv[1];
    Param* P = new Parser(infile);

    //P->extract("option type", type);
    P->extract("maturity", T);
    P->extract("option size", size);
    P->extract("spot", spot, size);
    P->extract("volatility", sigma, size);
    P->extract("interest rate", r);
    if (P->extract("dividend rate", divid, size, true) == false) {
    divid = pnl_vect_create_from_zero(size);
    }
    P->extract("strike", strike);
    P->extract("sample number", n_samples);
    P->extract("correlation", correlation);
    P->extract("timestep number", nbTimeStep);

    printf("nbTimeStep = %d \n ", nbTimeStep);
    P->extract("payoff coefficients", payoff_coefficients, size); 


    BasketOption *pBasketOption1 = new BasketOption(T, nbTimeStep, size, payoff_coefficients, strike);

    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    // long M = 1E5;
    // int dim = 40;
    pnl_rng_sseed(rng, time(NULL));

    BlackScholesModel *blackScholesModel1 = new BlackScholesModel(size, r, correlation, sigma, spot);

    MonteCarlo *monteCarlo1 = new MonteCarlo(blackScholesModel1, pBasketOption1, rng, T/nbTimeStep, 50000);
    double price = 0.0;
    double stdev = 0.0;
    monteCarlo1->price(price, stdev);
    printf("Le prix de l'option est : %f \n", price);

    return 0;
}