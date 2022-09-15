#include <iostream>
#include <ctime>
#include <string>
//#include "gtest/gtest.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "BasketOption.hpp"
#include "AsianOption.hpp"
#include "PerformanceOption.hpp"
#include "BlackScholesModel.hpp"
#include "MonteCarlo.hpp"
#include "jlparser/parser.hpp"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"

int main(int argc, char **argv)
{

    double T, r, strike, correlation;
    PnlVect *spot, *sigma, *divid, *payoff_coefficients;
    int size, nbTimeStep;
    size_t n_samples;
    std::string type;
    double price = 0.0;
    double stdev = 0.0;

    char *infile = argv[1];
    Param *P = new Parser(infile);

    P->extract("option type", type);
    P->extract("maturity", T);
    P->extract("option size", size);
    P->extract("spot", spot, size);
    P->extract("volatility", sigma, size);
    P->extract("interest rate", r);
    if (P->extract("dividend rate", divid, size, true) == false)
    {
        divid = pnl_vect_create_from_zero(size);
    }
    P->extract("sample number", n_samples);
    P->extract("correlation", correlation);
    P->extract("timestep number", nbTimeStep);

    printf("nbTimeStep = %d \n ", nbTimeStep);
    P->extract("payoff coefficients", payoff_coefficients, size);

    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    BlackScholesModel *blackScholesModel1 = new BlackScholesModel(size, r, correlation, sigma, spot);
    if (type == "basket")
    {
        P->extract("strike", strike);
        BasketOption *pBasketOption1 = new BasketOption(T, nbTimeStep, size, payoff_coefficients, strike);

        MonteCarlo *monteCarlo1 = new MonteCarlo(blackScholesModel1, pBasketOption1, rng, T / nbTimeStep, 50000);
        monteCarlo1->price(price, stdev);
        delete (pBasketOption1);
        delete (monteCarlo1);
    }
    else if (type == "asian")
    {
        P->extract("strike", strike);
        AsianOption *pAsianOption1 = new AsianOption(T, nbTimeStep, size, strike, payoff_coefficients);
        MonteCarlo *monteCarlo1 = new MonteCarlo(blackScholesModel1, pAsianOption1, rng, T / nbTimeStep, 50000);
        monteCarlo1->price(price, stdev);
        delete (pAsianOption1);
        delete (monteCarlo1);
    }
    else
    {

        PerformanceOption *pPerfOption1 = new PerformanceOption(T, nbTimeStep, size, payoff_coefficients);

        MonteCarlo *monteCarlo1 = new MonteCarlo(blackScholesModel1, pPerfOption1, rng, T / nbTimeStep, 50000);
        monteCarlo1->price(price, stdev);
        delete (pPerfOption1);
        delete (monteCarlo1);
    }

    delete (blackScholesModel1);
    delete (P);

    pnl_vect_free(&divid);
    pnl_vect_free(&payoff_coefficients);
    pnl_rng_free(&rng);
    std::cout << "le prix de l'option " << type << " est " << price << std::endl;
    std::cout << "largeur de l'intervalle " << type << " est " << stdev << std::endl;

    return 0;
}