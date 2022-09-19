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
    double price = 0;
    double stdev = 1.0;
    double fdStep = 0.1;
    double zero = 0.0;

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
    double t = 0.12;

    P->extract("payoff coefficients", payoff_coefficients, size);
    PnlVect *vect_stdev = pnl_vect_create_from_zero(size);
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    BlackScholesModel *blackScholesModel1 = new BlackScholesModel(size, r, correlation, sigma, spot);
    BlackScholesModel *blackScholesModel2 = new BlackScholesModel(size, r, correlation, sigma, spot);
    if (type == "basket")
    {
        P->extract("strike", strike);
        BasketOption *pBasketOption1 = new BasketOption(T, nbTimeStep, size, payoff_coefficients, strike);

        MonteCarlo *monteCarlo1 = new MonteCarlo(blackScholesModel1, pBasketOption1, rng, fdStep, 50000);
        PnlVect *vect_stdev = pnl_vect_create_from_zero(size);
        PnlVect *delta1 = pnl_vect_create(size);

        // monteCarlo1->price(price, stdev);
        // pnl_vect_print(delta1);
        PnlMat *past = pnl_mat_create_from_zero(size, nbTimeStep + 1);
        blackScholesModel2->asset(past, T, nbTimeStep, rng);
        pnl_mat_resize(past, size, floor(t * nbTimeStep) + 1);
        pnl_mat_set_col(past, spot, 0);
        // PnlMat *shift_path = pnl_mat_create(past->n, past->m);
        monteCarlo1->delta(past, t, delta1,vect_stdev);
        //monteCarlo1->delta(delta1, vect_stdev);
        pnl_vect_print(delta1);
        // blackScholesModel2->shiftAsset(shift_path, past, 1, 9, 0, T / nbTimeStep);
        // pnl_mat_print(shift_path);
        //  pnl_mat_resize(past, size, floor(t * nbTimeStep) + 1);
        //  monteCarlo1->price(past, t, price, stdev);
        delete (pBasketOption1);
        delete (monteCarlo1);
    }
    else if (type == "asian")
    {
        P->extract("strike", strike);
        PnlVect *deltaT = pnl_vect_create(size);
        printf("Calcul de l'option Asiatique \n");
        AsianOption *pAsianOption1 = new AsianOption(T, nbTimeStep, size, strike, payoff_coefficients);
        MonteCarlo *monteCarlo1 = new MonteCarlo(blackScholesModel1, pAsianOption1, rng, T / nbTimeStep, 50000);

        //PnlMat *past = pnl_mat_create_from_zero(size, nbTimeStep + 1);
        //blackScholesModel2->asset(past, T, nbTimeStep, rng);

        //pnl_mat_resize(past, size, floor(t * nbTimeStep) + 1);
        monteCarlo1->price(price, stdev);
        //monteCarlo1->price(past, t, price, stdev);
        // printf("Le prix est : %f\n", price);
        // printf("Liste des deltats en t:\n");
        //monteCarlo1->delta(past, t, deltaT, vect_stdev);
        monteCarlo1->delta(deltaT, vect_stdev);
        pnl_vect_print(deltaT);
        delete (pAsianOption1);
        delete (monteCarlo1);
    }
    else
    {

        PerformanceOption *pPerfOption1 = new PerformanceOption(T, nbTimeStep, size, payoff_coefficients);

        MonteCarlo *monteCarlo1 = new MonteCarlo(blackScholesModel1, pPerfOption1, rng, T / nbTimeStep, 50000);

        PnlMat *past = pnl_mat_create_from_zero(size, nbTimeStep + 1);
        blackScholesModel2->asset(past, T, nbTimeStep, rng);
        pnl_mat_resize(past, size, floor(t * nbTimeStep) + 1);
        monteCarlo1->price(past, t, price, stdev);
        PnlVect *deltaT = pnl_vect_create(size);
        monteCarlo1->delta(deltaT, vect_stdev);
        pnl_vect_print(deltaT);
        delete (pPerfOption1);
        delete (monteCarlo1);
    }

    delete (blackScholesModel1);
    delete (P);

    pnl_vect_free(&divid);
    pnl_vect_free(&payoff_coefficients);
    pnl_rng_free(&rng);
    std::cout << "le prix de l'option " << type << " est " << price << std::endl;
    //std::cout << "largeur de l'intervalle " << type << " est " << stdev << std::endl;

    return 0;
}