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
#include <chrono>
#include "assert.h"
#include "PricingResults.hpp"

using namespace std;
int main(int argc, char **argv)
{

    double T, r, strike, correlation;
    PnlVect *spot, *spot2, *sigma, *sigma2, *divid, *payoff_coefficients;
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
    P->extract("spot", spot2, size);
    P->extract("volatility", sigma2, size);
    P->extract("interest rate", r);
    if (P->extract("dividend rate", divid, size, true) == false)
    {
        divid = pnl_vect_create_from_zero(size);
    }
    P->extract("sample number", n_samples);
    P->extract("correlation", correlation);
    P->extract("timestep number", nbTimeStep);
    double t = 0;

    P->extract("payoff coefficients", payoff_coefficients, size);
    PnlVect *vect_stdev = pnl_vect_create_from_zero(size);
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    BlackScholesModel *blackScholesModel1 = new BlackScholesModel(size, r, correlation, sigma, spot);
    BlackScholesModel *blackScholesModel2 = new BlackScholesModel(size, r, correlation, sigma2, spot2);
    if (type == "basket")
    {
        P->extract("strike", strike);
        BasketOption *pBasketOption1 = new BasketOption(T, nbTimeStep, size, payoff_coefficients, strike);

        MonteCarlo *monteCarlo1 = new MonteCarlo(blackScholesModel1, pBasketOption1, rng, fdStep, n_samples);
        PnlVect *vect_stdev1 = pnl_vect_create(size);
        PnlVect *vect_stdev2 = pnl_vect_create(size);
        PnlVect *delta1 = pnl_vect_create_from_zero(size);
        PnlVect *delta2 = pnl_vect_create_from_zero(size);

        //monteCarlo1->price(price, stdev);
        // pnl_vect_print(delta1);
        PnlMat *past = pnl_mat_create_from_zero(size, nbTimeStep + 1);
        // blackScholesModel2->asset(past, T, nbTimeStep, rng);
        // pnl_mat_resize(past, size, floor(t * nbTimeStep) + 1);
        // pnl_mat_set_col(past, spot, 0);
         PnlMat *shift_path = pnl_mat_create(past->n, past->m);
        // monteCarlo1->delta(past, t, delta1, vect_stdev);
         //monteCarlo1->delta(delta1, vect_stdev1);
        //  pnl_vect_print(delta1);
        //  pnl_vect_print(vect_stdev1);
        //  blackScholesModel2->shiftAsset(shift_path, past, 1, 9, 0, T / nbTimeStep);
        //  pnl_mat_print(shift_path);
        //  pnl_mat_resize(past, size, floor(t * nbTimeStep) + 1);
        // monteCarlo1->price(past, t, price, stdev);

        // auto startP = std::chrono::high_resolution_clock::now();
        // monteCarlo1->paralleldelta(delta1, vect_stdev1);
        // auto finishP = std::chrono::high_resolution_clock::now();
        // std::chrono::duration<double> elapsedP = finishP - startP;
        // std::cout << "Time for parallel loop : " << elapsedP.count() << std::endl;
        // pnl_vect_print(delta1);
        // pnl_vect_free(&delta1);
        // pnl_vect_free(&vect_stdev1);

        // auto start = std::chrono::high_resolution_clock::now();
        // monteCarlo1->delta(delta2, vect_stdev2);
        // auto finish = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = finish - start;
        //std::cout << "Time for normal loop : " << elapsed.count() << std::endl;
        //pnl_vect_print(delta2);
        // // pnl_vect_print(vect_stdev);
        double PL1 = 0;
         int H = 3;
         PnlMat *matTot = pnl_mat_create_from_zero(size, H+1);
         blackScholesModel2->asset(matTot, T, nbTimeStep, rng);
         monteCarlo1->PL(matTot, PL1);
         printf(" portf = %f\n", PL1);

        //PricingResults res = PricingResults(price, stdev, delta1, vect_stdev1);
        //cout << res << endl;
        // pnl_vect_print(delta1);
        pnl_vect_free(&vect_stdev1); 
        pnl_vect_free(&vect_stdev2); 
        pnl_vect_free(&delta1); 
        pnl_vect_free(&delta2); 
        pnl_mat_free(&shift_path);
        pnl_mat_free(&past);
        delete (pBasketOption1);
        delete (monteCarlo1);
    }
    else if (type == "asian")
    {
        P->extract("strike", strike);
        PnlVect *deltaT = pnl_vect_create(size);
        //printf("Calcul de l'option Asiatique \n");
        AsianOption *pAsianOption1 = new AsianOption(T, nbTimeStep, size, strike, payoff_coefficients);
        MonteCarlo *monteCarlo1 = new MonteCarlo(blackScholesModel1, pAsianOption1, rng, T / nbTimeStep, n_samples);

        monteCarlo1->price(price, stdev);
        monteCarlo1->delta(deltaT, vect_stdev);

        PricingResults res = PricingResults(price, stdev, deltaT, vect_stdev);
        cout << res << endl;


        delete (pAsianOption1);
        delete (monteCarlo1);
    }
    else
    {

        PerformanceOption *pPerfOption1 = new PerformanceOption(T, nbTimeStep, size, payoff_coefficients);

        MonteCarlo *monteCarlo1 = new MonteCarlo(blackScholesModel1, pPerfOption1, rng, T / nbTimeStep, n_samples);

        PnlMat *past = pnl_mat_create_from_zero(size, nbTimeStep + 1);
        blackScholesModel2->asset(past, T, nbTimeStep, rng);
        pnl_mat_resize(past, size, floor(t * nbTimeStep) + 1);
        monteCarlo1->price(past, t, price, stdev);
        PnlVect *deltaT = pnl_vect_create(size);
        monteCarlo1->delta(deltaT, vect_stdev);
        
        PricingResults res = PricingResults(price, stdev, deltaT, vect_stdev);
        cout << res << endl;

        delete (pPerfOption1);
        delete (monteCarlo1);
    }
    pnl_vect_free(&divid);
    pnl_vect_free(&payoff_coefficients);
    delete (blackScholesModel1);
    delete (blackScholesModel2);
    delete (P);
    pnl_rng_free(&rng);
    pnl_vect_free(&vect_stdev);
    //std::cout << "le prix de l'option " << type << " est " << price << std::endl;
    // std::cout << "largeur de l'intervalle " << type << " est " << stdev << std::endl;

    //return 0;
}