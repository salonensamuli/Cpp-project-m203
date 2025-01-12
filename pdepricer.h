// Paris Dauphine-PSL M203 introduction to c++ 
// course home assignment by Samuli Salonen (22301060)
// main source: https://rtraba.com/wp-content/uploads/2015/05/solving_pdes_in_c.pdf

#ifndef PDEPRICER_H
#define PDEPRICER_H

#include <vector>
#include "matrix.h"
using namespace std;

class PDEPricer {
private:
    double T_; // time to maturity in years
    int n_; // number of time steps
    int m_; // number of space steps
    double r_; // risk-free rate
    double sigma_; // constant volatility
    double Smax_; // max underlying price (i put this 2x strike)
    double K_; // strike price

    vector<double> timeGrid_; // (chapter 7)
    vector<double> spaceGrid_; // (chapter 7)
    vector<double> boundaryLower_; 
    vector<double> boundaryUpper_; 
    vector<double> terminalCond_;  

public:
    // constructor with similar parameters as above
    PDEPricer(double T, int n, int m, double r, double sigma, double Smax, double K);
    
    // solve function
    double solve();
};

#endif