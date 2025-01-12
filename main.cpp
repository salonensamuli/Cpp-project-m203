// Paris Dauphine-PSL M203 introduction to c++ 
// course home assignment by Samuli Salonen (22301060)
// main source: https://rtraba.com/wp-content/uploads/2015/05/solving_pdes_in_c.pdf

#include <iostream>
#include <iomanip>
#include <cmath> // for abs
#include <vector> // for matrices
#include <stdexcept> // for errors
#include <chrono> // for timer
#include "pdepricer.h"
#include "matrix.h"
using namespace std;

double cdfNormal(double x) {
    // estimate using the erfc (complementary error function)
    return 0.5 * erfc(-x / sqrt(2.0));
}

double blackScholesPricer(double S, double K, double T, double r, double sigma, int type = 1) {
    // type = 1 for call and type = -1 for put
    double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    if(type == 1) {
        return S*cdfNormal(d1) - K*exp(-r*T)*cdfNormal(d2);
    } else {
        return K*exp(-r*T)*cdfNormal(-d2) - S*cdfNormal(-d1);
    }
}

int main() {

    // market parameters
    double T = 1.0; // time to maturity in years
    double r = 0.05; // risk-free rate
    double sigma = 0.2; // annual volatility
    double K = 100.0; // option strike price
    double S0 = 100.0; // underlying spot price
    int type = 1; // option type, 1 for call and -1 for put

    // Black-Scholes closed form pricer
    double bsPrice = blackScholesPricer(S0, K, T, r, sigma, type);

    // PDE pricer
    vector<int> nmValues = {25, 50, 75, 100, 125, 150, 175, 200}; // equivalent time and space steps
    double Smax = 2.0 * K; // upper boundary in S

    cout << fixed << setprecision(3);
    cout << "Black-Scholes price: " << bsPrice << "\n\n";
    cout << "n=m   PDE-price      Error   Time (s)\n";
    cout << "--------------------------------------\n";

    // loop to study the convergence of the PDE pricer
    for (int nm : nmValues) {
        int n = nm;
        int m = nm;
        auto start = chrono::high_resolution_clock::now();
        PDEPricer pricer(T, n, m, r, sigma, K, Smax);
        double pdePrice = pricer.solve();
        auto end = chrono::high_resolution_clock::now();
        double elapsedTime = chrono::duration<double>(end - start).count();
        double abs_error = abs(bsPrice - pdePrice);
        cout << setw(3) << nm << 
                setw(12) << pdePrice << " " << 
                setw(10) << abs_error << " " << 
                setw(10) << elapsedTime << "\n";
    }
    return 0;
}