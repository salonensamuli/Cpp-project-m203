// Paris Dauphine-PSL M203 introduction to c++ 
// course home assignment by Samuli Salonen (22301060)
// main source: https://rtraba.com/wp-content/uploads/2015/05/solving_pdes_in_c.pdf

#include "pdepricer.h"
#include <algorithm> // for max
#include <cmath> // for exp and aabs
#include <vector> // for resize
#include <stdexcept> // for errors
using namespace std;

// constructor
PDEPricer::PDEPricer(double T, int n, int m, double r, double sigma, double K, double Smax)
    : T_(T), n_(n), m_(m), r_(r), sigma_(sigma), Smax_(Smax), K_(K)
{
    // time grid
    double dt = T_ / static_cast<double>(n_); // equal step sizes
    timeGrid_.resize(n_+1);
    for(int i=0; i<=n_; ++i) {
        timeGrid_[i] = i*dt; // from dt to T
    }

    // space grid
    double dS = Smax_ / static_cast<double>(m_-1); // equal step sizes
    spaceGrid_.resize(m_);
    for(int j=0; j<m_; ++j) {
        spaceGrid_[j] = j * dS; // from 0 to Smax
    }

    // boundary conditions for the price
    boundaryLower_.resize(n_+1, 0.0); // lower boundary is 0
    boundaryUpper_.resize(n_+1, 0.0);
    for(int i=0; i<=n_; ++i) {
        double t = timeGrid_[i];
        double tau = (T_ - t);
        boundaryLower_[i] = 0.0; 
        boundaryUpper_[i] = Smax_ - K_ * exp(-r_ * tau); // upper boundary is Smax - K*e^(-r*tau) where tau is T-t 
    }

    // terminal condition at maturity
    terminalCond_.resize(m_);
    for(int j=0; j<m_; ++j) {
        double S = spaceGrid_[j];
        terminalCond_[j] = max(S - K_, 0.0); // V(S,T) = max(S-K, 0)
    }
}

// solve the PDE with finite-differences (chapters 7 and 8)
double PDEPricer::solve() {
    // V_prev stores the option value at the current time step (then swapped with V_curr)
    // V_curr will hold the option value at the next time step
    vector<double> V_prev(m_), V_curr(m_);

    // initialize the option value at maturity (terminal condition)
    for (int j = 0; j < m_; ++j) {
        V_prev[j] = terminalCond_[j]; // Payoff at t = T
    }

    // time and space step sizes
    double dt = T_ / static_cast<double>(n_);
    double dS = spaceGrid_[1] - spaceGrid_[0];

    // matrix A for the implicit scheme
    Matrix A(m_, m_);

    // time-stepping backward from t = T to t = 0
    for (int step = n_ - 1; step >= 0; --step) {
        double t = timeGrid_[step]; // current time t
        double tau = T_ - t; // time to maturity
        // initialize the matrix A to zero
        for (int i = 0; i < m_; ++i) {
            for (int j = 0; j < m_; ++j) {
                A(i, j) = 0.0;
            }
        }
        // build the tri-diagonal matrix A and the right-hand side (rhs) vector
        // PDE: dV/dt = 0.5*sigma^2*S^2*V' + r*S*V' - r*V
        // rearrange for implicit scheme: A * V^(n+1) = V^n
        vector<double> rhs = V_prev; // rhs starts with V_prev

        for (int j = 1; j < m_ - 1; ++j) {
            double S = spaceGrid_[j]; // current spot

            // coefficients of the PDE terms
            double alpha = 0.5 * sigma_ * sigma_ * S * S / (dS * dS); // diffusion term
            double beta  = 0.5 * r_ * S / dS; // convection term
            double gamma = r_; // discount term

            // tri-diagonal coefficients
            double a = -dt * (alpha - beta); // sub-diagonal (j-1)
            double b = 1.0 + dt * (2.0 * alpha + gamma); // diagonal (j)
            double c = -dt * (alpha + beta); // super-diagonal (j+1)

            // populate the matrix A
            A(j, j - 1) = a; // sub-diagonal
            A(j, j) = b; // main diagonal
            A(j, j + 1) = c; // super-diagonal
        }

        // lower boundary (j = 0): V(0, t) = 0 for a call option
        A(0, 0) = 1.0; // set row to enforce V = 0
        for (int c = 1; c < m_; ++c) {
            A(0, c) = 0.0;
        }
        rhs[0] = boundaryLower_[step]; // boundary value

        // upper boundary (j = m-1): V(Smax, t) = Smax - K*exp(-r*tau)
        A(m_ - 1, m_ - 1) = 1.0; // set row to enforce boundary
        for (int c = 0; c < m_ - 1; ++c) {
            A(m_ - 1, c) = 0.0;
        }
        rhs[m_ - 1] = boundaryUpper_[step]; // boundary value

        // solve the linear system A * V_curr = rhs
        // check if invertible
        if (!A.isInvertible()) {
            throw runtime_error("Matrix not invertible! (singular?)");
        }

        Matrix Ainv = A.inverse(); // inverse of A

        // multiply Ainv by rhs to get V_curr
        for (int i = 0; i < m_; ++i) {
            double sum = 0.0;
            for (int j = 0; j < m_; ++j) {
                sum += Ainv(i, j) * rhs[j];
            }
            V_curr[i] = sum; // update the solution for the current time step
        }
        V_prev = V_curr; // to prepare for the next iteration V_curr becomes V_prev
    }
    // return the option value at S=K or closest to K
    double minDiff = 1e10; // min difference between S and K
    int idxStar = 0; // index of the closest grid point to K
    for (int j = 0; j < m_; ++j) {
        double diff = abs(spaceGrid_[j] - K_);
        if (diff < minDiff) {
            minDiff = diff;
            idxStar = j; // update the closest index
        }
    }
    return V_curr[idxStar]; // option price at S=K
}