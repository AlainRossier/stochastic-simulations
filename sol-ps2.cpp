#include <iostream>
#include <cstring>
#include <cmath>
#include <fstream>
#include <random>
#include <functional>

#include "utils.h"
#include "ps2.h"
#include "mlmc.h"


using namespace std;


void solutions_ps2() {

    // Initialize the random number generator
    default_random_engine rng;

    // Initialize the output file
    std::ofstream out;
    out.open(ABS_PATH + "sols/ps2/solutions_ps2.txt");

    // Print up to the 5th most significant digit
    out.precision(5);
    out.setf(ios::fixed);
    out.setf(ios::showpoint);



    // 2.a. Weak convergence of SDE
    out << "\n2.a. Weak convergence of SDE" << endl;
    double rate(0.05), sigma(0.2), maturity(1.0), initial_value(100.0), strike(100.0);
    size_t n_paths(1e6), n_levels(6);
    /*
    out << "Analysis of the weak convergence of GBM SDE with N = " << n_paths << + " paths and " << n_levels << " levels." << endl;
    european_call_path(rate, sigma, maturity, initial_value, strike, n_paths, n_levels,
                       ABS_PATH + "data/ps_2_2a_weak_convergence_sde.data", rng);

    // 2.b. Weak convergence of SDE with steps 2h
    out << "\n2.b. Weak convergence of SDE with steps 2h" << endl;
    out << "Analysis of the weak convergence of GBM SDE with N = " << n_paths << + " paths and " << n_levels << " levels, with the 2h trick." << endl;
    european_call_path_2h(rate, sigma, maturity, initial_value, strike, n_paths, n_levels,
                          ABS_PATH + "data/ps_2_2b_weak_convergence_sde_2h.data", rng);

    // 2.c. Strong convergence of SDE with steps 2h
    out << "\n2.c. Strong convergence of SDE with steps 2h" << endl;
    out << "Analysis of the strong convergence of GBM SDE with N = " << n_paths << + " paths and " << n_levels << " levels, with the 2h trick." << endl;
    gbm_strong_error(rate, sigma, maturity, initial_value, strike, n_paths, n_levels,
                     ABS_PATH + "data/ps_2_2c_strong_convergence_sde_2h.data", rng);

    */


    // 3. Ornstein-Uhlenbeck process
    out << "\n4. Strong convergence of the Ornstein Uhlenbeck process with steps 2h." << endl;
    out << "Analysis of the strong convergence of OU SDE with N = " << n_paths << + " paths and " << n_levels << " levels, with the 2h trick." << endl;
    double theta(110), kappa(2), sigma(0.5);
    ornstein_uhlenbeck_strong_error(rate, maturity, initial_value, theta, kappa, n_paths, sigma, n_levels, "data/ps_3_ohlenstein_uhlenbeck.data", rng);

    // 4. Heston stochastic volatility model
    out << "\n4. Strong convergence of the Heston stochastic volatility model with steps 2h." << endl;
    out << "Analysis of the strong convergence of GBM SDE with N = " << n_paths << + " paths and " << n_levels << " levels, with the 2h trick." << endl;
    double initial_vola(0.25), theta(0.25), kappa(2), zeta(0.5), rho(-0.1);
    heston_stochastic_vola_strong_error(rate, initial_vola, maturity, initial_value, theta, kappa, zeta, rho,
                                        n_paths, n_levels, ABS_PATH + "data/ps_4_strong_convergence_heston.data", rng);


    // 5. MLMC results

    /*
    % These are similar to the MLMC tests for the MCQMC06 paper % using a Milstein discretisation with 2^l timesteps on level l
    %
    % The figures are slightly different due to
    % -- change in MSE split
    % -- change in cost calculation
    % -- different random number generation
    % -- switch to S_0=100
    */

    // Main code

    int M  = 2;     // refinement cost factor
    int N0 = 200;   // initial samples on each level
    int Lmin = 2;   // minimum refinement level
    int Lmax = 5; //10;  // maximum refinement level

    int   N, L;
    float Eps[11];

    // loop over different payoff

    for (int option(1); option<2; option++) {
        rng.seed(1234);

        if (option==1) {
            cout << "\n ---- option " << option << ": European call ----" << endl;
            N      = 20000;    // samples for convergence tests
            L      = 8;        // levels for convergence tests
            float Eps2[] = { 0.005, 0.01, 0.02, 0.05, 0.1, 0.0 };
            memcpy(Eps,Eps2,sizeof(Eps2));
        }
        else if (option==2) {
            cout << "\n ---- option " << option << ": Asian call ----" << endl;
            N      = 20000;    // samples for convergence tests
            L      = 8;        // levels for convergence tests
            float Eps2[] = { 0.005, 0.01, 0.02, 0.05, 0.1, 0.0 };
            memcpy(Eps,Eps2,sizeof(Eps2));
        }
        else if (option==3) {
            cout << "\n ---- option " << option << ": lookback call ----" << endl;
            N      = 20000;    // samples for convergence tests
            L      = 10;       // levels for convergence tests
            float Eps2[] = { 0.005, 0.01, 0.02, 0.05, 0.1, 0.0 };
            memcpy(Eps,Eps2,sizeof(Eps2));
        }
        else if (option==4) {
            cout << "\n ---- option " << option << ": digital call ----" << endl;
            N      = 200000;   // samples for convergence tests
            L      = 8;        // levels for convergence tests
            float Eps2[] = { 0.01, 0.02, 0.05, 0.1, 0.2, 0.0 };
            memcpy(Eps,Eps2,sizeof(Eps2));
        }
        else if (option==5) {
            cout << "\n ---- option " << option << ": barrier call ----" << endl;
            N      = 200000;   // samples for convergence tests
            L      = 8;        // levels for convergence tests
            float Eps2[] = { 0.005, 0.01, 0.02, 0.05, 0.1, 0.0 };
            memcpy(Eps,Eps2,sizeof(Eps2));
        }

        mlmc_test(mcqmc06_l, M, N,L, N0,Eps, Lmin,Lmax, out, rng, option);

        // print exact analytic value, based on S0=K

        float T, r, sig, sig2, K, B, d1, d2, val, k, d3, d4;

        T   = 1.0f;
        r   = 0.05f;
        sig = 0.2f;
        sig2 = sig*sig;
        K   = 100.0f;
        B   = 0.85f*K;
        d1  = (r+0.5f*sig2)*T / (sig*sqrtf(T));
        d2  = (r-0.5f*sig2)*T / (sig*sqrtf(T));

        if (option==1)
            val = K*( ncff(d1) - expf(-r*T)*ncff(d2) );
        else if (option==2)
            val = nanf("");
        else if (option==3) {
            k   = 0.5f*sig2/r;
            val = K*( ncff(d1) - ncff(-d1)*k - exp(-r*T)*(ncff(d2) - ncff(d2)*k) );
        }
        else if (option==4)
            val = K*expf(-r*T)*ncff(d2);
        else {
            k   = 0.5f*sig2/r;
            d3  = (2.0f*logf(B/K) + (r+0.5f*sig2)*T) / (sig*sqrtf(T));
            d4  = (2.0f*logf(B/K) + (r-0.5f*sig2)*T) / (sig*sqrtf(T));
            val = K*(ncff(d1) - expf(-r*T)*ncff(d2) -
                     powf(K/B,1.0f-1.0f/k)*((B*B)/(K*K)*ncff(d3) - expf(-r*T)*ncff(d4)));
        }

    // Do 100 MLMC calcs in parallel

    // mlmc_test_100(mcqmc06_l, val, N0,Eps,Lmin,Lmax, out);

  }





    out.close();

}

