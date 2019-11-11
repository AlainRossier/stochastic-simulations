#include <iostream>
#include <cmath>
#include <fstream>

#include "utils.h"
#include "ps2.h"


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
    out << "Analysis of the weak convergence of GBM SDE with N = " << n_paths << + " paths and " << n_levels << " levels." << endl;
    european_call_path(rate, sigma, maturity, initial_value, strike, n_paths, n_levels,
                       ABS_PATH + "data/ps_2_2a_weak_convergence_sde.data", rng);

    // 2.b. Weak convergence of SDE with steps 2h
    out << "\n2.b. Weak convergence of SDE with steps 2h" << endl;
    out << "Analysis of the weak convergence of GBM SDE with N = " << n_paths << + " paths and " << n_levels << " levels, with the 2h trick." << endl;
    european_call_path_2h(rate, sigma, maturity, initial_value, strike, n_paths, n_levels,
                          ABS_PATH + "data/ps_2_2b_weak_convergence_sde_2h.data", rng);






    out.close();

}
