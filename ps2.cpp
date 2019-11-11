#include <iostream>
#include <cmath>
#include "ps1.h"
#include <fstream>


using namespace std;


// Problem 2

void european_call_path(double rate, double sigma, double maturity, double initial_value, double strike,
                        size_t n_paths, size_t n_levels,
                        string path, default_random_engine& rng) {

    // Initialize
    size_t incr(1);
    double step(0.0);
    double sum, sumsq, price, dW, final_price, mean, sd;

    // Analytical price
    double value_european_call = analytical_european_call(rate, sigma, maturity, initial_value, strike, "value");

    // Random number generator
    normal_distribution<double> normal(0.0f, 1.0f);
    auto next_normal = bind(ref(normal), ref(rng));

    // Saving
    std::ofstream outfile;
    outfile.open(path);

    for (size_t p(1); p <= n_levels; p++) {
        incr *= 2;
        step = maturity / incr;
        sum = 0.0; sumsq = 0.0;
        for (size_t m(0); m < n_paths; m++) {
            price = initial_value;
            for (size_t i(0); i < incr; i++) {
                dW = sqrt(step) * next_normal();
                price *= (1 + rate*step + sigma*dW);
            }
            final_price = exp(-rate*maturity) * max(price-strike, 0.0);
            sum += final_price;
            sumsq += pow(final_price, 2);
        }

        mean = sum/n_paths;
        sd = sqrt((sumsq/n_paths - pow(mean, 2)) / (n_paths-1));
        outfile << step << " " << mean - value_european_call << " " << 3*sd << "\n";
    }

    outfile.close();

}


void european_call_path_2h(double rate, double sigma, double maturity, double initial_value, double strike,
                           size_t n_paths, size_t n_levels,
                           string path, default_random_engine& rng) {

    // Initialize
    size_t incr(1);
    double step(0.0);
    double sum, sumsq, price, price_2h, dW, dW_2h, final_price, final_price_2h, mean, sd;

    // Analytical price
    double value_european_call = analytical_european_call(rate, sigma, maturity, initial_value, strike, "value");

    // Random number generator
    normal_distribution<double> normal(0.0f, 1.0f);
    auto next_normal = bind(ref(normal), ref(rng));

    // Saving
    std::ofstream outfile;
    outfile.open(path);

    for (size_t p(1); p <= n_levels; p++) {
        incr *= 2;
        step = maturity / incr;
        sum = 0.0; sumsq = 0.0;
        for (size_t m(0); m < n_paths; m++) {
            price = initial_value;
            price_2h = initial_value;
            for (size_t i(0); i < incr; i++) {
                dW = sqrt(step) * next_normal();
                price *= (1 + rate*step + sigma*dW);
                dW_2h = sqrt(step) * next_normal();
                price_2h *= (1 + rate*2*step + sigma*(dW + dW_2h));
            }
            final_price = exp(-rate*maturity) * max(price-strike, 0.0);
            final_price_2h = exp(-rate*maturity) * max(price_2h-strike, 0.0);

            sum += (final_price - final_price_2h);
            sumsq += pow(final_price - final_price_2h, 2);
        }

        mean = sum/n_paths;
        sd = sqrt((sumsq/n_paths - pow(mean, 2)) / (n_paths-1));
        outfile << step << " " << mean - value_european_call << " " << 3*sd << "\n";
    }

    outfile.close();

}
