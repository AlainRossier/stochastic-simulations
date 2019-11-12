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
        outfile << step << " " << abs(mean - value_european_call) << " " << 3*sd << "\n";
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
                dW_2h = sqrt(step) * next_normal();
                price *= ((1 + rate*step + sigma*dW) * (1 + rate*step + sigma*dW_2h));
                price_2h *= (1 + rate*2*step + sigma*(dW + dW_2h));
            }
            final_price = exp(-rate*maturity) * max(price-strike, 0.0);
            final_price_2h = exp(-rate*maturity) * max(price_2h-strike, 0.0);

            sum += (final_price - final_price_2h);
            sumsq += pow(final_price - final_price_2h, 2);
        }

        mean = sum/n_paths;
        sd = sqrt((sumsq/n_paths - pow(mean, 2)) / (n_paths-1));
        outfile << step << " " << abs(mean) << " " << 3*sd << "\n";
    }

    outfile.close();

}


void gbm_strong_error(double rate, double sigma, double maturity, double initial_value, double strike,
                      size_t n_paths, size_t n_levels,
                      string path, default_random_engine& rng) {

    // Initialize
    size_t incr(1);
    double step(0.0);
    double sum, sumsq, sum_2h, sumsq_2h;
    double price, price_2h, dW1, dW2, brownian;
    double exact_price, err, err_2h;
    double mean, sd, mean_2h, sd_2h;

    // Random number generator
    normal_distribution<double> normal(0.0f, 1.0f);
    auto next_normal = bind(ref(normal), ref(rng));

    // Saving
    std::ofstream outfile;
    outfile.open(path);

    for (size_t p(1); p <= n_levels; p++) {
        incr *= 2;
        step = maturity / incr;
        sum = 0.0; sumsq = 0.0; sum_2h = 0.0; sumsq_2h = 0.0;

        for (size_t m(0); m < n_paths; m++) {
            price = initial_value;
            price_2h = initial_value;
            brownian = 0.0;

            for (size_t i(0); i < incr/2; i++) {
                dW1 = sqrt(step) * next_normal();
                dW2 = sqrt(step) * next_normal();
                price *= ((1 + rate*step + sigma*dW1) * (1 + rate*step + sigma*dW2));
                price_2h *= (1 + rate*2*step + sigma*(dW1+dW2));
                brownian += (dW1+dW2);
            }
            exact_price = initial_value * exp((rate - 0.5*pow(sigma, 2))*maturity + sigma*brownian);

            err = pow(price - exact_price, 2);
            sum += err; sumsq += pow(err, 2);

            err_2h = pow(price_2h - exact_price, 2);
            sum_2h += err_2h; sumsq_2h += pow(err_2h, 2);
        }

        mean = sum/n_paths;
        sd = sqrt((sumsq/n_paths - pow(mean, 2)) / (n_paths-1));
        mean_2h = sum_2h/n_paths;
        sd_2h = sqrt((sumsq_2h/n_paths - pow(mean_2h, 2)) / (n_paths-1));

        outfile << step << " " << sqrt(mean) << " " << (0.5/sqrt(mean)) * 3*sd << " ";
        outfile << sqrt(mean_2h) << " " << (0.5/sqrt(mean_2h)) * 3*sd_2h << endl;
    }

}
