#include <iostream>
#include <random>
#include <functional>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/distributions/normal.hpp>


using namespace std;
using namespace boost::numeric::ublas;


// Problem 1

void populate(std::vector<double>& samples, string method, default_random_engine& rng)
{
    if (method == "uniform")
    {
        uniform_real_distribution<double> uniform(0.0f, 1.0f);
        auto next_uniform = bind(ref(uniform), ref(rng));
        generate(samples.begin(), samples.end(), next_uniform);
    }
    else if (method == "normal")
    {
        normal_distribution<double> normal(0.0f, 1.0f);
        auto next_normal = bind(ref(normal), ref(rng));
        generate(samples.begin(), samples.end(), next_normal);
    }
    else
    {
        cerr << "This method is not currently supported." << endl;
    }
}


matrix<double> cholesky(matrix<double> input) {
    if (input.size1() != input.size2() || input.size1() == 0) {
        cerr << "The input matrix is not square." << endl;
    }
    else {
        size_t dim = input.size1();
        matrix<double> output(dim, dim, 0);
        output(0, 0) = sqrt(input(0, 0));
        for (size_t j(1); j < dim; j++) {
            output(j, 0) = input(j, 0) / output(0, 0);
        }
        for (size_t i(1); i < dim; i++) {
            double acc(0.0);
            for (size_t k(0); k < i; k++) {acc += pow(output(i, k), 2);}
            output(i, i) = sqrt(input(i, i) - acc);
            for (size_t j(i+1); j < dim; j++) {
                acc = 0.0;
                for (size_t k(0); k < i; k++) {acc += output(i, k) * output(j, k);}
                output(j, i) = (input(j, i) - acc) / output(i, i);
            }
        }
        return output;
    }
}


// Problem 2
double linear_cos(double unif) {
    return unif * cos(M_PI * unif);
}

double empiricalMean(size_t n_sim, const std::function<double(double)>& f, default_random_engine& rng) {
    double acc(0.0);
    size_t incr(0);
    uniform_real_distribution<double> uniform(0.0f, 1.0f);
    auto next_uniform = bind(ref(uniform), ref(rng));
    while (incr < n_sim) {
        acc += f(next_uniform());
        incr += 1;
    }

    return acc / n_sim;
}


double norm_cdf(double normal) {
    return 0.5 * erfc(-normal * M_SQRT1_2);
}

void confianceIntervals(std::vector<size_t> n_sims, const std::function<double(double)>& f, double true_mean,
                        string path, std::default_random_engine& rng) {
    // Initialization
    size_t n_experiments = n_sims.size();
    std::vector<double> empirical_mean(n_experiments);
    std::vector<double> empirical_sd(n_experiments);
    std::sort(n_sims.begin(), n_sims.end());
    size_t max_sim = n_sims.back();
    double sum(0.0), sumsq(0.0), acc(0.0);
    size_t incr(0), k(0);

    // Generation
    uniform_real_distribution<double> uniform(0.0f, 1.0f);
    auto next_uniform = bind(ref(uniform), ref(rng));
    while (incr < max_sim) {
        acc = f(next_uniform());
        sum += acc; sumsq += pow(acc, 2);
        incr += 1;
        if (incr == n_sims[k]) {
            empirical_mean[k] = sum / incr;
            empirical_sd[k] = sqrt((sumsq - incr * pow(empirical_mean[k], 2)) / (incr-1));
            k += 1;
        }
    }

    // Saving
    std::ofstream outfile_ci;
    outfile_ci.open(path);
    for (size_t i(0); i < n_experiments; i++) {
        double n = n_sims[i];
        double m = empirical_mean[i];
        double sd = empirical_sd[i];
        outfile_ci << n << " " << m << " " << m-3*sd/sqrt(n) << " " << m+3*sd/sqrt(n) << " " << true_mean << "\n";
    }
    outfile_ci.close();
}


double analytical_european_call(double rate, double sigma, double maturity, double initial_value, double strike, string type) {

    initial_value = max(initial_value, pow(10, -100));
    strike = max(strike, pow(10, -100));
    double dplus = (log(initial_value) - log(strike) + (rate + 0.5 * pow(sigma, 2)) * maturity) / (sigma * sqrt(maturity));
    double dminus = (log(initial_value) - log(strike) + (rate - 0.5 * pow(sigma, 2)) * maturity) / (sigma * sqrt(maturity));

    if (type == "value") {
        return initial_value * norm_cdf(dplus) - exp(-rate * maturity) * strike * norm_cdf(dminus);
    } else if (type == "delta") {
        return norm_cdf(dplus);
    } else if (type == "gamma") {
        return exp(-0.5*pow(dplus, 2)) / (sigma * sqrt(2*M_PI*maturity) * initial_value);
    } else if (type == "vega") {
        return initial_value * exp(-0.5*pow(dplus, 2)) / sqrt(2*M_PI) * (sqrt(maturity) - dplus/sigma) - \
            exp(-rate * maturity) * strike * exp(-0.5*pow(dminus, 2)) / sqrt(2*M_PI) * (-sqrt(maturity) - dminus/sigma);
    } else {
        cerr << "This greek is not available." << endl;
        return 0.0;
    }
}


double payoff_european_call(double unif, double rate, double sigma, double maturity,
                            double initial_value, double strike, boost::math::normal& dist) {
    double normal = quantile(dist, unif);
    double price = initial_value * exp((rate - 0.5*pow(sigma, 2)) * maturity + sigma * sqrt(maturity) * normal);
    return exp(-rate * maturity) * max(price - strike, 0.0);
}


// Problem 4

double antitheticVariables(std::vector<size_t> n_sims, const std::function<double(double)>& f,
                        string path, std::default_random_engine& rng) {
    // Initialization
    size_t n_experiments = n_sims.size();
    std::vector<double> empirical_variance(n_experiments);
    std::vector<double> empirical_variance_anti(n_experiments);
    std::sort(n_sims.begin(), n_sims.end());
    size_t max_sim = n_sims.back();
    double sample_plus(0.0), sample_minus(0.0);
    double sum(0.0), sumsq(0.0), acc(0.0);
    double sum_anti(0.0), sumsq_anti(0.0), acc_anti(0.0);
    double corr(0.0);
    size_t incr(0), k(0);

    // Generation
    uniform_real_distribution<double> uniform(0.0f, 1.0f);
    auto next_uniform = bind(ref(uniform), ref(rng));
    while (incr < max_sim) {
        double u = next_uniform();
        sample_plus = f(u), sample_minus = f(1-u);
        corr += sample_plus * sample_minus;
        sum += sample_plus; sumsq += pow(sample_plus, 2);
        sum_anti += 0.5 * (sample_plus + sample_minus);
        sumsq_anti += pow(0.5 * (sample_plus + sample_minus), 2);
        incr += 1;
        if (incr == n_sims[k]) {
            empirical_variance[k] = (sumsq - incr * pow(sum/incr, 2)) / (incr*(incr-1));
            empirical_variance_anti[k] = (sumsq_anti - incr * pow(sum_anti/incr, 2)) / (incr*(incr-1));
            k += 1;
        }
    }

    corr = (corr/incr - pow(sum/incr, 2)) / (incr * empirical_variance[k-1]);

    // Saving
    std::ofstream outfile_ci;
    outfile_ci.open(path);
    for (size_t i(0); i < n_experiments; i++) {
        outfile_ci << n_sims[i] << " " << empirical_variance[i]  << " " << empirical_variance_anti[i] << "\n";
    }
    outfile_ci.close();

    return corr;
}

