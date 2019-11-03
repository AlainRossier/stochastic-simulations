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
        return matrix<double>(0, 0);
    }
    else {
        size_t dim = input.size1();
        matrix<double> output(dim, dim, 0.0);
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


matrix<double> pca2d(matrix<double> input) {
    if (input.size1() != 2 || input.size2() != 2) {
        cerr << "The input is not a 2x2 square matrix." << endl;
        return matrix<double>(0, 0);
    }
    else if (input(0, 1) != input(1, 0) || input(0,0)*input(1,1) - input(0, 1)*input(1, 0) < 0) {
        cerr << "The input is not symmetric positive definite." << endl;
        return matrix<double>(0, 0);
    }
    else if (input(1, 0) == 0) {
        matrix<double> output(2, 2, 0.0);
        output(0, 0) = sqrt(input(0, 0));
        output(1, 1) = sqrt(input(1, 1));
        return output;
    }
    else {
        double c = input(1, 0);
        double trace = input(0, 0) + input(1, 1);
        double delta = sqrt(pow(input(0, 0) - input(1, 1), 2) + 4*pow(c, 2));
        double l0 = (trace + delta) / 2;
        double l1 = (trace - delta) / 2;
        double v00 = abs(c) / sqrt(pow(l0-input(0, 0), 2) + pow(c, 2));
        double v01 = (l0 - input(0, 0)) * v00 / c;
        double v10 = abs(c) / sqrt(pow(l1-input(0, 0), 2) + pow(c, 2));
        double v11 = (l1 - input(0, 0)) * v10 / c;

        matrix<double> output(2, 2, 0.0);

        output(0, 0) = sqrt(l0) * v00;
        output(1, 0) = sqrt(l0) * v01;
        output(0, 1) = sqrt(l1) * v10;
        output(1, 1) = sqrt(l1) * v11;

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
                        std::vector<double>& empirical_mean, std::vector<double>& empirical_sd,
                        string path, std::default_random_engine& rng) {
    // Initialization
    size_t n_experiments = n_sims.size();
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
            empirical_sd[k] = sqrt((sumsq - incr * pow(empirical_mean[k], 2)) / (incr*(incr-1)));
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
        outfile_ci << n << " " << m << " " << m-3*sd << " " << m+3*sd << " " << true_mean << "\n";
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


double payoff_european_call(double normal, double rate, double sigma, double maturity,
                            double initial_value, double strike) {
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
    double sum(0.0), sumsq(0.0);
    double sum_anti(0.0), sumsq_anti(0.0);
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


double controlVariate(std::vector<size_t> n_sims, const std::function<double(double)>& f,
                      const std::function<double(double)>& g, double true_mean_g,
                      string path, std::default_random_engine& rng) {
    // Initialization
    size_t n_experiments = n_sims.size();
    std::vector<double> empirical_variance(n_experiments);
    std::vector<double> empirical_variance_control(n_experiments);
    std::sort(n_sims.begin(), n_sims.end());
    size_t max_sim = n_sims.back();
    double sample_f(0.0), sample_g(0.0);
    double sum_f(0.0), sumsq_f(0.0), var_f(0.0);
    double sum_g(0.0), sumsq_g(0.0), var_g(0.0);
    double sum_control(0.0), sumsq_control(0.0), acc_control(0.0);
    double covar(0.0), corr(0.0), lambda(1.0);
    size_t incr(0), k(0);

    // Generation
    uniform_real_distribution<double> uniform(0.0f, 1.0f);
    auto next_uniform = bind(ref(uniform), ref(rng));
    double u(0.0);
    while (incr < max_sim) {
        u = next_uniform();
        sample_f = f(u), sample_g = g(u);
        sum_f += sample_f; sumsq_f += pow(sample_f, 2);
        sum_g += sample_g; sumsq_g += pow(sample_g, 2);
        acc_control = sample_f - lambda*(sample_g - true_mean_g);
        sum_control += acc_control; sumsq_control += pow(acc_control, 2);
        incr += 1;

        // Update the lambda for the first 1000 iterations
        if (incr <= 1000) {
            covar += sample_f * sample_g;
            if (incr >= 10) {
                var_f = (sumsq_f - incr * pow(sum_f/incr, 2)) / (incr-1);
                var_g = (sumsq_g - incr * pow(sum_g/incr, 2)) / (incr-1);
                lambda = (covar/incr - sum_f/incr * sum_g/incr) / var_g;
                corr = lambda * sqrt(var_g/var_f);
            }
        }
        if (incr == n_sims[k]) {
            empirical_variance[k] = (sumsq_f - incr * pow(sum_f/incr, 2)) / (incr*(incr-1));
            empirical_variance_control[k] = (sumsq_control - incr * pow(sum_control/incr, 2)) / (incr*(incr-1));
            k += 1;
        }
    }

    // Saving
    std::ofstream outfile_ci;
    outfile_ci.open(path);
    for (size_t i(0); i < n_experiments; i++) {
        outfile_ci << n_sims[i] << " " << empirical_variance[i]  << " " << empirical_variance_control[i] << "\n";
    }
    outfile_ci.close();

    return corr;
}


// Problem 5

double payoff_digital_put(double normal, double rate, double sigma, double maturity,
                          double initial_value, double strike) {
    double price = initial_value * exp((rate - 0.5*pow(sigma, 2)) * maturity + sigma * sqrt(maturity) * normal);
    return exp(-rate * maturity) * (strike - price > 0.0);
}


// Problem 6
void finite_difference_bumping(std::vector<double> bumps, size_t n_sim, const std::function<double(double, double)>& f,
                               double true_sensi, std::string path, std::default_random_engine& rng) {

     // Initialization
    size_t n_experiments = bumps.size();
    std::vector<double> sum_sensi_same_random(n_experiments);
    std::vector<double> sum_sensi_diff_random(n_experiments);
    std::vector<double> sumsq_sensi_same_random(n_experiments);
    std::vector<double> sumsq_sensi_diff_random(n_experiments);
    std::sort(bumps.begin(), bumps.end());
    double acc_same_random(0.0), acc_diff_random(0.0);

    // Generation
    uniform_real_distribution<double> uniform(0.0f, 1.0f);
    auto next_uniform = bind(ref(uniform), ref(rng));
    double u(0.0), v(0.0);
    for (size_t i(0); i < n_sim; i++) {
        u = next_uniform(); v = next_uniform();
        for (size_t k(0); k < n_experiments; k++) {
            acc_same_random = (f(u, bumps[k]) - f(u, -bumps[k])) / (2*bumps[k]);
            acc_diff_random = (f(u, bumps[k]) - f(v, -bumps[k])) / (2*bumps[k]);
            sum_sensi_same_random[k] += acc_same_random;
            sum_sensi_diff_random[k] += acc_diff_random;
            sumsq_sensi_same_random[k] += pow(acc_same_random, 2);
            sumsq_sensi_diff_random[k] += pow(acc_diff_random, 2);
        }
    }

    // Saving
    std::ofstream outfile_fd;
    outfile_fd.open(path);
    double mean_same(0.0), std_same(0.0), mean_diff(0.0), std_diff(0.0);
    for (size_t k(0); k < n_experiments; k++) {
        mean_same = sum_sensi_same_random[k]/n_sim;
        std_same = sqrt((sumsq_sensi_same_random[k] - n_sim * pow(mean_same, 2)) / (n_sim*(n_sim-1)));
        mean_diff = sum_sensi_diff_random[k]/n_sim;
        std_diff = sqrt((sumsq_sensi_diff_random[k] - n_sim * pow(mean_diff, 2)) / (n_sim*(n_sim-1)));

        outfile_fd << bumps[k] << " " << mean_same << " " << mean_same-3*std_same << " " << mean_same+3*std_same << " ";
        outfile_fd << mean_diff << " " << mean_diff-3*std_diff << " " << mean_diff+3*std_diff << " ";
        outfile_fd << true_sensi << "\n";
    }
    outfile_fd.close();
}
