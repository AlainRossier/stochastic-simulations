#include <iostream>
#include <algorithm>
#include <random>
#include <fstream>

#include "utils.h"
#include "ps1.h"
#include <boost/numeric/ublas/io.hpp>


using namespace std;
using namespace boost::numeric::ublas;

std::string ABS_PATH = "/home/alain/Documents/phd/lectures/stochastic_simulations/stochastic-simulations/";


int main()
{

    default_random_engine rng;

    // 1.a. Random variable generation

    clock_t start;
    start = clock();

    size_t n_sim(1000000);

    std::vector<double> sample_uniform(n_sim);
    populate(sample_uniform, "uniform", rng);
    cout << "CPU time elapsed : " << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms. \n" << endl;
    cout << "Empirical mean of the uniform distribution : " << mean(sample_uniform) << endl;
    cout << "Empirical variance of the uniform distribution : " << variance(sample_uniform) << "\n" << endl;


    // You can see that for n=10^7, the empirical variance is quite far away from the
    // theoretical one (order 10^(-2)). But it is fast (1.1s)
    // The reason is that we are using float, which are represented by 32 bits.
    // Using double, represented by 64 bits, the absolute error on the variance is of order 10^(-3),
    // at the expense of speed (1.4s).

    std::vector<double> sample_normal(n_sim);
    populate(sample_normal, "normal", rng);

    cout << "Empirical mean of the normal distribution : " << mean(sample_normal) << endl ;
    cout << "Empirical variance of the normal distribution : " << variance(sample_normal) << "\n" << endl ;


    // 1.b. Cholesky decomposition
    matrix<double> A(2, 2);
    A(0, 0) = 4; A(0, 1) = 1; A(1, 0) = 1; A(1, 1) = 4;
    matrix<double> L = cholesky(A);

    std::vector<double> sample_normal_coor(2*n_sim);
    populate(sample_normal_coor, "normal", rng);
    // Copy the elements in a matrix
    unbounded_array<double> storage(2*n_sim);
    std::copy(sample_normal_coor.begin(), sample_normal_coor.end(), storage.begin());
    matrix<double> multivariate_normal(n_sim, 2, storage);
    matrix<double> corr_multivariate_normal = prod(multivariate_normal, trans(L));

    cout << "Empirical mean of the multivariate normal distribution : " << mean(corr_multivariate_normal) << endl ;
    cout << "Empirical covariance of the multivariate normal distribution : " << covariance(corr_multivariate_normal) << "\n" << endl ;


    // 1.c. PCA factorization


    // 1.d


    // 2.a. Analytical computation
    double f_bar = -2/pow(M_PI, 2);
    double std_dev = sqrt(1.0/6 - 4/pow(M_PI, 4) + 1.0/(4*pow(M_PI, 2)));
    std::function<double(double)> linear_cos_lambda = [=](double x) -> double {return x;};

    // 2.b. Monte-Carlo simulations
    size_t inner_sim(1000);
    size_t outer_sim(1000);
    std::vector<double> cdf(outer_sim);
    for (size_t i(0); i < outer_sim; i++) {
        cdf[i] = empiricalMean(inner_sim, linear_cos_lambda, rng);
    }

    // 2.c. Sorting and plotting
    std::sort(cdf.begin(), cdf.end());
    std::ofstream outfile_2c;
    outfile_2c.open(ABS_PATH + "data/ps_1_2c_empirical_cdf.data");
    size_t count(0);
    for (const auto &e : cdf) {
        outfile_2c << e << " " << (count + 0.5) / outer_sim << " " << norm_cdf(sqrt(outer_sim) * (e-f_bar) / std_dev) << "\n";
        count++;
    }
    outfile_2c.close();


    // 2.d. Confidence intervals
    std::vector<size_t> n_sims(11);
    n_sims[0] = 1024;
    for (size_t i(1); i < n_sims.size(); i++) n_sims[i] = 2*n_sims[i-1];
    confianceIntervals(n_sims, linear_cos_lambda, f_bar, ABS_PATH + "/data/ps_1_2d_confidence_intervals.data", rng);



    // 3.a. Analytical computation
    double rate(0.05), sigma(0.2), maturity(1.0), initial_value(100.0), strike(100.0);
    f_bar = analytical_european_call(rate, sigma, maturity, initial_value, strike, "value");
    boost::math::normal normal_dist(0.0, 1.0);
    std::function<double(double)> payoff_european_call_lambda = [&](double x) -> double {
        return payoff_european_call(x, rate, sigma, maturity, initial_value, strike, normal_dist);
    };

    // 3.b. Monte-Carlo simulations
    for (size_t i(0); i < outer_sim; i++) {
        cdf[i] = empiricalMean(inner_sim, payoff_european_call_lambda, rng);
    }

    // 3.c. Sorting and plotting
    std::sort(cdf.begin(), cdf.end());
    std::ofstream outfile_3c;
    outfile_3c.open(ABS_PATH + "data/ps_1_3c_empirical_cdf.data");
    count = 0;
    for (const auto &e : cdf) {
        outfile_3c << e << " " << (count + 0.5) / outer_sim << " " << norm_cdf(sqrt(outer_sim) * (e-f_bar) / std_dev) << "\n";
        count++;
    }
    outfile_3c.close();


    // 3.d. Confidence intervals
    confianceIntervals(n_sims, payoff_european_call_lambda, f_bar,
                       ABS_PATH + "/data/ps_1_3d_confidence_intervals.data", rng);


    // 4.a. Antithetic variables
    double corr = antitheticVariables(n_sims, payoff_european_call_lambda,
                                      ABS_PATH + "/data/ps_1_4a_antithetic_variance.data", rng);

    cout << "Correlation between antithetic variables : " << corr << endl;

























    return 0;
}
