#include <iostream>
#include <random>
#include <functional>
#include <algorithm>
#include <cmath>
#include <fstream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "utils.h"
#include "ps1.h"


using namespace std;
using namespace boost::numeric::ublas;


// Write the solution file


void solutions_ps1() {

    // Initialize the random number generator
    default_random_engine rng;

    // Initialize the output file
    ofstream out;
    out.open(ABS_PATH + "sols/ps1/solutions_ps1.txt");

    // Print up to the 5th most significant digit
    out.precision(5);
    out.setf(ios::fixed);
    out.setf(ios::showpoint);

    // 1.a. Random variable generation
    out << "1.a. Random variable generation" << endl;

    clock_t start;
    start = clock();

    size_t n_sim(1000000);

    std::vector<double> sample_uniform(n_sim);
    populate(sample_uniform, "uniform", rng);
    out << "CPU time elapsed to generate " << n_sim << " doubles : " << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms." << endl;
    out << "Empirical mean of the uniform distribution : " << mean(sample_uniform) << endl;
    out << "Empirical variance of the uniform distribution : " << variance(sample_uniform) << endl;


    out << "\nIf we are using floats, for n=10^7, the empirical variance is still quite far away";
    out << "from the theoretical one (order 10^(-2)). But it is faster (1.1s)." << endl;
    out << "Using double, represented by 64 bits, the absolute error on the variance is of order 10^(-3) at the expense of speed (1.4s)." << endl;

    std::vector<double> sample_normal(n_sim);
    populate(sample_normal, "normal", rng);

    out << "\nEmpirical mean of the normal distribution : " << mean(sample_normal) << endl ;
    out << "Empirical variance of the normal distribution : " << variance(sample_normal) << endl ;


    // 1.b. Cholesky decomposition
    out << "\n1.b. Cholesky decomposition" << endl;

    matrix<double> A(2, 2);
    A(0, 0) = 4; A(0, 1) = 1; A(1, 0) = 1; A(1, 1) = 4;
    matrix<double> Lchol = cholesky(A);
    out << "With Cholesky, L = " << Lchol << endl;

    std::vector<double> sample_normal_coor(2*n_sim);
    populate(sample_normal_coor, "normal", rng);
    // Copy the elements in a matrix
    unbounded_array<double> storage(2*n_sim);
    std::copy(sample_normal_coor.begin(), sample_normal_coor.end(), storage.begin());
    matrix<double> multivariate_normal(n_sim, 2, storage);
    matrix<double> corr_multivariate_normal_chol = prod(multivariate_normal, trans(Lchol));

    out << "Empirical mean of the multivariate normal distribution with Cholesky : " << mean(corr_multivariate_normal_chol) << endl ;
    out << "Empirical covariance of the multivariate normal distribution with Cholesky : " << covariance(corr_multivariate_normal_chol) << endl ;


    // 1.c. PCA factorization
    out << "\n1.c. PCA factorisation" << endl;

    matrix<double> Lpca = pca2d(A);
    out << "With PCA, L = " << Lpca << endl;

    matrix<double> corr_multivariate_normal_pca = prod(multivariate_normal, trans(Lpca));

    out << "Empirical mean of the multivariate normal distribution with PCA : " << mean(corr_multivariate_normal_pca) << endl ;
    out << "Empirical covariance of the multivariate normal distribution with PCA : " << covariance(corr_multivariate_normal_pca) << endl ;

    // 1.d. Run simulations
    out << "\n1.d. Run simulations" << endl;
    uniform_real_distribution<double> uniform(0.0f, 1.0f);
    auto next_uniform = std::bind(std::ref(uniform), std::ref(rng));
    start = clock();
    bool cont(true);
    size_t count_samples(0);
    while (cont) {
        next_uniform();
        count_samples++;
        if (count_samples % 10000000 == 0) {
            if ((clock() - start) / (double)(CLOCKS_PER_SEC) > 60) {
                cont = false;
            }
        }
    }

    out << "Number of uniform samples generated in 1 minute : " << count_samples << endl;
    out << "It seems to be on par with Python's numpy library. I wonder why I can't gain an edge over it." << endl;

    // 2.a. Analytical computation
    out << "\n\n2.a. Analytical computation" << endl;

    double f_bar = -2/pow(M_PI, 2);
    double std_dev = sqrt(1.0/6 - 4/pow(M_PI, 4) + 1.0/(4*pow(M_PI, 2)));
    out << "f_bar = " << f_bar << endl;
    out << "sigma = " << std_dev << endl;
    std::function<double(double)> linear_cos_lambda = [=](double x) -> double {return x*cos(M_PI*x);};

    // 2.b. Monte-Carlo simulations
    out << "\n2.b. Monte-Carlo simulations" << endl;

    size_t inner_sim(1000);
    size_t outer_sim(1000);
    std::vector<double> cdf(outer_sim);
    for (size_t i(0); i < outer_sim; i++) {
        cdf[i] = empiricalMean(inner_sim, linear_cos_lambda, rng);
    }

    // 2.c. Sorting and plotting
    out << "\n2.c. Sorting and plotting" << endl;

    std::sort(cdf.begin(), cdf.end());
    std::ofstream outfile_2c;
    outfile_2c.open(ABS_PATH + "data/ps_1_2c_empirical_cdf.data");
    size_t count(0);
    double max_error(0);
    double normal_cdf;
    for (const auto &e : cdf) {
        normal_cdf = norm_cdf(sqrt(outer_sim) * (e-f_bar) / std_dev);
        if (abs((count + 0.5) / outer_sim - normal_cdf) > max_error) {
            max_error = abs((count + 0.5) / outer_sim - normal_cdf);
        }
        outfile_2c << e << " " << (count + 0.5) / outer_sim << " " << normal_cdf << "\n";
        count++;
    }
    outfile_2c.close();

    out << "The CLT guarantees a convergence in distribution, meaning that the empirical CDF and the normal one converges pointwise everywhere." << endl;
    out << "The maximum error between the empirical CDF and the normal CDF is " << max_error << " for N = " << inner_sim << " samples." << endl;
    out << "We observed that the latter decays at a rate of 1/sqrt(n), which is guaranteed by the Berry–Esseen theorem." << endl;


    // 2.d. Confidence intervals
    out << "\n2.d. Confidence intervals" << endl;

    size_t n_experiments(11);
    std::vector<size_t> n_sims(n_experiments);
    std::vector<double> empirical_mean_2d(n_experiments);
    std::vector<double> empirical_sd_2d(n_experiments);
    n_sims[0] = 1024;
    for (size_t i(1); i < n_sims.size(); i++) n_sims[i] = 2*n_sims[i-1];
    confianceIntervals(n_sims, linear_cos_lambda, f_bar,
                       empirical_mean_2d, empirical_sd_2d,
                       ABS_PATH + "/data/ps_1_2d_confidence_intervals.data", rng);



    // 3.a. Analytical computation
    out << "\n\n3.a. Analytical computation" << endl;
    double rate(0.05), sigma(0.2), maturity(1.0), initial_value(100.0), strike(100.0);
    f_bar = analytical_european_call(rate, sigma, maturity, initial_value, strike, "value");
    boost::math::normal normal_dist(0.0, 1.0);
    std::function<double(double)> payoff_european_call_lambda = [&](double x) -> double {
        return payoff_european_call(quantile(normal_dist, x), rate, sigma, maturity, initial_value, strike);
    };

    out << "f_bar = " << f_bar << endl;

    // 3.b. Monte-Carlo simulations
    out << "\n3.b. Monte-Carlo simulations" << endl;
    for (size_t i(0); i < outer_sim; i++) {
        cdf[i] = empiricalMean(inner_sim, payoff_european_call_lambda, rng);
    }

    std_dev = sqrt(variance(cdf));

    // 3.c. Sorting and plotting
    out << "\n3.c. Sorting and plotting" << endl;
    std::sort(cdf.begin(), cdf.end());
    std::ofstream outfile_3c;
    outfile_3c.open(ABS_PATH + "data/ps_1_3c_empirical_cdf.data");
    count = 0;
    max_error = 0;
    for (const auto &e : cdf) {
        normal_cdf = norm_cdf((e-f_bar) / std_dev);
        if (abs((count + 0.5) / outer_sim - normal_cdf) > max_error) {
            max_error = abs((count + 0.5) / outer_sim - normal_cdf);
        }
        outfile_3c << e << " " << (count + 0.5) / outer_sim << " " << normal_cdf << "\n";
        count++;
    }
    outfile_3c.close();

    out << "The maximum error between the empirical CDF and the normal CDF is " << max_error << " for N = " << inner_sim << " samples." << endl;


    // 3.d. Confidence intervals
    out << "\n3.d. Confidence intervals" << endl;

    std::vector<double> empirical_mean_3d(n_experiments);
    std::vector<double> empirical_sd_3d(n_experiments);
    confianceIntervals(n_sims, payoff_european_call_lambda, f_bar,
                       empirical_mean_3d, empirical_sd_3d,
                       ABS_PATH + "/data/ps_1_3d_confidence_intervals.data", rng);


    // 4.a. Antithetic variables
    out << "\n\n4.a. Antithetic variables" << endl;
    double corr_anti = antitheticVariables(n_sims, payoff_european_call_lambda,
                                           ABS_PATH + "/data/ps_1_4a_antithetic_variance.data", rng);

    out << "Correlation between antithetic variables : " << corr_anti << endl;
    out << "Expected factor of variance reduction : " << 2/(1 + corr_anti) << endl;

    // 4.b. Control variate
    out << "\n\n4.b. Control variate" << endl;
    std::function<double(double)> control_variate_lambda = [&](double x) -> double {
        return payoff_european_call(quantile(normal_dist, x), rate, sigma, maturity, initial_value, 0);
    };

    double corr_control = controlVariate(n_sims, payoff_european_call_lambda,
                                         control_variate_lambda, initial_value,
                                         ABS_PATH + "/data/ps_1_4b_control_variate.data", rng);

    out << "Correlation between control variates : " << corr_control << endl;
    out << "Expected factor of variance reduction : " << 1/(1 - pow(corr_control, 2))  << endl;


    // 5.a. Digital put option
    out << "\n\n5.a. Digital put option" << endl;
    rate = 0.05; sigma = 0.2; maturity = 1.0; initial_value = 100.0; strike = 50.0;
    std::function<double(double)> payoff_digital_put_lambda = [&](double x) -> double {
        return payoff_digital_put(quantile(normal_dist, x), rate, sigma, maturity, initial_value, strike);
    };

    std::vector<double> empirical_mean_5a(n_experiments);
    std::vector<double> empirical_sd_5a(n_experiments);

    confianceIntervals(n_sims, payoff_digital_put_lambda, DBL_MAX,
                       empirical_mean_5a, empirical_sd_5a,
                       ABS_PATH + "/data/ps_1_5a_digital_put_value.data", rng);

    out << "Without importance sampling" << endl;
    for (size_t i(0); i < n_experiments; i++) {
        if (i >= 6) {
            out << "With " << n_sims[i] << " samples, the value is " << empirical_mean_5a[i] << " and is correct within ";
            out << 100*(empirical_sd_5a[i] / max(abs(empirical_mean_5a[i]), pow(10, -30))) << "%." << endl;
        }
    }

    out << "We see that roughly 1'000'000 samples are needed to have a 10% confidence bound." << endl;


    // 5.b. With importance sampling
    out << "\n5.b. With importance sampling" << endl;
    out << "We have p_1(x) = 1/sqrt(2*pi) * exp(-1/2*x^2), f(x) = exp(-rT) * H(K - S0*exp((r-1/2*sigma^2)*T + sigma*sqrt(T)*x))" << endl;
    out << "We let p_2(x) = 1/sqrt(2*pi) * exp(-1/2*(x+mu)^2), so that R(x) = exp(-mu*x-1/2*mu^2)" << endl;
    out << "We want to pick mu s.t. K = S0*exp((r-1/2*sigma^2)*T + sigma*sqrt(T)*mu)" << endl;

    double mu = (log(strike / initial_value) - (rate - 0.5*pow(sigma, 2))*maturity) / (sigma*sqrt(maturity));

    std::function<double(double)> payoff_digital_put_is_lambda = [&](double x) -> double {
        double normal = quantile(normal_dist, x);
        return payoff_digital_put(mu + normal, rate, sigma, maturity, initial_value, strike) * \
               exp(-mu*normal - 0.5*pow(mu, 2));
    };

    std::vector<double> empirical_mean_5b(n_experiments);
    std::vector<double> empirical_sd_5b(n_experiments);
    confianceIntervals(n_sims, payoff_digital_put_is_lambda, DBL_MAX,
                       empirical_mean_5b, empirical_sd_5b,
                       ABS_PATH + "/data/ps_1_5b_digital_put_value_importance_sampling.data", rng);

    out << "With importance sampling" << endl;
    for (size_t i(0); i < n_experiments; i++) {
        if (i < 3) {
            out << "With " << n_sims[i] << " samples, the value is " << empirical_mean_5b[i] <<  " and is correct within ";
            out << 100*(empirical_sd_5b[i] / max(abs(empirical_mean_5b[i]), pow(10, -30))) << "%." << endl;
        }
    }

    out << "We see that less than 1'000 samples are needed to have a 10% confidence bound, an improvement of 3 orders of magnitude." << endl;

    // 6.a. Bumping
    out << "\n\n6.a. Bumping" << endl;
    out << "We are using 100'000 samples to compute the finite differences." << endl;
    rate = 0.05; sigma = 0.2; maturity = 1.0; initial_value = 100.0; strike = 100.0;
    std::function<double(double, double)> delta_european_call_bump = [&](double x, double bump) -> double {
        return payoff_european_call(quantile(normal_dist, x), rate, sigma, maturity, initial_value*exp(bump), strike) / initial_value;
    };
    std::function<double(double, double)> vega_european_call_bump = [&](double x, double bump) -> double {
        return payoff_european_call(quantile(normal_dist, x), rate, sigma*exp(bump), maturity, initial_value, strike) / sigma;
    };

    std::vector<double> bumps(n_experiments);
    for (int i(0); i < n_experiments; i++) {bumps[i] = pow(2, -i);}

    double true_delta = analytical_european_call(rate, sigma, maturity, initial_value, strike, "delta");
    double true_vega = analytical_european_call(rate, sigma, maturity, initial_value, strike, "vega");

    finite_difference_bumping(bumps, 100000, delta_european_call_bump, true_delta,
                              ABS_PATH + "/data/ps_1_6a_bumping_delta_european_call.data", rng);

    finite_difference_bumping(bumps, 100000, vega_european_call_bump, true_vega,
                              ABS_PATH + "/data/ps_1_6a_bumping_vega_european_call.data", rng);


    // 6.b. Pathwise sensitivity method
    out << "\n6.b. Pathwise sensitivity method" << endl;
    out << "For the delta, we have f(S0, Z) = exp(-rT) * max(S0*exp((r-1/2*sigma^2)*T + sigma*sqrt(T)*Z) - K, 0)" << endl;
    out << "so that del(f)/del(S0) = exp(-rT) * exp((r-1/2*sigma^2)*T + sigma*sqrt(T)*Z) * I(S0*exp((r-1/2*sigma^2)*T + sigma*sqrt(T)*Z) > K)" << endl;

    std::function<double(double)> delta_european_call_ipa = [&](double x) -> double {
        double normal = quantile(normal_dist, x);
        double price = initial_value * exp((rate - 0.5*pow(sigma, 2)) * maturity + sigma * sqrt(maturity) * normal);
        return exp(-rate * maturity) * price / initial_value * (price > strike);
    };

    std::vector<double> empirical_mean_6b_delta(n_experiments);
    std::vector<double> empirical_sd_6b_delta(n_experiments);

    confianceIntervals(n_sims, delta_european_call_ipa, true_delta,
                       empirical_mean_6b_delta, empirical_sd_6b_delta,
                       ABS_PATH + "/data/ps_1_6b_ipa_delta_european_call.data", rng);

    out << "\nFor the vega, we have del(f)/del(sigma) = exp(-rT) * S0*exp((r-1/2*sigma^2)*T + sigma*sqrt(T)*Z) *" << endl;
    out << "(-sigma*T + sqrt(T)*Z) * I(S0*exp((r-1/2*sigma^2)*T + sigma*Z) > K)" << endl;

    std::function<double(double)> vega_european_call_ipa = [&](double x) -> double {
        double normal = quantile(normal_dist, x);
        double price = initial_value * exp((rate - 0.5*pow(sigma, 2)) * maturity + sigma * sqrt(maturity) * normal);
        return exp(-rate * maturity) * price * (-sigma * maturity + sqrt(maturity) * normal) * (price > strike);
    };

    std::vector<double> empirical_mean_6b_vega(n_experiments);
    std::vector<double> empirical_sd_6b_vega(n_experiments);

    confianceIntervals(n_sims, vega_european_call_ipa, true_vega,
                       empirical_mean_6b_vega, empirical_sd_6b_vega,
                       ABS_PATH + "/data/ps_1_6b_ipa_vega_european_call.data", rng);

    out.close();

}

