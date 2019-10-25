#include <iostream>
#include <algorithm>
#include <random>
#include "helpers.h"
#include "problem1.h"
#include <boost/numeric/ublas/io.hpp>


using namespace std;
using namespace boost::numeric::ublas;


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
    // f_bar = -2/pi^2 = -0.20264
    // sigma^2 = 1/6 - 4/pi^4 + 1/(4*pi^2)

    // 2.b. Monte-Carlo simulations
    size_t inner_sim(1000);
    size_t outer_sim(1000);
    std::vector<double> cdf(outer_sim);
    for (size_t i(0); i < outer_sim; i++) {
        cdf[i] = empiricalMean(inner_sim, rng);
    }
    cout << cdf[1] << cdf[2] << cdf[3] << endl;

    // 2.c. Sorting and plotting
    std::sort(cdf.begin(), cdf.end());

    // 2.d. Confidence intervals
    std::vector<size_t> n_sims = {1000, 10000, 100000, 1000000};






















    return 0;
}
