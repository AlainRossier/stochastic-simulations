#include <iostream>
#include <random>
#include <functional>
#include <algorithm>
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>


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
double f(double x) {
    return x * cos(M_PI * x);
}

double empiricalMean(size_t n_sim, default_random_engine& rng) {
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

void confianceIntervals(std::vector<size_t> n_sims, std::default_random_engine& rng) {
    // Initializations
    size_t n_experiments = n_sims.size();
    std::vector<double> empirical_means(n_experiments);
    std::vector<double> empirical_sd(n_experiments);
    std::sort(n_sims.begin(), n_sims.end());
    size_t max_sim = n_sims.back();
    double sum(0.0), sumsq(0.1), acc(0.0);
    size_t incr(0), k(0);

    // Generation
    uniform_real_distribution<double> uniform(0.0f, 1.0f);
    auto next_uniform = bind(ref(uniform), ref(rng));
    while (incr < max_sim) {
        acc = f(next_uniform());
        sum += acc; sumsq += pow(acc, 2);
        incr += 1;
        if (incr == n_sims[k]) {
            empirical_means[k] = sum / incr;
            empirical_sd[k] = sqrt((sumsq - incr * empirical_means[k]) / (incr-1));
            k += 1;
        }
    }

}

