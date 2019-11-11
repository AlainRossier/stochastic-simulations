#include <numeric>
#include <functional>
#include <algorithm>
#include "utils.h"

#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace boost::numeric::ublas;

string ABS_PATH = "/home/alain/Documents/phd/lectures/stochastic_simulations/stochastic-simulations/";

double mean(std::vector<double> v)
{
    double sum = std::accumulate(begin(v), end(v), 0.0);
    return sum / v.size();
}


double variance(std::vector<double> v)
{
    double accum = 0.0;
    double m(mean(v));
    std::for_each (std::begin(v), std::end(v),
                   [&](const double d) {accum += (d - m) * (d - m);});
    return accum / (v.size()-1);
}


// Perform column-wise mean
matrix<double> mean(matrix<double> m)
{
    matrix<double> ones(1, m.size1(), 1);
    return prod(ones, m) / m.size1();
}

// Compute the covariance matrix
matrix<double> covariance(matrix<double> m)
{
    matrix<double> means = mean(m);
    return prod(trans(m), m) / m.size1() - prod(trans(means), means);
}
