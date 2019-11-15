#include <numeric>
#include <functional>
#include <algorithm>
#include "utils.h"

#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace boost::numeric::ublas;

string ABS_PATH = "/home/alain/Documents/phd/lectures/stochastic_simulations/stochastic-simulations/";



std::vector<double> operator+(std::vector<double> const& v , std::vector<double> const& w) {
    if (v.size() != w.size()) {
        cerr << "You're shit." << endl;
        return std::vector<double>(0);

    } else {
        size_t n = v.size();
        std::vector<double> output(n);
        for (size_t i(0); i < n; ++i) {
            output[i] = v[i] + w[i];
        }
        return output;
    }
}


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
