#ifndef UTILS
#define UTILS

#include <vector>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>

extern std::string ABS_PATH;

std::vector<double> operator+(std::vector<double> const&, std::vector<double> const&);

double mean(std::vector<double>);
double variance(std::vector<double>);

boost::numeric::ublas::matrix<double> mean(boost::numeric::ublas::matrix<double> m);
boost::numeric::ublas::matrix<double> covariance(boost::numeric::ublas::matrix<double> m);


#endif
