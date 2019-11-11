#ifndef PS1
#define PS1

#include <string>
#include <random>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/distributions/normal.hpp>


// Problem 1
void populate(std::vector<double>&, std::string, std::default_random_engine &);
boost::numeric::ublas::matrix<double> cholesky(boost::numeric::ublas::matrix<double>);
boost::numeric::ublas::matrix<double> pca2d(boost::numeric::ublas::matrix<double>);


// Problem 2
double linear_cos(double);
double empiricalMean(size_t, const std::function<double(double)>&, std::default_random_engine &);
double norm_cdf(double);
void confianceIntervals(std::vector<size_t>, const std::function<double(double)>&, double,
                        std::vector<double>&, std::vector<double>&,
                        std::string, std::default_random_engine &);

// Problem 3
double analytical_european_call(double, double, double, double, double, std::string);
double payoff_european_call(double, double, double, double, double, double);

// Problem 4
double antitheticVariables(std::vector<size_t>, const std::function<double(double)>&,
                           std::string, std::default_random_engine&);
double controlVariate(std::vector<size_t>, const std::function<double(double)>&,
                      const std::function<double(double)>&, double,
                      std::string, std::default_random_engine&);

// Problem 5
double payoff_digital_put(double, double, double, double, double, double);

// Problem 6
void finite_difference_bumping(std::vector<double>, size_t, const std::function<double(double, double)>&,
                               double, std::string, std::default_random_engine&);


#endif
