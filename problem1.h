#ifndef PROBLEM1
#define PROBLEM1

#include <string>
#include <random>
#include <boost/numeric/ublas/matrix.hpp>

// Problem 1
void populate(std::vector<double>&, std::string, std::default_random_engine &);
boost::numeric::ublas::matrix<double> cholesky(boost::numeric::ublas::matrix<double>);

// Problem 2
double f(double);
double empiricalMean(size_t, std::default_random_engine &);
void confianceIntervals(std::vector<size_t>, std::default_random_engine &);


#endif
