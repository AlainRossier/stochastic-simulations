#ifndef HELPERS
#define HELPERS

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>


double mean(std::vector<double>);
double variance(std::vector<double>);

boost::numeric::ublas::matrix<double> mean(boost::numeric::ublas::matrix<double> m);
boost::numeric::ublas::matrix<double> covariance(boost::numeric::ublas::matrix<double> m);


#endif
