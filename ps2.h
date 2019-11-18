#ifndef PS1
#define PS1

#include <string>
#include <random>


// Problem 2
void european_call_path(double, double, double, double, double, size_t, size_t,
                        std::string, std::default_random_engine&);

void european_call_path_2h(double, double, double, double, double, size_t, size_t,
                           std::string, std::default_random_engine&);

void gbm_strong_error(double, double, double, double, double, size_t, size_t,
                      std::string, std::default_random_engine&);

// Problem 3


// Problem 4

void heston_stochastic_vola_strong_error(double, double, double, double, double,
                                         double, double, double, size_t, size_t,
                                         std::string, std::default_random_engine&);

// Problem 5

void mcqmc06_l(int l, int N, double *sums, std::default_random_engine& rng, int option);


#endif

