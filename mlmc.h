#ifndef MLMC
#define MLMC

#include <fstream>
#include <random>

float mlmc(int Lmin, int Lmax, int N0, float eps,
           void (*mlmc_l)(int, int, double *, std::default_random_engine&, int),
           float alpha_0, float beta_0, float gamma_0, int *Nl, float *Cl, std::default_random_engine&, int args);

void mlmc_test(void (*mlmc_l)(int, int, double *, std::default_random_engine&, int),
               int M,int N,int L, int N0, float *Eps, int Lmin, int Lmax,
               std::ofstream&, std::default_random_engine&, int args);


#endif
