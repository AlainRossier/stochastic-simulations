#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <time.h>
#include <random>

#include "utils.h"


using namespace std;


/*
   P = mlmc(Lmin,Lmax,N0,eps, mlmc_l, alpha,beta,gamma, Nl,Cl)

   multilevel Monte Carlo control routine

   Lmin  = minimum level of refinement       >= 2
   Lmax  = maximum level of refinement       >= Lmin
   N0    = initial number of samples         > 0
   eps   = desired accuracy (rms error)      > 0

   mlmc_l(l,N,sums)   low-level function
        l       = level
        N       = number of paths
        sums[0] = sum(cost)
        sums[1] = sum(Y)
        sums[2] = sum(Y.^2)
        where Y are iid samples with expected value:
        E[P_0]           on level 0
        E[P_l - P_{l-1}] on level l>0

   alpha -> weak error is  O(2^{-alpha*l})
   beta  -> variance is    O(2^{-beta*l})
   gamma -> sample cost is O(2^{gamma*l})

   if alpha, beta, gamma are not positive then they will be estimated

   P   = value
   Nl  = number of samples at each level
   Cl  = average cost of samples at each level

   rng = random number generator
   args = additional arguments
*/


float mlmc(int Lmin, int Lmax, int N0, float eps,
           void (*mlmc_l)(int, int, double *, std::default_random_engine&, int),
           float alpha_0, float beta_0, float gamma_0,
           int *Nl, float *Cl, std::default_random_engine& rng, int args) {

  double sums[7], suml[3][21];
  float  ml[21], Vl[21], NlCl[21], x[21], y[21],
         alpha, beta, gamma, sum, theta;
  int    dNl[21], L, converged;

  int    diag = 0;  // diagnostics, set to 0 for none

  //
  // check input parameters
  //

  if (Lmin<2) {
    cerr << "error: needs Lmin >= 2" << endl;
    exit(1);
  }
  if (Lmax<Lmin) {
    cerr << "error: needs Lmax >= Lmin" << endl;
    exit(1);
  }

  if (N0<=0 || eps<=0.0f) {
    cerr << "error: needs N>0, eps>0" << endl;
    exit(1);
  }

  //
  // initialisation
  //

  alpha = fmax(0.0f,alpha_0);
  beta  = fmax(0.0f,beta_0);
  gamma = fmax(0.0f,gamma_0);
  theta = 0.25f;             // MSE split between bias^2 and variance

  L = Lmin;
  converged = 0;

  for(int l=0; l<=Lmax; l++) {
    Nl[l]   = 0;
    Cl[l]   = powf(2.0f,(float)l*gamma);
    NlCl[l] = 0.0f;

    for(int n=0; n<3; n++) suml[n][l] = 0.0;
  }

  for(int l=0; l<=Lmin; l++) dNl[l] = N0;

  //
  // main loop
  //

  while (!converged) {

    //
    // update sample sums
    //

    for (int l=0; l<=L; l++) {
      if (diag) {cout << dNl[l] << endl;}

      if (dNl[l]>0) {
        mlmc_l(l,dNl[l],sums, rng, args);
        suml[0][l] += (float) dNl[l];
        suml[1][l] += sums[1];
        suml[2][l] += sums[2];
        NlCl[l]    += sums[0];  // sum total cost
      }
    }
    if (diag) printf(" \n");

    //
    // compute absolute average, variance and cost,
    // correct for possible under-sampling,
    // and set optimal number of new samples
    //

    sum = 0.0f;

    for (int l=0; l<=L; l++) {
      ml[l] = fabs(suml[1][l]/suml[0][l]);
      Vl[l] = fmaxf(suml[2][l]/suml[0][l] - ml[l]*ml[l], 0.0f);
      if (gamma_0 <= 0.0f) Cl[l] = NlCl[l] / suml[0][l];

      if (l>1) {
        ml[l] = fmaxf(ml[l],  0.5f*ml[l-1]/powf(2.0f,alpha));
        Vl[l] = fmaxf(Vl[l],  0.5f*Vl[l-1]/powf(2.0f,beta));
      }

      sum += sqrtf(Vl[l]*Cl[l]);
    }

    for (int l=0; l<=L; l++) {
      dNl[l] = ceilf( fmaxf( 0.0f,
                       sqrtf(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)
                     - suml[0][l] ) );
    }

    //
    // use linear regression to estimate alpha, beta, gamma if not given
    //

    if (alpha_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(ml[l]);
      }
      regression(L,x,y,alpha,sum);
      if (diag) {cout << " alpha = " << alpha << endl;}
    }

    if (beta_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(Vl[l]);
      }
      regression(L,x,y,beta,sum);
      if (diag) {cout << " beta = " << beta << endl;}
    }

     if (gamma_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = log2f(Cl[l]);
      }
      regression(L,x,y,gamma,sum);
      if (diag) {cout << " gamma = " << gamma << endl;}
    }

    //
    // if (almost) converged, estimate remaining error and decide
    // whether a new level is required
    //

    sum = 0.0;
      for (int l=0; l<=L; l++)
        sum += fmaxf(0.0f, (float)dNl[l]-0.01f*suml[0][l]);

    if (sum==0) {
      if (diag) cout << "Achieved variance target" << endl;

      converged = 1;
      float rem = ml[L] / (powf(2.0f,alpha)-1.0f);

      if (rem > sqrtf(theta)*eps) {
        if (L==Lmax) {
          cout << "*** failed to achieve weak convergence ***" << endl;
        } else {
          converged = 0;
          L++;
          Vl[L] = Vl[L-1]/powf(2.0f,beta);
          Cl[L] = Cl[L-1]*powf(2.0f,gamma);

          if (diag) {cout << "L = " << L << endl;}

          sum = 0.0f;
          for (int l=0; l<=L; l++) sum += sqrtf(Vl[l]*Cl[l]);
          for (int l=0; l<=L; l++)
            dNl[l] = ceilf( fmaxf( 0.0f,
                            sqrtf(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)
                          - suml[0][l] ) );
        }
      }
    }
  }

  //
  // finally, evaluate multilevel estimator and set outputs
  //

  float P = 0.0f;
  for (int l=0; l<=L; l++) {
    P    += suml[1][l]/suml[0][l];
    Nl[l] = suml[0][l];
    Cl[l] = NlCl[l] / Nl[l];
  }

  return P;
}



/*

  mlmc_test(mlmc_l, M, N,L, N0,Eps,Lmin,Lmax, out)

  multilevel Monte Carlo test routine

   mlmc_l(l,N,sums)     low-level routine

   inputs:  l = level
            N = number of paths

   output: sums[0] = sum(cost)
           sums[1] = sum(Pf-Pc)
           sums[2] = sum((Pf-Pc).^2)
           sums[3] = sum((Pf-Pc).^3)
           sums[4] = sum((Pf-Pc).^4)
           sums[5] = sum(Pf)
           sums[6] = sum(Pf.^2)

   M      = refinement cost factor (if zero, gamma is estimated)

   N      = number of samples for convergence tests
   L      = number of levels for convergence tests

   N0     = initial number of samples
   Eps    = desired accuracy array (terminated by value 0)
   Lmin   = minimum level of refinement
   Lmax   = maximum level of refinement

   out     = file handle for output
*/




void mlmc_test(void (*mlmc_l)(int, int, double *, std::default_random_engine&, int),
               int M,int N,int L, int N0, float *Eps, int Lmin, int Lmax,
               std::ofstream &out, std::default_random_engine& rng, int args) {

    // first, convergence tests

    // current date/time based on current system
    time_t now = time(NULL);
    char *date = ctime(&now);

    out << "**********************************************************" << endl;
    out << "*** MLMC file version 0.9     produced by              ***" << endl;
    out << "*** C++ mlmc_test on " << date << "         ***" << endl;
    out << "**********************************************************" << endl;
    out << endl;
    out << "**********************************************************" << endl;
    out << "*** Convergence tests, kurtosis, telescoping sum check ***" << endl;
    out << "*** using N = " << N << " samples                 ***\n" << endl;
    out << "**********************************************************" << endl;
    out << "l   ave(Pf-Pc)    ave(Pf)   var(Pf-Pc)    var(Pf)" << endl;
    out << "    kurtosis     check        cost \n--------------------------";
    out << "-------------------------------------------------------------" << endl;

    double sums[7];
    float *cost = (float *)malloc((L+1)*sizeof(float));
    float *del1 = (float *)malloc((L+1)*sizeof(float));
    float *del2 = (float *)malloc((L+1)*sizeof(float));
    float *var1 = (float *)malloc((L+1)*sizeof(float));
    float *var2 = (float *)malloc((L+1)*sizeof(float));
    float *chk1 = (float *)malloc((L+1)*sizeof(float));
    float *kur1 = (float *)malloc((L+1)*sizeof(float));

    for (int l=0; l<=L; l++) {
        mlmc_l(l,N,sums,rng,args);

        for (int m=0; m<7; m++) {sums[m] = sums[m]/N;}

        if (M>0) {cost[l] = powf((float)M,(float)l);}
        else {cost[l] = sums[0];}

        del1[l] = sums[1];
        del2[l] = sums[5];
        var1[l] = fmax(sums[2]-sums[1]*sums[1], 1e-10);
        var2[l] = fmax(sums[6]-sums[5]*sums[5], 1e-10);

        kur1[l]  = (sums[4] - 4.0*sums[3]*sums[1] + 6.0*sums[2]*sums[1]*sums[1]
                    - 3.0*sums[1]*sums[1]*sums[1]*sums[1]) / (var1[l]*var1[l]);

        if (l==0) {
            chk1[l] = 0.0f;
        } else {
            chk1[l] = sqrtf((float) N) * fabsf(del1[l] + del2[l-1] - del2[l])
            / (3.0f*(sqrtf(var1[l]) + sqrtf(var2[l-1]) + sqrtf(var2[l])));

            out << l << del1[l] << del2[l] << var1[l] << var2[l] << kur1[l] << chk1[l] << cost[l] << endl;
        }

        // print out a warning if kurtosis or consistency check looks bad

        if (kur1[L] > 100.0f) {
            out << "\n WARNING: kurtosis on finest level = " << kur1[L] << endl;
            out << " indicates MLMC correction dominated by a few rare paths; \n";
            out << " for information on the connection to variance of sample variances,\n";
            out << " see http://mathworld.wolfram.com/SampleVarianceDistribution.html \n";
        }

        float max_chk = 0.0f;
        for (int l=0; l<=L; l++) max_chk = fmaxf(max_chk,chk1[l]);
        if (max_chk > 1.0f) {
            out << "\n WARNING: maximum consistency error = " << max_chk << endl;
            out << " indicates identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied \n";
        }

        // use linear regression to estimate alpha, beta, gamma

        float alpha, beta, gamma, foo;
        float *x = (float *)malloc(L*sizeof(float));
        float *y = (float *)malloc(L*sizeof(float));

        for (int l=1; l<=L; l++) {
            x[l-1] = l;
            y[l-1] = - log2f(fabsf(del1[l]));
        }
        regression(L,x,y,alpha,foo);

        for (int l=1; l<=L; l++) {
            x[l-1] = l;
            y[l-1] = - log2f(var1[l]);
        }
        regression(L,x,y,beta,foo);

        for (int l=1; l<=L; l++) {
            x[l-1] = l;
            y[l-1] = log2f(cost[l]);
        }
        regression(L,x,y,gamma,foo);

        out << "\n******************************************************\n";
        out << "*** Linear regression estimates of MLMC parameters ***\n";
        out << "******************************************************\n";
        out << "\n alpha = " << alpha << " (exponent for MLMC weak convergence)\n";
        out << " beta  = " << beta << " (exponent for MLMC variance) \n";
        out << " gamma = " << gamma << " (exponent for MLMC cost) \n";

        // Second, mlmc complexity tests

        out << endl;
        out << "***************************** \n";
        out << "*** MLMC complexity tests *** \n";
        out << "***************************** \n\n";
        out << "  eps       value   mlmc_cost   std_cost  savings     N_l \n";
        out << "--------------------------------------------------------- \n";

        int i=0;
        int   *Nl = (int *)malloc((Lmax+1)*sizeof(int));
        float *Cl = (float *)malloc((Lmax+1)*sizeof(float));

        while (Eps[i]>0) {
            float eps = Eps[i++];

            float P = mlmc(Lmin,Lmax,N0,eps,mlmc_l, alpha,beta,gamma, Nl,Cl,rng, args);

            float std_cost = 0.0f, mlmc_cost = 0.0f, theta=0.25f;

            for (int l=0; l<=Lmax; l++) {
                if (Nl[l]>0) {
                    // printf(" l, Cl, cost = %d  %f  %f \n",l,Cl[l],cost[l]);
                    mlmc_cost += Nl[l]*Cl[l];
                    if (l<=L) {
                        std_cost = var2[l]*cost[l] / ((1.0f-theta)*eps*eps);
                    } else {
                        std_cost = var2[L]*Cl[l] / ((1.0f-theta)*eps*eps);
                    }
                }
            }

            out << eps << P << mlmc_cost << std_cost << std_cost/mlmc_cost << endl;
            for (int l=0; Nl[l]>0; l++) {out << Nl[l] << endl;}
        }

        out << endl;
    }
}
