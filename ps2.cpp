#include <iostream>
#include <cmath>
#include <fstream>

#include "ps1.h"
#include "utils.h"

#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace boost::numeric::ublas;


// Problem 2

void european_call_path(double rate, double sigma, double maturity, double initial_value, double strike,
                        size_t n_paths, size_t n_levels,
                        string path, default_random_engine& rng) {

    // Initialize
    size_t incr(1);
    double step(0.0);
    double sum, sumsq, price, dW, final_price, mean, sd;

    // Analytical price
    double value_european_call = analytical_european_call(rate, sigma, maturity, initial_value, strike, "value");

    // Random number generator
    normal_distribution<double> normal(0.0f, 1.0f);
    auto next_normal = bind(ref(normal), ref(rng));

    // Saving
    std::ofstream outfile;
    outfile.open(path);

    for (size_t p(1); p <= n_levels; p++) {
        incr *= 2;
        step = maturity / incr;
        sum = 0.0; sumsq = 0.0;
        for (size_t m(0); m < n_paths; m++) {
            price = initial_value;
            for (size_t i(0); i < incr; i++) {
                dW = sqrt(step) * next_normal();
                price *= (1 + rate*step + sigma*dW);
            }
            final_price = exp(-rate*maturity) * max(price-strike, 0.0);
            sum += final_price;
            sumsq += pow(final_price, 2);
        }

        mean = sum/n_paths;
        sd = sqrt((sumsq/n_paths - pow(mean, 2)) / (n_paths-1));
        outfile << step << " " << abs(mean - value_european_call) << " " << 3*sd << "\n";
    }

    outfile.close();

}


void european_call_path_2h(double rate, double sigma, double maturity, double initial_value, double strike,
                           size_t n_paths, size_t n_levels,
                           string path, default_random_engine& rng) {

    // Initialize
    size_t incr(1);
    double step(0.0);
    double sum, sumsq, price, price_2h, dW, dW_2h, final_price, final_price_2h, mean, sd;

    // Analytical price
    double value_european_call = analytical_european_call(rate, sigma, maturity, initial_value, strike, "value");

    // Random number generator
    normal_distribution<double> normal(0.0f, 1.0f);
    auto next_normal = bind(ref(normal), ref(rng));

    // Saving
    std::ofstream outfile;
    outfile.open(path);

    for (size_t p(1); p <= n_levels; p++) {
        incr *= 2;
        step = maturity / incr;
        sum = 0.0; sumsq = 0.0;
        for (size_t m(0); m < n_paths; m++) {
            price = initial_value;
            price_2h = initial_value;
            for (size_t i(0); i < incr; i++) {
                dW = sqrt(step) * next_normal();
                dW_2h = sqrt(step) * next_normal();
                price *= ((1 + rate*step + sigma*dW) * (1 + rate*step + sigma*dW_2h));
                price_2h *= (1 + rate*2*step + sigma*(dW + dW_2h));
            }
            final_price = exp(-rate*maturity) * max(price-strike, 0.0);
            final_price_2h = exp(-rate*maturity) * max(price_2h-strike, 0.0);

            sum += (final_price - final_price_2h);
            sumsq += pow(final_price - final_price_2h, 2);
        }

        mean = sum/n_paths;
        sd = sqrt((sumsq/n_paths - pow(mean, 2)) / (n_paths-1));
        outfile << step << " " << abs(mean) << " " << 3*sd << "\n";
    }

    outfile.close();

}


void gbm_strong_error(double rate, double sigma, double maturity, double initial_value, double strike,
                      size_t n_paths, size_t n_levels,
                      string path, default_random_engine& rng) {

    // Initialize
    size_t incr(1);
    double step(0.0);
    double sum, sumsq, sum_2h, sumsq_2h;
    double price, price_2h, dW1, dW2, brownian;
    double exact_price, err, err_2h;
    double mean, sd, mean_2h, sd_2h;

    // Random number generator
    normal_distribution<double> normal(0.0f, 1.0f);
    auto next_normal = bind(ref(normal), ref(rng));

    // Saving
    std::ofstream outfile;
    outfile.open(path);

    for (size_t p(1); p <= n_levels; p++) {
        incr *= 2;
        step = maturity / incr;
        sum = 0.0; sumsq = 0.0; sum_2h = 0.0; sumsq_2h = 0.0;

        for (size_t m(0); m < n_paths; m++) {
            price = initial_value;
            price_2h = initial_value;
            brownian = 0.0;

            for (size_t i(0); i < incr/2; i++) {
                dW1 = sqrt(step) * next_normal();
                dW2 = sqrt(step) * next_normal();
                price *= ((1 + rate*step + sigma*dW1) * (1 + rate*step + sigma*dW2));
                price_2h *= (1 + rate*2*step + sigma*(dW1+dW2));
                brownian += (dW1+dW2);
            }
            exact_price = initial_value * exp((rate - 0.5*pow(sigma, 2))*maturity + sigma*brownian);

            err = pow(price - exact_price, 2);
            sum += err; sumsq += pow(err, 2);

            err_2h = pow(price_2h - exact_price, 2);
            sum_2h += err_2h; sumsq_2h += pow(err_2h, 2);
        }

        mean = sum/n_paths;
        sd = sqrt((sumsq/n_paths - pow(mean, 2)) / (n_paths-1));
        mean_2h = sum_2h/n_paths;
        sd_2h = sqrt((sumsq_2h/n_paths - pow(mean_2h, 2)) / (n_paths-1));

        outfile << step << " " << sqrt(mean) << " " << (0.5/sqrt(mean)) * 3*sd << " ";
        outfile << sqrt(mean_2h) << " " << (0.5/sqrt(mean_2h)) * 3*sd_2h << endl;
    }

}


// Problem 3



// Problem 4

void heston_stochastic_vola_strong_error(double rate, double initial_vola, double maturity, double initial_value,
                                         double theta, double kappa, double zeta, double rho,
                                         size_t n_paths, size_t n_levels,
                                         string path, default_random_engine& rng) {

    // Initialize
    size_t incr(1);
    double step(0.0);
    double sum, sumsq;
    double price_h, price_2h, vola_h, vola_2h;
    double norm_0, norm_1, dW1_1, dW1_2, dW2_1, dW2_2;
    double err, mean, sd;

    // Covariance transform
    matrix<double> Sigma(2, 2);
    Sigma(0, 0) = 1; Sigma(0, 1) = -0.1; Sigma(1, 0) = -0.1; Sigma(1, 1) = 1;
    matrix<double> L = cholesky(Sigma);

    // Random number generator
    normal_distribution<double> normal(0.0f, 1.0f);
    auto next_normal = bind(ref(normal), ref(rng));

    // Saving
    std::ofstream outfile;
    outfile.open(path);

    for (size_t p(1); p <= n_levels; p++) {
        incr *= 2;
        step = maturity / incr;
        sum = 0.0; sumsq = 0.0;

        for (size_t m(0); m < n_paths; m++) {
            price_h = initial_value;
            price_2h = initial_value;
            vola_h = initial_vola;
            vola_2h = initial_vola;

            for (size_t i(0); i < incr/2; i++) {
                // Draw correlated normals
                norm_0 = next_normal(); norm_1 = next_normal();
                dW1_1 = sqrt(step) * (L(0,0) * norm_0  + L(0,1) * norm_1) ;
                dW1_2 = sqrt(step) * (L(1,0) * norm_0  + L(1,1) * norm_1) ;
                norm_0 = next_normal(); norm_1 = next_normal();
                dW2_1 = sqrt(step) * (L(0,0) * norm_0  + L(0,1) * norm_1) ;
                dW2_2 = sqrt(step) * (L(1,0) * norm_0  + L(1,1) * norm_1) ;

                // Update path h
                vola_h += kappa*(theta - vola_h)*step + zeta*sqrt(abs(vola_h))*dW1_2;
                price_h *= (1 + rate*step + sqrt(abs(vola_h))*dW1_1);
                vola_h += kappa*(theta - vola_h)*step + zeta*sqrt(abs(vola_h))*dW2_2;
                price_h *= (1 + rate*step + sqrt(abs(vola_h))*dW2_1);

                // Update path 2h
                vola_2h += kappa*(theta - vola_2h)*2*step + zeta*sqrt(abs(vola_2h))*(dW1_2 + dW2_2);
                price_2h *= (1 + rate*2*step + sqrt(abs(vola_2h))*(dW1_1 + dW2_1));
            }

            err = pow(price_2h - price_h, 2);
            sum += err; sumsq += pow(err, 2);
        }

        mean = sum/n_paths;
        sd = sqrt((sumsq/n_paths - pow(mean, 2)) / (n_paths-1));

        outfile << step << " " << sqrt(mean) << " " << (0.5/sqrt(mean)) * 3*sd << endl;
    }

}



// Problem 5


// level l estimator

void mcqmc06_l(int l, int N, double *sums, std::default_random_engine& rng, int option) {

    // Define generators
    std::normal_distribution<float> normal(0.0f,1.0f);
    std::lognormal_distribution<float> lognormal(0.0f,1.0f);
    auto next_normal = std::bind(std::ref(normal), std::ref(rng));
    auto next_lognormal = std::bind(std::ref(lognormal), std::ref(rng));

    int M, nf, nc;
    float T, r, sig, B, hf, hc, X0, Xf, Xc, Af, Ac, Mf, Mc, Bf, Bc,
    Xf0, Xc0, Xc1, vf, vc, dWc, ddW, Pf, Pc, dP, K;

    float dWf[2], dIf[2], Lf[2];

    // model parameters

    K   = 100.0f;
    T   = 1.0f;
    r   = 0.05f;
    sig = 0.2f;
    B   = 0.85f*K;

    nf = 1<<l;
    nc = nf/2;

    hf = T / ((float) nf);
    hc = T / ((float) nc);

    for (int k=0; k<7; k++) sums[k] = 0.0;

    for (int np = 0; np<N; np++) {
        X0 = K;
        Xf = X0;
        Xc = Xf;

        Af  = 0.5f*hf*Xf;
        Ac  = 0.5f*hc*Xc;

        Mf  = Xf;
        Mc  = Xc;

        Bf  = 1.0f;
        Bc  = 1.0f;

        if (l==0) {
            dWf[0] = sqrt(hf) * next_normal();
            Lf[0]  = - next_lognormal();
            dIf[0] = sqrt(hf/12.0f)*hf * next_normal();

            Xf0 = Xf;
            Xf  = Xf + r*Xf*hf + sig*Xf*dWf[0]
                    + 0.5f*sig*sig*Xf*(dWf[0]*dWf[0]-hf);
            vf  = sig*Xf0;
            Af  = Af + 0.5f*hf*Xf + vf*dIf[0];
            Mf  = fminf(Mf, 0.5f*(Xf0+Xf-sqrtf((Xf-Xf0)*(Xf-Xf0)-2.0f*hf*vf*vf*Lf[0])));
            Bf  = Bf*(1.0f-expf(-2.0f*fmaxf(0.0f,(Xf0-B)*(Xf-B)/(hf*vf*vf))));

        } else {
            for (int n=0; n<nc; n++) {
                dWf[0] = sqrt(hf) * next_normal();
                dWf[1] = sqrt(hf) * next_normal();
                Lf[0]  = - next_lognormal();
                Lf[1]  = - next_lognormal();
                dIf[0] = sqrt(hf/12.0f)*hf * next_normal();
                dIf[1] = sqrt(hf/12.0f)*hf * next_normal();

                for (int m=0; m<2; m++) {
                    Xf0 = Xf;
                    Xf  = Xf + r*Xf*hf + sig*Xf*dWf[m]
                    + 0.5f*sig*sig*Xf*(dWf[m]*dWf[m]-hf);
                    vf  = sig*Xf0;
                    Af  = Af + hf*Xf + vf*dIf[m];
                    Mf  = fminf(Mf,
                            0.5f*(Xf0+Xf-sqrtf((Xf-Xf0)*(Xf-Xf0)-2.0f*hf*vf*vf*Lf[m])));
                    Bf  = Bf*(1.0f-expf(-2.0f*fmaxf(0.0f,(Xf0-B)*(Xf-B)/(hf*vf*vf))));
                }

                dWc = dWf[0] + dWf[1];
                ddW = dWf[0] - dWf[1];

                Xc0 = Xc;
                Xc  = Xc + r*Xc*hc + sig*Xc*dWc + 0.5f*sig*sig*Xc*(dWc*dWc-hc);

                vc  = sig*Xc0;
                Ac  = Ac + hc*Xc + vc*(dIf[0]+dIf[1] + 0.25f*hc*ddW);
                Xc1 = 0.5f*(Xc0 + Xc + vc*ddW);
                Mc  = fminf(Mc, 0.5f*(Xc0+Xc1-sqrtf((Xc1-Xc0)*(Xc1-Xc0)-2.0f*hf*vc*vc*Lf[0])));
                Mc  = fminf(Mc, 0.5f*(Xc1+Xc -sqrtf((Xc -Xc1)*(Xc -Xc1)-2.0f*hf*vc*vc*Lf[1])));
                Bc  = Bc *(1.0f-expf(-2.0f*fmaxf(0.0f,(Xc0-B)*(Xc1-B)/(hf*vc*vc))));
                Bc  = Bc *(1.0f-expf(-2.0f*fmaxf(0.0f,(Xc1-B)*(Xc -B)/(hf*vc*vc))));
            }

            Af = Af - 0.5f*hf*Xf;
            Ac = Ac - 0.5f*hc*Xc;
        }

        if (option==1) {
            Pf  = fmaxf(0.0f,Xf-K);
            Pc  = fmaxf(0.0f,Xc-K);
        }
        else if (option==2) {
            Pf  = fmaxf(0.0f,Af-K);
            Pc  = fmaxf(0.0f,Ac-K);
        }
        else if (option==3) {
            Pf  = Xf - Mf;
            Pc  = Xc - Mc;
        }
        else if (option==4) {
            Pf  = K*ncff((Xf0+r*Xf0*hf-K)/(sig*Xf0*sqrt(hf)));
            if (l==0) {Pc = Pf;}
            else {Pc = K * ncff((Xc0+r*Xc0*hc+sig*Xc0*dWf[0]-K)/(sig*Xc0*sqrt(hf)));}
        }
        else if (option==5) {
            Pf  = Bf*fmaxf(0.0f,Xf-K);
            Pc  = Bc*fmaxf(0.0f,Xc-K);
        }

        dP  = exp(-r*T)*(Pf-Pc);
        Pf  = exp(-r*T)*Pf;

        if (l==0) dP = Pf;

        sums[0] += nf;     // add number of timesteps as cost
        sums[1] += dP;
        sums[2] += dP*dP;
        sums[3] += dP*dP*dP;
        sums[4] += dP*dP*dP*dP;
        sums[5] += Pf;
        sums[6] += Pf*Pf;
    }
}


