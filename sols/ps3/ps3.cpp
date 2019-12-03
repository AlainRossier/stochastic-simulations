/*

*/

#include <iostream>
#include "mlmc_test.cpp"
#include "poissinv.h"
#include <random>
#include <functional>

using namespace std;

int M = 4; // refinement cost factor 
int n0 = 4; // number of steps in layer 0

void ps3_l(int, int, double *);

// declare generator and output distributions

std::default_random_engine rng;
std::uniform_real_distribution<float> uniform(0.0f,1.0f);
auto next_uniform = std::bind(std::ref(uniform), std::ref(rng));

//
// main code
//

int main(int argc, char **argv) {
    
    int N0 = 100;   // initial samples on each level
    int Lmin = 2;   // minimum refinement level
    int Lmax = 15;  // maximum refinement level
    int N = 20000; // samples for convergence tests
    int L = 5;      // levels for convergence tests

    float Eps[6];
    float Eps2[] = {1, 2, 5, 10, 20, 50};
    memcpy(Eps,Eps2,sizeof(Eps2));
    
    char  filename[32];
    FILE *fp;

    printf("\n ---- Biochemical kinetics ---- \n");

    sprintf(filename, "biochemical_kinetics.txt");
    fp = fopen(filename,"w");

    rng.seed(1234);
    mlmc_test(ps3_l, M, N, L, N0, Eps, Lmin, Lmax, fp);

    fclose(fp);
}


// Coupling between level l-1 and level l
void ps3_l(int l, int N, double *sums) {

    // initalization
    float T, hf, hc;
    int nc, nf;
    int XGf, XMf, XPf, XDf, XGc, XMc, XPc, XDc;
    float df, dc;
    int lf[5], lc[5];
    int poiss_min, poiss_abs;
    int changeG[] = {0, 0, 0, 0, 0};
    int changeM[] = {1, 0, 0, -1, 0};
    int changeP[] = {0, 1, -2, 0, -1};
    int changeD[] = {0, 0, 1, 0, 0};
    


    int initG[] = {1, 0, 0, 0, 0};
    int initM[] = {0, 1, 0, 1, 0};
    int initP[] = {0, 0, 2, 0, 1};
    int initD[] = {0, 0, 0, 0, 0};
    int finalG[] = {1, 0, 0, 0, 0};
    int finalM[] = {1, 1, 0, 0, 0};
    int finalP[] = {0, 1, 0, 0, 0};
    int finalD[] = {0, 0, 1, 0, 0};

    
    float Pf, Pc, dP;
    
    // model parameters

    nf = n0 * (int)pow(M, l);
    nc = nf / M;
    T = 1;
    hf = T / ((float) nf);
    hc = T / ((float) nc);
    
    // simulate

    for (size_t k(0); k<7; k++) sums[k] = 0.0;

    for (size_t np(0); np<N; np++) {
	XGf = 1; XMf = 0; XPf = 0; XDf = 0;
	XGc = 1; XMc = 0; XPc = 0; XDc = 0;
	
	for (size_t n(0); n<nf; n++) {
	    lf[0] = 25*XGf; 
	    lf[1] = 1000*XMf; 
	    lf[2] = 0.001*XPf*(XPf-1); 
	    lf[3] = 0.1*XMf; 
	    lf[4] = XPf; 

	    if (n % M == 0) {
		lc[0] = 25*XGc;
		lc[1] = 1000*XMc;
		lc[2] = 0.001*XPc*(XPc-1);
		lc[3] = 0.1*XMc;
		lc[4] = XPc;
	    }

	    for (size_t k(0); k<5; k++) {

		poiss_abs = poissinv(next_uniform(), hf*abs(lf[k]-lc[k]));
		poiss_min = poissinv(next_uniform(), hf*fminf(lf[k], lc[k]));
		    
		if (lf[k] > lc[k]) {
		    df = poiss_abs + poiss_min;
		    dc = poiss_min;
		    
		} else {
		    dc = poiss_abs + poiss_min;
		    df = poiss_min;
		}
	    
		// Reaction on finest level
		XGf += df*changeG[k];
		XMf += df*changeM[k];
		XPf += df*changeP[k];
		XDf += df*changeD[k];

		// Reaction on coarest level
		XGc += dc*changeG[k];
		XMc += dc*changeM[k];
		XPc += dc*changeP[k];
		XDc += dc*changeD[k];


		/*
		
		// Reaction on finest level
		df = poiss1 + poiss2;
		
		if (initG[k] > 0) {
		    df = fminf(df, XGf/initG[k]);
		}
		if (initM[k] > 0) {
		    df = fminf(df, XMf/initM[k]);
		}
		if (initP[k] > 0) {
		    df = fminf(df, XPf/initP[k]);
		}
		if (initG[k] > 0) {
		    df = fminf(df, XGf/initG[k]);
		}
		    
		    XGf += df*(finalG[k] - initG[k]);
		    XMf += df*(finalM[k] - initM[k]);
		    XPf += df*(finalP[k] - initP[k]);
		    XDf += df*(finalD[k] - initD[k]);

		    

		// Reaction on coarsest level

		dc = poiss1 + poiss3;

		if (initG[k] > 0) {
		    dc = fminf(dc, XGc/initG[k]);
		}
		if (initM[k] > 0) {
		    dc = fminf(dc, XMc/initM[k]);
		}
		if (initP[k] > 0) {
		    dc = fminf(dc, XPc/initP[k]);
		}
		if (initG[k] > 0) {
		    dc = fminf(dc, XGc/initG[k]);
		}
		    
		    XGc += dc*(finalG[k] - initG[k]);
		    XMc += dc*(finalM[k] - initM[k]);
		    XPc += dc*(finalP[k] - initP[k]);
		    XDc += dc*(finalD[k] - initD[k]);

		*/
	      
	
    


	    }

	    XGf = fmaxf(XGf, 0);
	    XMf = fmaxf(XMf, 0);
	    XPf = fmaxf(XPf, 0);
	    XDf = fmaxf(XDf, 0);

	    XGc = fmaxf(XGc, 0);
	    XMc = fmaxf(XMc, 0);
	    XPc = fmaxf(XPc, 0);
	    XDc = fmaxf(XDc, 0);

	    

	    
	}

	Pf = (float)XDf;
	Pc = (float)XDc;
	dP = (float)(Pf - Pc);

	sums[0] += nf;     // add number of timesteps as cost
	sums[1] += dP;
	sums[2] += dP*dP;
	sums[3] += dP*dP*dP;
	sums[4] += dP*dP*dP*dP;
	sums[5] += Pf;
	sums[6] += Pf*Pf;

    }
}

