# Stochastic Simulation
Solutions in C++ of the problem sheets of the Stochastic Simulation lecture by Prof. Mike Giles for the CDT of Mathematics of Random Systems, Michaelmas 2019.

## Installation and set-up

### Requirements
You will need 
- C++ compiler with C++14 flag (I use the GNU GCC Compiler built inside the IDE Code::Blocks),
- the C++ library Boost, for general purpose computations such as matrix multiplication,
- gnuplot, to plot the figure using the commands below.

### Execute the program
The exectuable file is in `./bin/Debug/Problem_sheets`  
Execute `./Problem_sheets` to write the solutions to the problem sheets.  

### Solutions and plots
You can find the solutions in text files under `./sols/`, together with the plots.


## Problem sheet 1

### Plotting
- Exercise 2.c.  
`plot "ps_1_2c_empirical_cdf.data" using 1:2 title "Empirical CDF", "ps_1_2c_empirical_cdf.data" using 1:3 with lines title "Theoretical CDF"`
- Exercise 2.d.  
`set logscale x`  
`plot "ps_1_2d_confidence_intervals.data" using 1:2:3:4 with yerrorlines title "Confidence intervals", "ps_1_2d_confidence_intervals.data" using 1:5 title "True mean"`
- Exercise 3.c.  
`plot "ps_1_3c_empirical_cdf.data" using 1:2 title "Empirical CDF", "ps_1_3c_empirical_cdf.data" using 1:3 with lines title "Theoretical CDF"`
- Exercise 3.d.  
`set logscale x`  
`plot "ps_1_3d_confidence_intervals.data" using 1:2:3:4 with yerrorlines title "Confidence intervals for European call option value", "ps_1_3d_confidence_intervals.data" using 1:5 title "True mean"`
- Exercise 4.a.  
`set logscale x`  
`plot "ps_1_4a_antithetic_variance.data" using 1:2 with lines title "Variance estimate for European call option value", "ps_1_4a_antithetic_variance.data" using 1:3 with lines title "With antithetic variables"`  
- Exercise 4.b.  
`set logscale x`  
`plot "ps_1_4b_control_variate.data" using 1:2 with lines title "Variance estimate for European call option value", "ps_1_4b_control_variate.data" using 1:3 with lines title "With control variate"`  
- Exercise 5.a.  
`set logscale x`  
`plot "ps_1_5a_digital_put_value.data" using 1:2:3:4 with yerrorlines title "Confidence intervals for digital put option value"`
- Exercise 5.b.  
`set logscale x`  
`plot "ps_1_5b_digital_put_value_importance_sampling.data" using 1:2:3:4 with yerrorlines title "Confidence intervals for digital put option value with importance sampling"`
- Exercise 6.a.  
`set logscale x`  
`plot "ps_1_6a_bumping_delta_european_call.data" using 1:2:3:4 with yerrorlines title "Confidence intervals for delta with the bumping method, with the same random numbers", "ps_1_6a_bumping_delta_european_call.data" using 1:8 title "True delta"`  
`plot "ps_1_6a_bumping_delta_european_call.data" using 1:5:6:7 with yerrorlines title "Confidence intervals for delta with bumping method, with different random numbers", "ps_1_6a_bumping_delta_european_call.data" using 1:8 title "True delta"`  
`plot "ps_1_6a_bumping_vega_european_call.data" using 1:2:3:4 with yerrorlines title "Confidence intervals for vega with the bumping method, with the same random numbers", "ps_1_6a_bumping_vega_european_call.data" using 1:8 title "True vega"`  
`plot "ps_1_6a_bumping_vega_european_call.data" using 1:5:6:7 with yerrorlines title "Confidence intervals for vega with bumping method, with different random numbers", "ps_1_6a_bumping_vega_european_call.data" using 1:8 title "True vega"`  
- Exercise 6.b.  
`plot "ps_1_6b_ipa_delta_european_call.data" using 1:2:3:4 with yerrorlines title "Confidence intervals for delta with the IPA method", "ps_1_6b_ipa_delta_european_call.data" using 1:5 title "True delta"`  
`plot "ps_1_6b_ipa_vega_european_call.data" using 1:2:3:4 with yerrorlines title "Confidence intervals for vega with the IPA method", "ps_1_6b_ipa_vega_european_call.data" using 1:5 title "True vega"`  


## Problem sheet 2

### Plotting
- Exercise 2.a.  
`plot "ps_2_2a_weak_convergence_sde.data" using 1:2:3 with yerrorlines title "Confidence intervals for weak convergence method of GBM SDE"`  
- Exercise 2.b.  
``


