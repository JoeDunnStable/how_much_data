# Observations on Taleb's "How much data do you need"

## Introduction

Nassim Nicholas Taleb has recently released a paper entitled *How much data do you need?
An operational metric for fat-tailedness*. to appear in the International Journal of Forecasting.  The
paper focuses on the preasymptotic behavior of the sum of independent identically distributed
random variables that are governed by various fat-tailed distributions.  He introduces a 
metric labeled "kappa" that's related to the growth rate of the mean absolute deviation (MAD) 
as the number of summands increases.  Kappa is defined by the following formula.

   kappa(n0,n) = 2- (log(n)-log(n0))/(log(MAD(n) - MAD(n0))
   
Normally distributed variables have a kappa of zero.  Thus the extent to which kappa is in 
excess of zero measures the "fat-tailedness".

## Purpose

Taleb analyzes a number of distributions to estimate their kappa.  He gives explicit formulae
for kappa(1,2), which he was able to derive analytically, but he also includes tables of
kappa(1,30) and kappa(1,100) for two distributions, namely the Pareto distribution and 
Student's t distribution.  This package is my effort to replicate those tables and perhaps add
others for the other distributions.  It's a work in process, as there are some items that I can't 
closely replicate.

## What's included

I've experimented with two algorithms both implemented in C++, one of which relies on 
monte carlo simulations and the other of which relies on the use of the discrete fourier 
transform of the characteristic function of the convolution of the distribution functions. 
The package includes four files with code.

*  convolution_test.cpp.  The main program to nun the convolution test.
*  monte_carlo_test.cpp.  The main program to run the monte carlo test.
*  pareto_distribution.h.  Contains a class for the pareto distribution, modeled on the classes in boost::random and
including items normally computed in boost::math::statistical_distributions
*  student_t_distribution.h. Similar to pareto but starts as a derived class from
boost::random::student_t_distribution

In order to improve portability, a meson.build file is included, which allows an easy port 
to other systems once the needed packages are installed.

## Observations

So far my results are close to Taleb's except for the cases where alpha is close to one.  As 
Taleb mentions in his paper such distributions require huge amounts of data to produce 
reasonable estimates of the MAD and this fact is mirrored in the number of monte carlo runs
or in the size the arrays used in the fast fourier transform.  I'm pushing the limit of
my computer's capability for the cases where alpha = 1.25 or alpha = 1.5.  The convolution 
is limited by the size of available memory and the monte carlo approach is limited by the
amount of time and the number of processors available.

I've experimented with other measures of scale, such as the 95% confidence interval spread,
for which the amount of computation needed is much more modest, but these results may not 
be relevant if MAD is the measure which best characterizes the uncertainty.

## To Do

* Extend to other distributions
* Add better documention using doxygen.
* Switch to the parallel version of the fftw.

## Acknowledgements

1. The package makes heavy use of the Boost C++ headers available at 
[boost](http://www.boost.org)).
2. The package uses the Eigen headers solely for the purpose of wrapping the fast fourier
transform code.  These are available at [Eigen](http://www.eigen.tuxfamily.org).
3. As distributed the Eigen header wraps the fftw3 available at [fftw]{http://fftw.org}.  The fftw uses
a GPL, so if you want a less restrictive license you should delete the line
#define EIGEN_FFTW_DEFAULT at the beginning of convolution, which well revert to the
kissfft.

## License
The code included here is covered the the MIT license.
