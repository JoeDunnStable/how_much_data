# Observations on Taleb's "How much data do you need"

## Introduction

Nassim Nicholas Taleb has recently released a paper entitled *How much data do you need?
An operational metric for fat-tailedness*. to appear in the International Journal of Forecasting.  The paper focuses on the preasymptotic behavior of the sum of independent identically distributedrandom variables that are governed by various fat-tailed distributions.  He introduces
a metric labeled "kappa" that's related to the growth rate of the mean absolute deviation (MAD) 
as the number of summands increases.  Kappa is defined by the following formula.

   kappa(n0,n) = 2- (log(n)-log(n0))/(log(MAD(n) - MAD(n0))
   
Normally distributed variables have a kappa of zero.  Thus the extent to which kappa is in 
excess of zero measures the "fat-tailedness".

## Purpose

Taleb analyzes a number of distributions to estimate their kappa.  He gives explicit formulae
for kappa(1,2), which he was able to derive analytically, but he also includes tables of
kappa(1,30) and kappa(1,100) for two distributions, namely the Pareto distribution and 
Student's t distribution.  This package is my effort to replicate those tables and perhaps add
others for the other distributions.  It's a work in progress, as there are some items that I can't 
closely replicate.  However, none of the differences are large enough to call into question
Taleb's results.

## What's included

I've experimented with two algorithms both implemented in C++, one of which relies on 
monte carlo simulations and the other of which relies on the use of the discrete fourier 
transform of the characteristic function of the convolution of the distribution functions. 
The package includes four files with code.

*  convolution_test.cpp.  The main program to nun the convolution test.
*  monte_carlo_test.cpp.  The main program to run the monte carlo test.
*  pareto_distribution.h.  Contains a class for the pareto distribution, modeled on the 
classes in boost::random and including items normally computed in
boost::math::statistical_distributions
*  student_t_distribution.h. Similar to pareto but starts as a derived class from
boost::random::student_t_distribution
* exponential_distribution.h. A derived class from boost::random::exponential_distribution
* lognormal_distribution.h.  A derived class from boost::random::lognormal_distribution.  The
class implementing the distribution includes several versions of the calculation of the 
characteristic function, some based on numerical integration from the definition, one based 
on a p-spline approximation to the more accurate but much slower integrals, and one based
on a approximation using Lambert W functions.
* lognormal_test.cpp.  A program to test the various versions of the calculation of the 
characteristic function.  So far the Lambert W version seems to be the best compromise
between speed and accuracy, but it's not without problems.

In order to improve portability, a meson.build file is included, which allows an easy port 
to other systems once the needed packages are installed.

## Further Documentation

The code has been documented using the doxygen system.  The result can be accessed in
the doc subdirectory of the output directory.  The html version is accessed through 
html/index.html

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

### Relationship with the Statble Distribution

Almost all of the distributions modeled have an associated stable distribution that is 
asymptotically equivalent to the modeled distribution.  The kappa for stable distributions
is always 2 - alpha.  However, in most cases the MAD for the associated stable distribution is significantly different from the MAD for the modelled distribution and the divergence of kappa
from 2 - alpha is associated with the rate of convergence of the MAD to the MAD of the 
associated stable distribution.

### Lognormal Distribution

As Taleb mentions in his paper, the lognormal distribution is a well behaved thin tail distriubtion
for small sigma, but becomes a fat-tailed distribution for sigmas much in excess of one.
The weakest part of the calculation here is associated with the large sigma runs of the
for the convolution test of the lognormal distribution.  These have three widely divergent scales
associated with them.
For instance, when mu=0 and sigma=5, the mode is about 1.4e-11,
the mean and MAD are about 2.7e+5 and the stardard deviation is about 7.2e+10.
The MAD for the normal distribution is sqrt(2/pi)  times it's standard deviation and therefore the the MAD normallized by n^.5 must go from 2.7e+5 to 5.7e+10 to reach it's asymptotic state.

## To Do

* Tune the number of threads used by the parallel version of fftw.  The number of threads is
currently set at 8, the number of processors on my machine, which produces no net improvement in the overall runtime vs the baseline of no multithreading.

## Acknowledgements

1. The package makes heavy use of the Boost C++ headers available at 
[boost](http://www.boost.org).  The version must be at least 1.68 in order to compute
the complex integrals.
2. The package uses the Eigen headers for the purpose of wrapping the fast fourier
transform code and holding the large arrays that are used by the fft.  It also uses unsupported Eigen headers for the calculation of splines.
These are available at [Eigen](http://www.eigen.tuxfamily.org).
3. As distributed the Eigen header wraps the fftw3 available at [fftw](http://fftw.org).  The fftw uses
a GPL, so if you want a less restrictive license you should delete the line
#define EIGEN_FFTW_DEFAULT at the beginning of convolution_test.cpp, which will revert to the
kissfft.  Warning: kissfft uses much more memory and causes a severe slow down on my
machine for the larger arrays because of swapping activiity.
4. One of the calculations of the characteristic function for the lognormal distribution uses 
Lambert W functions.  I've used the C++ code for the complex Lambert W function from 
[Istvan Mezo's web page](https://sites.google.com/site/istvanmezo81/others).

## License
The code included here is covered the the MIT license.
