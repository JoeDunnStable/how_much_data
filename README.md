# Observations on Taleb's "How much data do you need"

## Introduction

Nassim Nicholas Taleb has recently released a paper entitled *How much data do you need?
An operational metric for fat-tailedness*. to appear in the International Journal of Forecasting.  The paper focuses on the preasymptotic behavior of the sum of independent identically distributed random variables that are governed by various fat-tailed distributions.  He introduces
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

This package contains two algorithms implemented in C++, one of which relies on 
monte carlo simulations, and a second one that takes advantage of an integral representation of the MAD given the derivative of characteristic function of the underlying distribution.  I experimented with a third algorithm that used a discrete fourier transform, but it's been
dropped from the package because it's much inferior to the integral representation.

The package includes several files developed for this project:

*  monte_carlo_test.cpp.  The main program to run the monte carlo test.  It's bullet-proof
but so slow that it's hard to obtain the required accuracy in a reasonable time frame 
even using multi-threading.
* pinelis_taleb_test.cpp.   The main program implementing the integral representation
of MAD.  This is my current method of choice.
* pinelis_taleb_graph.R. An R script to graph the results of a run of pinelis_taleb_test.
*  pareto_distribution.h.  Contains a class for the pareto distribution, modeled on the 
classes in boost::random and including items normally computed in
boost::math::statistical_distributions
*  student_t_distribution.h. Similar to pareto but starts as a derived class from
boost::random::student_t_distribution
* exponential_distribution.h. A derived class from boost::random::exponential_distribution
* lognormal_distribution.h.  A derived class from boost::random::lognormal_distribution.  The
class implementing the distribution includes several versions of the calculation of the 
characteristic function: some based on numerical integration from the definition, one based on a asymptotic series for small angular frequency and one based on a approximation using Lambert W functions.
* normal_switch_mean.h. A class implementing a 50/50 mixture of two normal distributions both
with sigma=1, one of which has mean +d and the second of which has mean -d.
* normal_switch_stddev.h. A class implements a mixture of two normal distributions, the first of which has mean=0 and sigma=1 occuring with probability (1-p) and the second of which has mean=0 and sigma=a occuring with probability p.

In addition, i've included the apparatus I developed for adaptive integration, which seems to 
work much better that the Boost version that was originally used.  The code is spread 
through the following files:

* adaptive_integration.h.  The file that defines the interface to the adaptive integration routines.
* adaptive_integration_impl.h.  Most of the actual implementation.
* gauss_kronrad.h.  Interface to the functions used to calculate the nodes and weights 
for the integration.
* gauss_kronrod_impl.h.  The implementing code for gauss_kronrod.
* myfloat.h.  Some basic definitions allowing various types of multi-precision numbers.  Multiprecision is not used in the current implementation of how_much_data.

Finally included is a test program.

* lognormal_test.cpp.  A program to test the various versions of the calculation of the 
characteristic function. .

In order to improve portability, a meson.build file is included, which allows an easy port 
to other systems once the needed packages are installed.

### Update Feb. 15, 2019

I noticed that Taleb covered the same topic in N. N. Taleb,[Statistical Consequences of Fat Tails](https://www.academia.edu/37221402/THE_STATISTICAL_CONSEQUENCES_OF_FAT_TAILS_TECHNICAL_INCERTO_COLLECTION_?auto=download).  The procedure outlined in Section 6.5 of that monograph using results of
Pinelis is of general applicability and is incorporated into pinelis_taleb_test.

## Further Documentation

The code has been documented using the doxygen system.  The result can be accessed in
the doc subdirectory of the output directory.  

## Observations

So far my results are close to Taleb's.  As 
Taleb mentions in his paper such distributions require huge amounts of data to produce 
reasonable estimates of the MAD and this fact is mirrored in the number of monte carlo runs
or in the size the arrays used in the fast fourier transform.  I'm pushing the limit of
my computer's capability for the cases where alpha = 1.25 or alpha = 1.5.  The convolution 
is limited by the size of available memory and the monte carlo approach is limited by the
amount of time and the number of processors available.  The pinelis_taleb integral bypasses
most of these problems, but it requires (1) a very accurate calculation of the characteristic
function and its derivative, and (2) some ingenuity in choosing the right contour for integration.

I've experimented with other measures of scale, such as the 95% confidence interval spread,
for which the amount of computation needed is much more modest, but these results may not 
be relevant if MAD is the measure which best characterizes the uncertainty.

### Relationship with the Stable Distribution

Almost all of the distributions modeled have an associated stable distribution that is 
asymptotically equivalent to the modeled distribution.  The kappa for stable distributions
is always 2 - alpha.  However, in most cases the MAD for the associated stable distribution is significantly different from the MAD for the modelled distribution and the divergence of kappa
from 2 - alpha is associated with the rate of convergence of the MAD to the MAD of the 
associated stable distribution.

### Lognormal Distribution

As Taleb mentions in his paper, the lognormal distribution is a well behaved thin tail distriubtion
for small sigma, but becomes a fat-tailed distribution for sigmas much in excess of one.
The weakest part of the calculation here is associated with the large sigma runs, These have three widely divergent scales associated with them.
For instance, when mu=0 and sigma=5, the mode is about 1.4e-11,
the mean and MAD are about 2.7e+5 and the stardard deviation is about 7.2e+10.
The MAD for the normal distribution is sqrt(2/pi)  times it's standard deviation and therefore the the MAD normallized by n^.5 must go from 2.7e+5 to 5.7e+10 to reach it's asymptotic state.

## Acknowledgements

1. The package makes heavy use of the Boost C++ headers available at 
[boost](http://www.boost.org). 
2. The routines to calculate the Gauss/Kronrod nodes and weights use the Eigen headers.
These are available at [Eigen](http://www.eigen.tuxfamily.org).
3. One of the calculations of the characteristic function for the lognormal distribution uses 
Lambert W functions.  I've used the C++ code for the complex Lambert W function from 
[Istvan Mezo's web page](https://sites.google.com/site/istvanmezo81/others).
4. The routines in the adaptive_integration routine started out life
as machine C++ translations of Fortran routines in QUADPACK, which is part of 
SLATEC and therefore in the public domain (http://en.wikipedia.org/wiki/QUADPACK).
The routines were then heavily modified to take advantage of the C++ language.
5. One of the modifications made to QUADPACK is the addition of the ability to calculate
the nodes and weights for the Gauss Kronrod integration on the fly.  For this purpose
Dirk Laurie's method is used [kronrod.ps](http://dip.sun.ac.za/~laurie/papers/kronrod/kronrod.ps).
The routines used here are C++ translations of Dirk Laurie's MATLAB code, which
is included in Walter Gautschi's OPQ suite [OPQ](https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html).

## License
The code included here is covered by the MIT license.
