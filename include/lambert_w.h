//
/// \file lambert_w.h
/// source: https://sites.google.com/site/istvanmezo81/others
/// \package how_much_data
//
/// \author Created by Joseph Dunn on 2/20/19.
/// source: https://sites.google.com/site/istvanmezo81/others
/// \copyright © 2019 Joseph Dunn. All rights reserved.
//

#ifndef lambert_w_h
#define lambert_w_h

#include <complex>
using std::complex;
#include <boost/math/constants/constants.hpp>

// functions used for Lambert's W

/// z * exp(z)
template<typename RealType>
complex<RealType> zexpz(complex<RealType> z)
{
  return z * exp(z);
}

///The derivative of z * exp(z) = exp(z) + z * exp(z)
template<typename RealType>
complex<RealType> zexpz_d(complex<RealType> z)
{
  return exp(z) + z * exp(z);
}

///The second derivative of z * exp(z) = 2. * exp(z) + z * exp(z)
template<typename RealType>
complex<RealType> zexpz_dd(complex<RealType> z)
{
  return 2. * exp(z) + z * exp(z);
}

///Determine the initial point for the root finding
template<typename RealType>
complex<RealType> InitPoint(complex<RealType> z, int k)
{
  const RealType pi = boost::math::constants::pi<RealType>();
  const RealType e{ boost::math::constants::e<RealType>()};
  complex<RealType> I{0, 1};
  complex<RealType> two_pi_k_I{ 0., 2. * pi * k };
  complex<RealType> ip{ log(z) + two_pi_k_I - log(log(z) + two_pi_k_I) };// initial point coming from the general asymptotic approximation
  complex<RealType> p{ sqrt( 2. * (e * z + 1.) ) };// used when we are close to the branch cut around zero and when k=0,-1
  
  if (abs(z - (-exp(-1.))) <= 1.) //we are close to the branch cut, the initial point must be chosen carefully
  {
    if (k == 0) ip = -1. + p - 1./3. * pow(p, 2) + 11./72. * pow(p, 3);
    if (k == 1 && z.imag() < 0.) ip = -1. - p - 1./3. * pow(p, 2) - 11. / 72. * pow(p, 3);
    if (k == -1 && z.imag() > 0.) ip = -1. - p - 1./3. * pow(p, 2) - 11./72. * pow(p, 3);
  }
  
  if (k == 0 && abs(z - .5) <= .5) ip = (0.35173371 * (0.1237166 + 7.061302897 * z)) / (2. + 0.827184 * (1. + 2. * z));// (1,1) Pade approximant for W(0,a)
  
  if (k == -1 && abs(z - .5) <= .5) ip = -(((2.2591588985 +
                                             4.22096*I) * ((-14.073271 - 33.767687754*I) * z - (12.7127 -
                                                                                                19.071643*I) * (1. + 2.*z))) / (2. - (17.23103 - 10.629721*I) * (1. + 2.*z)));// (1,1) Pade approximant for W(-1,a)
  
  return ip;
}

/// Labert's W function for complex arguments.
template<typename RealType>
complex<RealType> LambertW(complex<RealType> z, int k = 0)
{
  //For some particular z and k W(z,k) has simple value:
  if (z == 0.) return (k == 0) ? 0. : -INFINITY;
  if (z == -exp(-1.) && (k == 0 || k == -1)) return -1.;
  if (z == exp(1.) && k == 0) return 1.;
  
  //Halley method begins
  complex<RealType> w{ InitPoint(z, k) }, wprev{ InitPoint(z, k) }; // intermediate values in the Halley method
  const unsigned int maxiter = 30; // max number of iterations. This eliminates improbable infinite loops
  unsigned int iter = 0; // iteration counter
  RealType prec = 1.E-30; // difference threshold between the last two iteration results (or the iter number of iterations is taken)
  
  do
  {
    wprev = w;
    w -= 2.*((zexpz(w) - z) * zexpz_d(w)) /
    (2.*pow(zexpz_d(w),2) - (zexpz(w) - z)*zexpz_dd(w));
    iter++;
  } while ((abs(w - wprev) > prec) && iter < maxiter);
  return w;
}

/// Derivative Labert's W function for complex arguments.
template<typename RealType>
complex<RealType> LambertWPrime(complex<RealType> z, int k = 0)
{
  return complex<RealType>(1)/(z + exp(LambertW(z,k)));
}



#endif /* lambert_w_h */
