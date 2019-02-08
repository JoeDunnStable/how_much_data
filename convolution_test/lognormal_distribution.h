//
/// \file lognormal_distribution.h
/// \package  how_much_data
//
/// \author Created by Joseph Dunn on 1/3/19.
/// \copyright © 2019 Joseph Dunn. All rights reserved.
//

#ifndef lognormal_distribution_h
#define lognormal_distribution_h

#include <iostream>
using std::cout;
using std::endl;
using std::ostream;
#include <iomanip>
using std::setw;
using std::setprecision;
using std::scientific;
using std::defaultfloat;

#include <complex>
using std::complex;
#include <string>
using std::string;
#include <memory>
using std::shared_ptr;
using std::make_shared;
#include <boost/random.hpp>
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
using boost::math::quadrature::gauss_kronrod;
#include <boost/math/quadrature/exp_sinh.hpp>
using boost::math::quadrature::exp_sinh;
#include "p_spline.h"

/** a class with functions related to the lognormal distribution.
 Describes a random variable distributed as exp(sigma X + mu)
 where X is normally distributed
*/
template<class RealType = double>
struct lognormal_distribution : public boost::random::lognormal_distribution<RealType> {
  /// the constructor for the distribution
  lognormal_distribution(RealType mu,        ///< [in] the mu parameter
                         RealType sigma,      ///< [in] the sigma parameter
                         int type = 1        ///< [in] type of cf calculaiton
                                             ///< 1 ln(x), 2 x, 3 w, 4 spline
                         ) : boost::random::lognormal_distribution<RealType>(mu, sigma),
                             ln_dist(mu, sigma), _raw_mean(exp(mu+pow(sigma,2)/2)),
                             _raw_stddev(sqrt(exp(pow(sigma,2))-1)*exp(mu+pow(sigma,2)/2)),
                             _cf_spline(512, 300, 128, this, type), _type(type) {}
  
  /// return a random number from the normalized distribution
  template <typename Engine>
  RealType operator() (Engine& eng)
  {
    return (boost::random::lognormal_distribution<RealType>::operator()(eng)-_raw_mean)/_raw_stddev;
    
  }
  
  /// return the cumulative distribution function
  RealType cdf(RealType x,                    ///< [in] the quantile variable
               bool lower_tail = true         ///< [in] flag indicating which tail to use
               ) const
  {
    RealType xx = _raw_stddev*x+_raw_mean;
    if (xx < 0)
      return lower_tail ? 0 : 1;
    else
      return lower_tail ? boost::math::cdf(ln_dist, xx)
                        : boost::math::cdf(boost::math::complement(ln_dist, xx));
  }
  
  /// return the probability density function
  RealType pdf(RealType x                      ///< the quantile variable
               ) const
  {
    RealType xx = _raw_stddev * x + _raw_mean;
    return _raw_stddev * boost::math::pdf(ln_dist, xx);
  }
  
  /// return the quantile corresponding to a given propability
  RealType quantile(RealType p                 ///< the target probability
                   ) const {
    return (boost::math::quantile(ln_dist, p)-_raw_mean)/_raw_stddev;
  }
  
  /// return the mu parameter of the distribution
  RealType mu() const { return boost::random::lognormal_distribution<>::m();}
  
  /// return sigma parameter of the distribution
  RealType sigma() const { return boost::random::lognormal_distribution<>::s();}
  
  /// return the alpha parameter of the asymptotically equivalent stable distribution
  RealType alpha_stable() const { return static_cast<RealType>(2);}
  
  /** Returns the smallest value that the distribution can produce. */
  RealType min () const
  { return -_raw_mean/_raw_stddev; }
  /** Returns the largest value that the distribution can produce. */
  RealType max () const
  { return (std::numeric_limits<RealType>::infinity)(); }
  
  /** Return mean of distribution */
  RealType mean() const {return RealType(0);}
  
  /** Return mean absolute deviation of the distribution */
  RealType mad() const {
    return (2*exp(mu() + pow(sigma(),2)/2.)*boost::math::erf(sigma()/(2.*sqrt(2))))/_raw_stddev;
  }
  
  /** Return the mad of the square of the distribution.
   this is the approximation used by Taleb
  */
  RealType mad2() const {
    return (4*exp(mu() + pow(sigma(),2)/2.)
           *boost::math::erf(sqrt(log(.5*(exp(pow(sigma(),2))+1)))/(2*sqrt(2))))/_raw_stddev;
  }
  
  /** Return the confidence interval of the distribution */
  RealType ci(RealType level = RealType(.05)) const
  {
    return (boost::math::quantile(ln_dist,1-level/2)-boost::math::quantile(ln_dist,level/2))/_raw_stddev;
  }
  
  /// return the characteristic function via Fourier integral in x domain
  complex<RealType> cf_fourier_x (RealType omega  ///< the angular frequency
                                           ) const
  {
    omega *= exp(mu())/_raw_stddev;
    complex<RealType> ret;
    if (omega == 0) {
      ret= static_cast<complex<RealType> >(1);
    }
    else {
      // First deal with the near Dirac delta near the origin
      RealType lo_x = 1e-9;
      ret = boost::math::cdf(ln_dist, lo_x);
      // cap the upper range when pdf is samll or frequency is high enough
      // to assure cancelation.
      const RealType eps = std::numeric_limits<RealType>::epsilon();
      const RealType tol = boost::math::tools::root_epsilon<RealType>();
      RealType hi_x = boost::math::quantile(boost::math::complement(ln_dist,eps));
      hi_x = std::min(hi_x, 1/(fabs(omega)*tol));
      // the integration range will be mapped to y in [0,1]
      RealType lo_y=sqrt(lo_x)/(1+sqrt(lo_x));
      RealType hi_y=sqrt(hi_x)/(1+sqrt(hi_x));
      CFIntegrandX cf_integrand( omega, sigma());
      const int kronrod_order = 15;
      unsigned max_depth = 16;
      ret += (hi_y > lo_y) ? gauss_kronrod<RealType, kronrod_order>::integrate(cf_integrand, lo_y, hi_y,
                                                   max_depth, tol)
                       :0;
    }
    return exp(complex<RealType>{0, -_raw_mean * omega}) * ret;
  }
  
  /// return the approximate characteristic function using Lambert W funciton
  /// Much faster than either integral but has problems when sigma > 1
  complex<RealType> cf_lambert_w (RealType omega  ///< the angular frequency
  ) const
  {
    omega *= exp(mu())/_raw_stddev;
    complex<RealType> ret;
    if (omega == 0) {
      ret= static_cast<complex<RealType> >(1);
    }
    else {
      const complex<RealType> one{1};
      const complex<RealType> two(2);
      complex<RealType> w = LambertW(complex<RealType>(0,-omega * pow(sigma(),2)));
      ret = exp(-(pow(w,2) + two*w)/(two * pow(sigma(),2)))/sqrt(one + w);
    }
    
    return exp(complex<RealType>{0, -_raw_mean* omega}) * ret;
  }
  
  /// Return the characteristic function via Fourier integral in ln(x) domain
  complex<RealType> cf_fourier_lnx(RealType omega  ///< the angular frequency
                       ) const
  {
    omega *= exp(mu())/_raw_stddev;
    complex<RealType> ret;
    if (omega == 0) {
      ret= static_cast<complex<RealType> >(1);
    }
    else {
      // Low frequency, oscillation is not a problem
      boost::math::normal_distribution<RealType> nd(0,sigma());
      CFIntegrandLnX cf_integrand( omega, nd);
      // exclude low end with insignificant pdf
      const RealType eps = std::numeric_limits<RealType>::epsilon();
      RealType lo = boost::math::quantile(nd, eps);
      // exclude high end with insifnificant pdf or high cancelation
      const RealType tol = boost::math::tools::root_epsilon<RealType>();
      RealType hi = boost::math::quantile(boost::math::complement(nd,eps));
      hi = std::min(hi, -log(fabs(omega))-log(tol));
      const int kronrod_order = 15;
      const unsigned max_depth = 16;
      
      ret = (hi > lo) ?gauss_kronrod<RealType,kronrod_order>::integrate(cf_integrand, lo, hi,
                                                             max_depth, tol)
                      :0;
    }
    return exp(complex<RealType>{0, -_raw_mean * omega}) * ret;
  }
  
  /// return the approximate characteristic function using precompluted
  /// cubic splines
  complex<RealType> characteristic_function (RealType omega  ///< [in] the angular frequency
                                            ) const
  {
    switch ( _type ) {
      case 1:
        return cf_fourier_lnx(omega);
      case 2:
        return cf_fourier_x(omega);
      case 3:
        return cf_lambert_w(omega);
      case 4:
        return _cf_spline(omega);
      default:
        return cf_fourier_lnx(omega);
    }
  }
  
  /**
   * Write distribution to std::ostream
   */
  template <class charT, class traits>
  friend std::basic_ostream<charT,traits>& operator<< ( std::basic_ostream<charT,traits>& os,
                                                       const lognormal_distribution& dist ) {
    os << "lognormal distriubution w mu = " << dist.mu()
       << ", & sigma = " << dist.sigma() << endl;
    os << "Mean   = " << scientific << setw(12)<< setprecision(5)<< dist.mean() << endl
    << "min    = " << scientific << setw(12)<< setprecision(5)<< dist.min() << endl
    << "MAD    = " << scientific << setw(12)<< setprecision(5)<< dist.mad() << endl
    << "MAD2   = " << scientific << setw(12)<< setprecision(5)<< dist.mad2() << endl
    << "kappa2 = " << defaultfloat << setw(12) << 2 - log(2)/(log(dist.mad2())-log(dist.mad())) << endl
    << "95% CI = " << scientific << setw(12)<< setprecision(5)<< dist.ci()
    << defaultfloat << endl;
    return os;
  }
private:
  
  // functions used for Lambert's W

  /// z * exp(z)
  static complex<RealType> zexpz(complex<RealType> z)
  {
    return z * exp(z);
  }
  
  ///The derivative of z * exp(z) = exp(z) + z * exp(z)
  static complex<RealType> zexpz_d(complex<RealType> z)
  {
    return exp(z) + z * exp(z);
  }
  
  ///The second derivative of z * exp(z) = 2. * exp(z) + z * exp(z)
  static complex<RealType> zexpz_dd(complex<RealType> z)
  {
    return 2. * exp(z) + z * exp(z);
  }
  
  ///Determine the initial point for the root finding
  static complex<RealType> InitPoint(complex<RealType> z, int k)
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
  /// Source: https://sites.google.com/site/istvanmezo81/others
  static complex<RealType> LambertW(complex<RealType> z, int k = 0)
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
  
  /// Compute the pdf for complex argurments.  Used for characterietic function
  /// Note mu is assumed zero in this routine
  static RealType pdf_real(RealType x, RealType sigma) {
    return RealType(1)/(exp(pow(log(x),2)/(2.*pow(sigma,2)))*
                        sqrt(2*pi<RealType>())*x*sigma);
  }
  
  /// class of integrand used for Fourier integral in x domain
  class CFIntegrandX {
  public:
    CFIntegrandX(RealType omega, RealType sigma) : omega(omega), sigma(sigma){}
    complex<RealType> operator() (RealType y) const {
      // map the interval [0,1] into [0,infinity
      RealType x = pow(y/(1-y),2);
      RealType dx = 2 * y/pow((1-y),3);
      return dx* exp(complex<RealType>(0,x * omega))
                * pdf_real(x,sigma);
    }
  private:
    RealType omega;
    RealType sigma;
  };

  /// class of integrand used for Foureir integral in ln(x) domain
  class CFIntegrandLnX {
  public:
    CFIntegrandLnX(RealType omega,
                const boost::math::normal_distribution<RealType> nd
                ) : omega(omega), nd(nd) {}
    complex<RealType> operator() (RealType ln_x) const {
      RealType omega_exp_x = omega * exp(ln_x);
      if (!isfinite(omega_exp_x))
        return 0;
      else
        return exp(complex<RealType>(0,omega_exp_x))* boost::math::pdf(nd,ln_x);
    }
  private:
    RealType omega;
    boost::math::normal_distribution<RealType> nd;
  };

  /// class using cubic splines to evalute characteristic function
  class CFSpline {
    typedef Spline<RealType, 2, 3> Spl;
    
  public:
    CFSpline (int m,                  ///< number of samples between 0 & omegamax
              RealType omegamax,      ///< the maximum omega
              int n,                  ///< the number of knots between o & omegamax
              const lognormal_distribution *ln,            ///< the distribution
              int type
              ) : _omegamax(omegamax)
    {
      if (type != 4) return;
      Matrix<double, 2, Dynamic> y(2,2*m+1);
      Matrix<double, Dynamic, 1> x(2*m+1);
      x.setLinSpaced(2*m+1, -omegamax, omegamax);
      y.col(m) << 1, 0;
      for (int i=1; i<=m; ++i) {
        complex<double> cf = ln->cf_fourier_lnx(x(m+i));
        y.col(m+i) << real(cf), imag(cf);
        y.col(m-i) << real(cf), -imag(cf);
      }
      x.array() += omegamax;
      Matrix<RealType, Dynamic, 1> knots(2*n+1);
      knots.setLinSpaced(2*n+1,0,2*omegamax);
      double penality = .000000;
      _spline = SplineFitPSpline<Spl, Matrix<RealType, 2, Dynamic> >(y, x, 3, knots, penality);
    }
    
    complex<RealType> operator() (RealType omega) const {
      Matrix<RealType, 2, 1> res = _spline(omega+_omegamax);
      return complex<RealType> {res(0), res(1)};
    }
    
  private:
    RealType _omegamax;
    Spl _spline;
  };
    
  RealType _raw_mean;
  RealType _raw_stddev;
  CFSpline _cf_spline;
  boost::math::lognormal_distribution<RealType> ln_dist;
  int _type;

};


#endif /* lognormal_distribution_h */
