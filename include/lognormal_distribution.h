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
using std::right;
using std::scientific;
using std::defaultfloat;

#include <complex>
using std::complex;
#include <string>
using std::string;
using std::to_string;
#include <memory>
using std::shared_ptr;
using std::make_shared;
#include <utility>
using std::pair;
using std::make_pair;
#include <boost/random.hpp>
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;
using boost::math::constants::one_div_root_two_pi;
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/special_functions/erf.hpp>

using boost::math::tools::root_epsilon;

#include <boost/math/tools/roots.hpp>
using boost::math::tools::newton_raphson_iterate;

#include "adaptive_integration.h"
using namespace adaptive_integration;

#include "lambert_w.h"

template <typename RealType>
complex<RealType> expm1(complex<RealType> z) {
  if (abs(z) > .5)
    return exp(z) - complex<RealType>{1};
  else {
    RealType x = z.real();
    RealType y = z.imag();
    return complex<RealType>(expm1(x)*cos(y) - 2*pow(sin(y/2),2),
                             exp(x)*sin(y));
  }
}

template <typename RealType>
RealType sign(RealType x) {
  const RealType zero{0};
  return static_cast<RealType>((x>zero)-(x<zero));
}

/** a class with functions related to the lognormal distribution.
 Describes a random variable distributed as exp(sigma X + mu) adjusted
 to zero mean  where X is normally distributed.
*/
template<class RealType = double>
struct lognormal_distribution : public boost::random::lognormal_distribution<RealType> {
  /// the constructor for the distribution
  lognormal_distribution(RealType mu,        ///< [in] the mu parameter
                         RealType sigma,      ///< [in] the sigma parameter
                         IntegrationController<RealType>& cf_ctl, ///< [in} reference controller
                         int type = 1        ///< [in] type of cf calculaiton
                                             ///< 1 ln(x), 2 bespoke, 3 w, 4 mixed
                         ) : boost::random::lognormal_distribution<RealType>(mu, sigma),
                             ln_dist(mu, sigma), _type(type), _cf_ctl(cf_ctl),
                             _raw_mean(exp(mu+sigma*sigma/2))
  {
    RealType a = exp(sigma*sigma/2);
    central_moments.push_back(1);
    central_moments.push_back(0);
    central_moments.push_back(pow(a,2)*(pow(a,2)-1));
    central_moments.push_back(2*pow(a,3)-3*pow(a,5)+pow(a,9));
    central_moments.push_back(-3*pow(a,4)+6*pow(a,6)-4*pow(a,10)+pow(a,16));
    central_moments.push_back(4*pow(a,5)-10*pow(a,7)+10*pow(a,11)-5*pow(a,17)+pow(a,25));
                              
  }
  
  /// return a random number from the normalized distribution
  template <typename Engine>
  RealType operator() (Engine& eng)
  {
    return boost::random::lognormal_distribution<RealType>::operator()(eng)-_raw_mean;
    
  }
  
  /// return the cumulative distribution function
  RealType cdf(RealType x,                    ///< [in] the quantile variable
               bool lower_tail = true         ///< [in] flag indicating which tail to use
               ) const
  {
    if (x < min())
      return lower_tail ? 0 : 1;
    else
      return lower_tail ? boost::math::cdf(ln_dist, x + _raw_mean)
                        : boost::math::cdf(boost::math::complement(ln_dist, x + _raw_mean));
  }
  
  /// return the probability density function
  RealType pdf(RealType x                      ///< the quantile variable
               ) const
  {
    return boost::math::pdf(ln_dist, x+_raw_mean);
  }
  
  /// return the quantile corresponding to a given propability
  RealType quantile(RealType p                 ///< the target probability
                   ) const {
    return boost::math::quantile(ln_dist, p) - _raw_mean;
  }
  
  /// return the mu parameter of the distribution
  RealType mu() const { return boost::random::lognormal_distribution<>::m();}
  
  /// return sigma parameter of the distribution
  RealType sigma() const { return boost::random::lognormal_distribution<>::s();}
  
  /// return the alpha parameter of the asymptotically equivalent stable distribution
  RealType alpha_stable() const { return static_cast<RealType>(2);}
  
  /** Returns the smallest value that the distribution can produce. */
  RealType min () const
  { return -_raw_mean; }
  /** Returns the largest value that the distribution can produce. */
  RealType max () const
  { return (std::numeric_limits<RealType>::infinity)(); }
  
  /** Return mean of distribution */
  RealType mean() const {return 0;}
  
  /** Return mean absolute deviation of the distribution */
  RealType mad() const {
    return (2*exp(mu() + pow(sigma(),2)/2.)*boost::math::erf(sigma()/(2.*sqrt(2))));
  }
  
  /** Return the mad of the square of the distribution.
   this is the approximation used by Taleb
  */
  RealType mad2() const {
    return (4*exp(mu() + pow(sigma(),2)/2.)
           *boost::math::erf(sqrt(log(.5*(exp(pow(sigma(),2))+1)))/(2*sqrt(2))));
  }
  
  /** Return the confidence interval of the distribution */
  RealType ci(RealType level = RealType(.05)) const
  {
    return (boost::math::quantile(ln_dist,1-level/2)-boost::math::quantile(ln_dist,level/2));
  }
  
  /// cf via asymptotic series for small omega
  complex<RealType> cf_series(complex<RealType> omega,  ///< [in] the angular frequency
                                  RealType* _error = nullptr, ///< [out] the integration error
                                  RealType* _l1_norm = nullptr, ///< [out] the l1 norm
                                  int* _neval = nullptr          ///< [out] the # of evaluations
                                )
  {
    complex<RealType> i{0,1};
    complex<RealType> ret{1};
    complex<RealType> fac = i * omega;
    RealType old_size{1};
    for (int n = 2; n< int(central_moments.size()); n++) {
      fac *= i * omega/RealType(n);
      complex<RealType> term = central_moments.at(n) * fac;
      if (abs(term) > old_size) break;
      ret += term;
      old_size = abs(term);
    }
    if(_error) *_error = old_size;
    if(_l1_norm) *_l1_norm = std::numeric_limits<RealType>::quiet_NaN();
    if(_neval) *_neval = 0;
    
    return ret;
  }
  
  /// derivative of cf via asymptotic series for small omega
  complex<RealType> cfprime_series(complex<RealType> omega,  ///< [in] the angular frequency
                                  RealType* _error = nullptr, ///< [out] the integration error
                                  RealType* _l1_norm = nullptr, ///< [out] the l1 norm
                                  int* _neval = nullptr          ///< [out] the # of evaluations
  )
  {
    complex<RealType> i{0,1};
    complex<RealType> ret{0};
    complex<RealType> fac = i;
    RealType old_size{1};
    for (int n = 2; n< int(central_moments.size()); n++) {
      fac *= i * omega/RealType(n);
      complex<RealType> term = RealType(n) * central_moments.at(n) * fac;
      if (abs(term) > old_size) break;
      ret += term;
      old_size = abs(term);
    }
    if(_error) *_error = old_size;
    if(_l1_norm) *_l1_norm = std::numeric_limits<RealType>::quiet_NaN();
    if(_neval) *_neval = 0;
    
    return ret;
  }
  
  //  /// return the approximate characteristic function using Lambert W funciton
  /// Much faster than the fourier integrals but has problems when sigma > 1
  complex<RealType> cf_lambert_w (complex<RealType> omega,  ///< [in] the angular frequency
                                  RealType* _error = nullptr, ///< [out] the integration error
                                  RealType* _l1_norm = nullptr, ///< [out] the l1 norm
                                  int* _neval = nullptr       ///< [out] the # of evaluations
                                 )
  {
    omega *= exp(mu());
    complex<RealType> ret;
    const complex<RealType> i{0,1};
    if (_neval) *_neval = 0;
    if (abs(omega) == 0) {
      if (_error) *_error=std::numeric_limits<RealType>::quiet_NaN();
      if (_l1_norm) *_l1_norm=std::numeric_limits<RealType>::quiet_NaN();
      ret= static_cast<complex<RealType> >(1);
    }
    else {
      const complex<RealType> i{0.,1.};
      const complex<RealType> one{1};
      const complex<RealType> two(2);
      RealType ps_sq = expm1(pow(sigma(),2));
      complex<RealType> w = LambertW(- i * omega * ps_sq);
      ret = exp(-(pow(w,2) + two*w)/(two * ps_sq))/sqrt(one + w);
    }
    if(_error) *_error = std::numeric_limits<RealType>::quiet_NaN();
    if(_l1_norm) *_l1_norm = std::numeric_limits<RealType>::quiet_NaN();

    return exp(-i * _raw_mean * omega) * ret;
  }
  
  /// return the approximate derivative of char. function using Lambert W funciton
  /// Much faster than the fourier integrals but has problems when sigma > 1
  complex<RealType> cfprime_lambert_w (complex<RealType> omega,  ///< [in] the angular frequency
                                       RealType* _error = nullptr, ///< [out] the integration error
                                       RealType* _l1_norm = nullptr, ///< [out] the l1 norm
                                       int* _neval = nullptr       ///< [out] the # of evaluations
                                      )
  {
    complex<RealType> ret;
    const complex<RealType> i{0,1};
    if (_neval) *_neval = 0;
    if (abs(omega) == 0) {
      if (_error) *_error=std::numeric_limits<RealType>::quiet_NaN();
      if (_l1_norm) *_l1_norm=std::numeric_limits<RealType>::quiet_NaN();
      ret= complex<RealType>(0.,mean());
    }
    else {
      const complex<RealType> one{1};
      const complex<RealType> two(2);
      const complex<RealType> i{0,1};
      RealType ps_sq = 2*expm1(pow(sigma(),2)/2);
      complex<RealType> arg {-i * omega * exp(mu()) * ps_sq};
      complex<RealType> darg_domega = -i * exp(mu()) * ps_sq;
      complex<RealType> w = LambertW(arg);
      complex<RealType> dw_domega = LambertWPrime(arg)*darg_domega;
      ret = exp(-(pow(w,2) + two*w)/(two * ps_sq))
            * ((-w - one)/ps_sq/sqrt(one+w) -.5/pow(one+w,-1.5))
            * dw_domega;
    }
    if(_error) *_error = std::numeric_limits<RealType>::quiet_NaN();
    if(_l1_norm) *_l1_norm = std::numeric_limits<RealType>::quiet_NaN();
    return exp(-i * _raw_mean * omega) * ret
           - i * _raw_mean * cf_lambert_w(omega);
  }
  
  ///Return the characteristic function along a designer countour
  complex<RealType> cf_fourier(complex<RealType> omega,  ///< [in] the angular frequency
                               RealType* _error = nullptr, ///< [out] the integration error
                               RealType* _l1_norm = nullptr, ///< [out] the l1 norm
                               int* _neval = nullptr          ///< [out] the # of evaluations
  )
  {
    return cf<ContourZeroExponentImag> (omega, _error, _l1_norm, _neval);
  }
  
  /// Return the derivative of the char. fun. wrt omega using designer contour
  complex<RealType> cfprime_fourier(complex<RealType> omega,  ///< the angular frequency
                                    RealType* _error = nullptr, ///< [out] the integration error
                                    RealType* _l1_norm = nullptr, ///< [out] the
                                    int* _neval = nullptr     ///< the number of evaluations
  )
  {
    return cfprime<ContourZeroExponentImag>(omega, _error, _l1_norm, _neval);
  }
  
  /// Return the characteristic function via Fourier integral in ln(x) domain
  complex<RealType> cf_fourier_lnx(complex<RealType> omega,  ///< the angular frequency
                                   RealType* _error = nullptr, ///< [out] the integration error
                                   RealType* _l1_norm = nullptr, ///< [out] the
                                   int* _neval = nullptr       ///< [out] the # of evaluations
                       )
  {
    return cf<ContourZeroImag>(omega, _error, _l1_norm, _neval);
  }
  
  /// Return the derivative of char. function via Fourier integral in ln(x) domain
  complex<RealType> cfprime_fourier_lnx(complex<RealType> omega,  ///< [in] the angular frequency
                                        RealType* _error = nullptr, ///< [out] the integration error
                                        RealType* _l1_norm = nullptr, ///< [out] the l1 norm
                                        int* _neval = nullptr       ///< [out] the # of evaluations
                                       )
  {
    return cfprime<ContourZeroImag>(omega, _error, _l1_norm, _neval);
  }
  
  /// Return the char. fun. using mixed Lambert and shifted contour
  complex<RealType> cf_fourier_mixed(complex<RealType> omega,  ///< [in] the angular frequency
                                     RealType* _error = nullptr, ///< [out] the integration error
                                     RealType* _l1_norm = nullptr, ///< [out] the l1 norm
                                     int* _neval = nullptr       ///< [out] the # of evaluations
                                    )
  {
    if (sigma() < 1/10.*(pi<RealType>()/2-arg(omega)) ) {
      return cf_lambert_w(omega, _error, _l1_norm, _neval);
    } else {
      return cf<ContourPi2Imag>(omega, _error, _l1_norm, _neval);
    }
  }
  
  /// The derivative of the char. fun.
  complex<RealType> cfprime_fourier_mixed(complex<RealType> omega,  ///< [in] the angular frequency
                                     RealType* _error = nullptr, ///< [out] the integration error
                                     RealType* _l1_norm = nullptr, ///< [out] the l1 norm
                                     int* _neval = nullptr       ///< [out] the # of evaluations
  )
  {
    if (sigma() < 1/10.*(pi<RealType>()/2-arg(omega)) ) {
      return cfprime_lambert_w(omega, _error, _l1_norm, _neval);
    } else {
      return cfprime<ContourPi2Imag>(omega, _error, _l1_norm, _neval);
    }
  }
  
/// return the approximate characteristic function
  complex<RealType> characteristic_function (complex<RealType> omega,  ///< [in] the angular frequency
                                             RealType* _error = nullptr, ///< [out] the integration error
                                             RealType* _l1_norm = nullptr, ///< [out] the l1 norm
                                             int* _neval = nullptr      ///< {out} the # of evaluations
                                            )
  {
    if (real(omega) < 0.) return conj(characteristic_function(-conj(omega), _error, _l1_norm, _neval));
    switch ( _type ) {
      case 1:
        return cf_fourier_lnx(omega, _error, _l1_norm, _neval);
      case 2:
        return cf_fourier(omega, _error, _l1_norm, _neval);
      case 3:
        return cf_lambert_w(omega, _error, _l1_norm, _neval);
      case 4:
        return cf_fourier_mixed(omega, _error, _l1_norm, _neval);
      case 5:
        return cf_series(omega, _error, _l1_norm, _neval);
      default:
        return cf_fourier(omega, _error, _l1_norm, _neval);
    }
  }
  
  /// return the derivative of the  characteristic function
  complex<RealType> characteristic_function_prime (complex<RealType> omega,  ///< [in] the angular frequency
                                                   RealType* _error = nullptr, ///< [out] the integration error
                                                   RealType* _l1_norm = nullptr, ///< [out] the
                                                   int* _neval = nullptr       ///< [out] the # of evaluations
                                            )
  {
    if (real(omega) < 0.) return -conj(characteristic_function_prime(-conj(omega), _error, _l1_norm));
    switch ( _type ) {
      case 1:
        return cfprime_fourier_lnx(omega, _error, _l1_norm, _neval);
      case 2:
        return cfprime_fourier(omega, _error, _l1_norm, _neval);
      case 3:
        return cfprime_lambert_w(omega, _error, _l1_norm, _neval);
      case 4:
        return cfprime_fourier_mixed(omega, _error, _l1_norm, _neval);
      case 5:
        return cfprime_series(omega, _error, _l1_norm, _neval);
      default:
        return cfprime_fourier(omega, _error, _l1_norm, _neval);
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
  
  /// print integrand of CFIntegrand
  static void print_fourier_integrand(ostream& os, int n,
                              complex<RealType> omega, RealType sigma) {
    CFIntegrand<ContourZeroExponentImag> integrand(omega, sigma);
    ContourZeroExponentImag contour(omega, sigma);
    os << setw(16) << "omega.real" << ","
       << setw(16) << "omega.imag" << ","
       << setw(16) << "sigma" << ","
       << setw(15) << "s" << ","
       << setw(15) << "x" << ","
       << setw(15) << "y" << ","
       << setw(15) << "dx" << ","
       << setw(15) << "dy" << ","
       << setw(16) << "integrand.real" << ","
       << setw(16) << "integrand.imag" << endl;
    for (int i=0; i<n; ++i) {
      RealType t = (i+.5)/n;
      RealType s = (1-t)*contour.points().front() + t*contour.points().back();
      pair<complex<RealType>, complex<RealType> > z_dz_ds = contour.z_dz_ds(s);
      os << scientific
         << setw(16) << omega.real() << ","
         << setw(16) << omega.imag() << ","
         << setw(16) << sigma << ","
         << defaultfloat
         << setw(15) << s << ","
         << setw(15) << z_dz_ds.first.real() << ","
         << setw(15) << z_dz_ds.first.imag() << ","
         << setw(15) << z_dz_ds.second.real() << ","
         << setw(15) << z_dz_ds.second.imag() << ","
         << scientific
         << setw(16) << integrand(s).real() << ","
         << setw(16) << integrand(s).imag() 
         << defaultfloat << endl;
    }
  }
  
  
private:
  
  /// Compute the pdf for real argurments.  Used for characterietic function
  /// Note mu is assumed zero in this routine
  static complex<RealType> pdf(complex<RealType> x, RealType sigma) {
    const complex<RealType> one{1};
    return one/(exp(pow(log(x),2)/(2.*pow(sigma,2)))*
                        sqrt(2*pi<RealType>())*x*sigma);
  }
  
  /// class to generate contour that follows x axis
  class ContourZeroImag {
  public:
    ContourZeroImag(complex<RealType> omega, RealType sigma)
    :_omega(omega), _sigma(sigma), _points{-1, 1}{ }
    
    /// return y and dy_dx given the parameter s
    pair<complex<RealType>, complex<RealType> > z_dz_ds(RealType s) const{
      pair<complex<RealType>, complex<RealType> > ret;
      RealType x = _sigma * s/(1-s*s);
      RealType dx = _sigma * (1+s*s)/pow(1-s*s,2);
      
      ret.first = x;
      ret.second = dx;
      return ret;
    }
    
    vector<RealType> points() const {
      return _points;
    }
    
  private:
    complex<RealType> _omega;
    RealType _sigma;
    vector<RealType> _points;
    
  };

  /// contour from -inf + (pi/2-arg(omega)) *i to +inf + (pi/2-arg(omega)) * i
  
  class ContourPi2Imag {
  public:
    ContourPi2Imag(complex<RealType> omega, RealType sigma)
    :_shift(pi<RealType>()/2 - arg(omega)), // omega should be in the first quadrant
     _sigma(sigma), _points{-1, 1}{ }
    
    /// return y and dy_dx given the parameter s
    pair<complex<RealType>, complex<RealType> > z_dz_ds(RealType s) const{
      const complex<RealType> i{0,1};
      pair<complex<RealType>, complex<RealType> > ret;
      RealType x = _sigma * s/(1-s*s);
      RealType dx = _sigma * (1+s*s)/pow(1-s*s,2);
      
      ret.first = x + _shift * i;
      ret.second = dx;
      return ret;
    }
    
    vector<RealType> points() const {
      return _points;
    }
    
  private:
    complex<RealType> _shift;
    RealType _sigma;
    vector<RealType> _points;
    
  };
  
  class Contour;

  /// functor passed to newton raphson to find imag(exponent) == 0
  class ZeroExponentImagFunctor {
  public:
    ZeroExponentImagFunctor (complex<RealType> omega, RealType sigma, RealType x)
    : omega(omega), sigma(sigma), x(x){}
    
    /// return the imaginary part of exponent and it's derivative wrt y
    pair<RealType, RealType> operator() (RealType y) const {
      complex<RealType> i{0,1};
      pair<RealType, RealType> ret;
      ret.first = -(y*x)/(sigma*sigma)+(i*exp(x)*exp(i*y)*omega).imag();
      ret.second = -x/(sigma*sigma)-(exp(x)*exp(i*y)*omega).imag();
      return ret;
    }
    
  private:
    complex<RealType> omega;
    RealType sigma;
    RealType x;
  };
  
  /// class to generate contour that give an integrand exponent
  /// mostly zero
  class ContourZeroExponentImag {
  public:
    ContourZeroExponentImag(complex<RealType> omega, RealType sigma)
                 :_omega(omega), _sigma(sigma)
    {
      
      _xlower = - sigma * sqrt(abs(omega));
      _ylower = y_calc(_xlower);
      _xupper = sigma * sqrt(abs(omega));
      RealType yupper = y_calc(_xupper);
      _mid_dy_dx = (yupper-_ylower)/(_xupper-_xlower);
      _points.push_back(RealType(-1));
      _points.push_back(s_from_x(_xlower));
      _points.push_back(s_from_x(_xupper));
      _points.push_back(RealType(1));
    }
    
    /// return y and dy_dx given the parameter s
    pair<complex<RealType>, complex<RealType> > z_dz_ds(RealType s) const{
      const complex<RealType> i{0,1}, one{1,0};
      pair<complex<RealType>, complex<RealType> > ret;
      RealType x = _sigma * s/(1-s*s);
      RealType dx = _sigma * (1+s*s)/pow(1-s*s,2);

      if (x < _xlower || x > _xupper) {
        ZeroExponentImagFunctor f(_omega, _sigma, x);
        RealType y = y_calc(x);
        ret.first = x + i*y;
        ret.second = dx * (one+i*(-y/(_sigma*_sigma)+(i*exp(x+i*y)*_omega).imag())
                         / (x/(_sigma*_sigma)+(exp(x+i*y)*_omega).imag()));
      } else {
        ret.first = x + i * (_ylower+_mid_dy_dx * (x-_xlower));
        ret.second = dx * (one + i * _mid_dy_dx);
      }
      return ret;
    }
    
    vector<RealType> points() const {
      return _points;
    }
    
  private:
    RealType s_from_x(RealType x) {
      x /= _sigma;
      RealType ret;
      if (abs(x) < .25 ) {
        RealType term = x;
        RealType fac = pow(2 * x, 2);
        ret = term;
        for (int n=2; n<=20; ++n) {
          term *= (RealType(.5) - n)/n * fac;
          ret += term;
          if (abs(term) < std::numeric_limits<RealType>::epsilon() * std::abs(ret))
            break;
        }
      } else if (x < 0)
        ret = -1./(2.*x) - sqrt((4*pow(x,2)+1)/pow(x,2))/2.;
      else
        ret = -1./(2.*x) + sqrt((4*pow(x,2)+1)/pow(x,2))/2.;
      return ret;
      
    }
    
    /// return y satifying imag(exponent(x+iy) == 0
    RealType y_calc(RealType x) const {
      const RealType eps = std::numeric_limits<RealType>::epsilon();
      if (x>0 && exp(x)*real(_omega) > 1/eps)
        return pi<RealType>()/2;
      else if (x<0 && exp(x)*real(_omega) < -1/eps)
        return RealType{0};
      else {
        RealType guess = 0;
        const int digits = std::numeric_limits<RealType>::digits - 8;
        uintmax_t max_iter = 50;
        ZeroExponentImagFunctor f(_omega, _sigma, x);
        RealType ret = newton_raphson_iterate(f, guess,
                                              ylimits().first, ylimits().second,
                                              digits, max_iter);
        if (isnan(ret) || max_iter==50)
          throw std::runtime_error("ContourZeroExponentImag.y_calc: iteration failed");
        return ret;
      }
    }
    
    pair<RealType, RealType> ylimits() const {
      const RealType pi2 = pi<RealType>()/2;
      return make_pair<RealType>(-pi2+arg(_omega), pi2-arg(_omega));
    }
    
    complex<RealType> _omega;
    RealType _sigma;
    RealType _xlower;
    RealType _xupper;
    vector<RealType> _points;
    RealType _ylower;
    RealType _mid_dy_dx;

  };
  
  /// class of integrand used for integral in mixed domains
  template<typename Contour>
  class CFIntegrand {
  public:
    CFIntegrand(complex<RealType> omega,
                     RealType sigma,
                     int* neval = nullptr
                     ) : omega(omega), sigma(sigma),
                         contour(omega, sigma)
    {
      _neval = neval;
      if(neval) *neval=0;
    }
    
    /// return the integrand
    complex<RealType> operator() (RealType s) const {
      // map s from -1 to 1 to lnx_re domain -inf to +inf
      const complex<RealType> i{0.,1.};
      const RealType two{2};
      if (_neval) (*_neval)++;
      pair<complex<RealType>, complex<RealType> > z_dz_ds=contour.z_dz_ds(s);
      complex<RealType> ln_x = z_dz_ds.first;
      complex<RealType> dln_x = z_dz_ds.second;
      if ( !isfinite(exp(ln_x.real())) ) return complex<RealType>(0);
      RealType c = one_div_root_two_pi<RealType>() / sigma;
      complex<RealType> ln_pdf = log(c) + -pow(ln_x/sigma,2)/two + log(dln_x);
      complex<RealType> ln_fac = i * omega * exp(ln_x);
      complex<RealType> ret = exp(ln_fac + ln_pdf);
      if (isnan(real(ret)) || isnan(imag(ret))
          || !isfinite(real(ret)) || !isfinite(imag(ret))) {
        throw std::runtime_error("CFIntegrandMixed: invalid return");
      }
      return ret;
    }
    
    vector<RealType> points() const {return contour.points();}
    
  private:
    
    complex<RealType> omega;
    RealType sigma;
    Contour contour;
    int* _neval;
  };
  
  /// class wrapper for CFIntegrand to deliver real part
  template<typename Dist>
  class CFIntegrandRe {
  public:
    CFIntegrandRe(Dist& cf_integrand)
    : _cf_integrand(cf_integrand) {}
    
    RealType operator() (RealType s) const {return real(_cf_integrand(s));}
    
    void operator() (std::vector<RealType>& xs) const
    {
      for (RealType& x : xs)
        x = operator() (x);
    }

  private:
    Dist& _cf_integrand;
  };
  
  /// class wrapper for CFIntegrand to deliver imaginary part
  template<typename Dist>
  class CFIntegrandIm {
  public:
    CFIntegrandIm(Dist& cf_integrand)
    : _cf_integrand(cf_integrand) {}
    
    RealType operator() (RealType s) const {return imag(_cf_integrand(s));}

    void operator() (std::vector<RealType>& xs) const
    {
      for (RealType& x : xs)
        x = operator() (x);
    }

  private:
    Dist& _cf_integrand;
  };
  
  /// class of integrand used for integral in mixed domains
  template<typename Contour>
  class CFPrimeIntegrand {
  public:
    CFPrimeIntegrand(complex<RealType> omega,
                          RealType sigma,
                          int* neval = nullptr
                          ) : omega(omega), sigma(sigma),
                              contour(omega, sigma)
    {
      _neval=neval;
      if (neval) *_neval = 0;
    }
    
    /// return integrand for computation of derivative of characteristic fun.
    complex<RealType> operator() (RealType s) const {
      // map s from 0 to 1 to lnx domain of -inf+pi/2*i to inf+pi/2*i
      const complex<RealType> i{0.,1.};
      const RealType two{2};
      if (_neval) (*_neval)++;
      pair<complex<RealType>, complex<RealType> > z_dz_ds=contour.z_dz_ds(s);
      complex<RealType> ln_x = z_dz_ds.first;
      complex<RealType> dln_x = z_dz_ds.second;
      if ( !isfinite(exp(ln_x.real())) ) return complex<RealType>(0);
      RealType c = one_div_root_two_pi<RealType>() / sigma;
      complex<RealType> ln_pdf = log(c) -pow(ln_x/sigma,2)/two + log(dln_x);
      complex<RealType> ln_fac = i * omega * exp(ln_x); // should be negative real
      complex<RealType> ret = i * exp(ln_x + ln_fac + ln_pdf);
      if (isnan(real(ret)) || isnan(imag(ret))
          || !isfinite(real(ret)) || !isfinite(imag(ret))) {
        throw std::runtime_error("CFIntegrandMixed: invalid return");
      }
      return ret;
    }

    vector<RealType> points() const {return contour.points();}

  private:
    complex<RealType> omega;
    RealType sigma;
    Contour contour;
    int* _neval;
  };
  
  /// Return the characteristic function via Fourier integral using d countour
  template<typename Contour>
  complex<RealType> cf(complex<RealType> omega,  ///< the angular frequency
                       RealType* _error = nullptr, ///< [out] the integration error
                       RealType* _l1_norm = nullptr, ///< [out] the l1 norm
                       int* _neval = nullptr          ///< [out] the # of evaluations
  )
  {
    static complex<RealType> last_omega{0., -1000.};
    static complex<RealType> last_ret;
    if (omega == last_omega) {
      return last_ret;
    }
    complex<RealType> ret;
    RealType error_series;
    complex<RealType> ret_series = cf_series(omega, &error_series);
    const complex<RealType> i{0,1};
    if (error_series < _cf_ctl.epsrel*abs(ret_series) || error_series < _cf_ctl.epsabs) {
      ret = ret_series;
      if (_error) *_error = error_series;
      if (_l1_norm) *_l1_norm=std::numeric_limits<RealType>::quiet_NaN();
      if (_neval) *_neval = 0;
    } else if (abs(omega) == 0) {
      ret= static_cast<complex<RealType> >(1);
      if (_error) *_error=std::numeric_limits<RealType>::quiet_NaN();
      if (_l1_norm) *_l1_norm=std::numeric_limits<RealType>::quiet_NaN();
      if (_neval) *_neval = 0;
    } else {
      CFIntegrand<Contour> cf_integrand( omega*exp(mu()), sigma());
      
      RealType error_re, error_im;
      int neval_re, neval_im;
      typename IntegrationController<RealType>::TerminationCode t_c_re, t_c_im;
      int last_re, last_im;
      
      RealType ret_re, ret_im;
      {
        CFIntegrandRe<CFIntegrand<Contour> > cf_integrand_re(cf_integrand);
        _cf_ctl.integrate(cf_integrand_re, cf_integrand.points(),
                          ret_re, error_re, neval_re, t_c_re, last_re);
      }
      {
        CFIntegrandIm<CFIntegrand<Contour> > cf_integrand_im(cf_integrand);
        _cf_ctl.integrate(cf_integrand_im, cf_integrand.points(),
                          ret_im, error_im, neval_im, t_c_im, last_im);
      }
      ret = (ret_re + i * ret_im) * exp(-i * _raw_mean * omega);
      if (_error) *_error = error_re+error_im;
      if (_l1_norm) *_l1_norm = std::numeric_limits<RealType>::quiet_NaN();
      if (_neval) *_neval = neval_re + neval_im;
      
    }
    
    last_omega = omega;
    last_ret = ret;
    return ret;
  }
  
  /// Return the characteristic function via Fourier integral on contour
  template <typename Contour>
  complex<RealType> cfprime(complex<RealType> omega,  ///< the angular frequency
                            RealType* _error = nullptr, ///< [out] the integration error
                            RealType* _l1_norm = nullptr, ///< [out] the
                            int* _neval = nullptr     ///< the number of evaluations
  )
  {
    complex<RealType> ret;
    const complex<RealType> i{0,1};
    RealType error_series;
    complex<RealType> ret_series = cfprime_series(omega, &error_series);
    if (error_series < _cf_ctl.epsrel*abs(ret_series) || error_series < _cf_ctl.epsabs) {
      ret = ret_series;
      if (_error) *_error = error_series;
      if (_l1_norm) *_l1_norm=std::numeric_limits<RealType>::quiet_NaN();
      if (_neval) *_neval = 0;
    } else     if (abs(omega) == 0) {
      if (_error) *_error=std::numeric_limits<RealType>::quiet_NaN();
      if (_l1_norm) *_l1_norm=std::numeric_limits<RealType>::quiet_NaN();
      if (_neval) *_neval = 0;
      ret = complex<RealType> (0,mean());
    } else {
      CFPrimeIntegrand<Contour> cf_prime_integrand( omega*exp(mu()), sigma());
      
      RealType error_re, error_im;
      int neval_re, neval_im;
      typename IntegrationController<RealType>::TerminationCode t_c_re, t_c_im;
      int last_re, last_im;
      RealType ret_re, ret_im;
      
      {
        CFIntegrandRe<CFPrimeIntegrand<Contour> > cf_prime_integrand_re(cf_prime_integrand);
        _cf_ctl.integrate(cf_prime_integrand_re, cf_prime_integrand.points(),
                          ret_re, error_re, neval_re, t_c_re, last_re);
      }
      {
        CFIntegrandIm<CFPrimeIntegrand<Contour> > cf_prime_integrand_im(cf_prime_integrand);
        _cf_ctl.integrate(cf_prime_integrand_im, cf_prime_integrand.points(),
                          ret_im, error_im, neval_im, t_c_im, last_im);
      }
      ret = (ret_re + i * ret_im)*exp(mu());
      if (_error) *_error = (error_re+error_im)*exp(mu());
      if (_l1_norm) *_l1_norm = std::numeric_limits<RealType>::quiet_NaN();
      if (_neval) *_neval = neval_re + neval_im;
      ret = exp(-i * _raw_mean * omega) * ret
             -i * _raw_mean * cf<Contour>(omega);
    }
    return ret;
  }
  

  /// class of integrand used for Foureir integral in ln(x) domain
  class CFIntegrandLnX {
  public:
    CFIntegrandLnX(complex<RealType> omega,
                const boost::math::normal_distribution<RealType> nd
                ) : omega(omega), nd(nd),
                    _points{-1,1} {}
    complex<RealType> operator() (RealType s) const {
      RealType ln_x = nd.scale() * s/(1-s*s);
      RealType dln_x = nd.scale() * (1+s*s)/pow(1-s*s,2);
      const complex<RealType> i{0,1};
      complex<RealType> omega_exp_x = omega * exp(ln_x);
      if (!isfinite(real(omega_exp_x)) || !isfinite(imag(omega_exp_x)))
        return 0;
      else
        return exp(i* omega_exp_x)* boost::math::pdf(nd,ln_x)*dln_x;
    }
    
    vector<RealType> points() {return _points;}
    
  private:
    complex<RealType> omega;
    boost::math::normal_distribution<RealType> nd;
    vector<RealType> _points;
  };

  /// class of integrand used for Foureir integral in ln(x) domain
  class CFPrimeIntegrandLnX {
  public:
    CFPrimeIntegrandLnX(complex<RealType> omega,
                   const boost::math::normal_distribution<RealType> nd
                   ) : omega(omega), nd(nd),
                       _points{-1,1}  {}
    complex<RealType> operator() (RealType s) const {
      const complex<RealType> i{0.,1.};
      RealType ln_x = nd.scale() * s/(1-s*s);
      RealType dln_x = nd.scale() * (1+s*s)/pow(1-s*s,2);
      RealType x = exp(ln_x);
      if (!isfinite(x))
        return 0;
      else
        return i*x*exp(i*omega*x)* boost::math::pdf(nd,ln_x)*dln_x;
    }
    
    vector<RealType> points() {return _points;}
    
  private:
    complex<RealType> omega;
    boost::math::normal_distribution<RealType> nd;
    vector<RealType> _points;
  };
  
  /// class of integrand used for integral in mixed domains
  class CFIntegrandMixed {
  public:
    CFIntegrandMixed(complex<RealType> omega,
                   RealType sigma
                   ) : omega(omega), sigma(sigma),
                      _points{0,1} {}
    complex<RealType> operator() (RealType s) const {
      // map s from 0 to 1 to lnx domain of -inf+pi/2*i to inf+pi/2*i
      const complex<RealType> i{0.,1.};
      const RealType pi2 = pi<RealType>()/2.;
      const RealType two{2};
      RealType ln_x_re = (1./(1-s)-1/s);
      RealType dln_x_re = 1/pow(1-s,2)+1/pow(s,2);
      complex<RealType> ln_x = ln_x_re + i * pi2;
      complex<RealType> dln_x = dln_x_re;
      RealType exp_ln_x_re = exp(ln_x_re);
      if (!isfinite(exp_ln_x_re)) return complex<RealType>(0);
      RealType c = one_div_root_two_pi<RealType>() / sigma;
      complex<RealType> ln_pdf = log(c) + -pow(ln_x/sigma,2)/two + log(dln_x);
      complex<RealType> ln_fac = i * omega * exp(ln_x);
      complex<RealType> ret = exp(ln_fac + ln_pdf);
      if (isnan(real(ret)) || isnan(imag(ret))
          || !isfinite(real(ret)) || !isfinite(imag(ret))) {
        throw std::runtime_error("CFIntegrandMixed: invalid return");
      }
      return ret;
    }
    
    vector<RealType> points() {return _points;}
    
  private:
    complex<RealType> omega;
    RealType sigma;
    vector<RealType> _points;
  };
  
  /// class of integrand used for integral in mixed domains
  class CFPrimeIntegrandMixed {
  public:
    CFPrimeIntegrandMixed(complex<RealType> omega,
                     RealType sigma
                     ) : omega(omega), sigma(sigma),
                     _points{0., 1.}{}
    complex<RealType> operator() (RealType s) const {
      // map s from 0 to 1 to lnx domain of -inf+pi/2*i to inf+pi/2*i
      const complex<RealType> i{0.,1.};
      const RealType pi2 = pi<RealType>()/2.;
      const RealType two{2};
      RealType ln_x_re = (1./(1-s)-1/s);
      RealType dln_x_re = 1/pow(1-s,2)+1/pow(s,2);
      complex<RealType> ln_x = ln_x_re + i * pi2;
      complex<RealType> dln_x = dln_x_re;
      RealType exp_ln_x_re = exp(ln_x_re);
      if (!isfinite(exp_ln_x_re)) return complex<RealType>(0.);
      RealType c = one_div_root_two_pi<RealType>() / sigma;
      complex<RealType> ln_pdf = log(c) -pow(ln_x/sigma,2)/two + log(dln_x);
      complex<RealType> ln_fac = i * omega * exp(ln_x); // should be negative real
      complex<RealType> ret = i * exp(ln_x + ln_fac + ln_pdf);
      if (isnan(real(ret)) || isnan(imag(ret))
          || !isfinite(real(ret)) || !isfinite(imag(ret))) {
        throw std::runtime_error("CFIntegrandMixed: invalid return");
      }
      return ret;
    }
    
    vector<RealType> points() {return _points;}
  private:
    complex<RealType> omega;
    RealType sigma;
    vector<RealType> _points;
  };
  
  boost::math::lognormal_distribution<RealType> ln_dist;
  int _type;
  IntegrationController<RealType> _cf_ctl;
  RealType _raw_mean;
  vector<RealType> central_moments;

};


#endif /* lognormal_distribution_h */
