//
/// \file pareto_distribution.h
/// \package how_murh_data
//
/// \author Created by Joseph Dunn on 12/31/18.
/// \copyright © 2018 Joseph Dunn. All rights reserved.
//

#ifndef pareto_distribution_h
#define pareto_distribution_h

#include <iostream>
using std::cout;
using std::end;

#include <sstream>
using std::ostringstream;

#include <complex>
using std::complex;

#include <random>
using std::uniform_real_distribution;

#include <boost/math/quadrature/gauss_kronrod.hpp>
using boost::math::quadrature::gauss_kronrod;
#include <boost/math/constants/constants.hpp>
using boost::math::constants::euler;
#include <string>
using std::to_string;

/**
 * Instantiations of class template pareto_distribution model a
 * Pareto Type 2 distribution. Such a distribution produces random numbers
 * with \f$\displaystyle 1-F(x) = (1+\frac{x-\mu}{\sigma})^{-\alpha}\f$
 * for \f$ x > \mu \f$.
 *
 */
template<class RealType = double>
class pareto_distribution
{
public:
  typedef RealType result_type;
  
  class param_type
  {
  public:
    
    typedef pareto_distribution distribution_type;
    
    /** Constructs the parameters of a pareto_distribution. */
    explicit param_type(RealType alpha_arg,
                        RealType mu_arg = RealType(0.0),
                        RealType sigma_arg = RealType(1.0))
    : _alpha(alpha_arg), _mu(mu_arg), _sigma(sigma_arg) {}
    
    /** Returns the "alpha" parameter of the distribution. */
    RealType alpha() const { return _alpha; }
    
    /** Returns the "mu" parameter of the distribution. */
    RealType mu() const { return _mu; }
    
    /** Returns the "sigma" parameter of the distribution. */
    RealType sigma() const { return _sigma; }
    
    /** Writes the parameters to a std::ostream. */
    template <class charT, class traits>
    friend std::basic_ostream<charT,traits>& operator<< ( std::basic_ostream<charT,traits>& os,
                                             const param_type& parm ) {
      os << parm._alpha << " " << parm._mu << " " << parm._sigma;
      return os;
    }
    
    /** Reads the parameters from a std::istream. */
    template <class charT, class traits>
    friend std::basic_istream<charT,traits>& operator>> ( std::basic_istream<charT,traits>& is,
                                             param_type& parm )    {
      is >> parm._alpha >> std::ws >> parm._mu >> std::ws >> parm._sigma;
      return is;
    }
    
    /** Returns true if the two sets of parameters are equal. */
    friend bool operator== (const param_type& lhs, const param_type& rhs)
    { return lhs.alpha==rhs.alpha && lhs._mu == rhs._mu && lhs._sigma == rhs._sigma; }
    
    /** Returns true if the two sets of parameters are different. */
    friend bool operator!= (const param_type& lhs, const param_type& rhs) {
      return ! (lhs==rhs);
    }
    
  private:
    RealType _alpha;
    RealType _mu;
    RealType _sigma;
  };
  
  /**
   * Constructs a pareto_distribution. @c alpha @c mu and @c sigma are the
   * parameters of the distribution.
   */
  explicit pareto_distribution(RealType alpha_arg,
                               RealType mu_arg = RealType(0.0),
                               RealType sigma_arg = RealType(1.0))
  : _alpha(alpha_arg), _mu(mu_arg), _sigma(sigma_arg) {}
  
  /**
   * Constructs a pareto_distribution from its parameters.
   */
  explicit pareto_distribution(const param_type& parm)
  : _alpha(parm.alpha()), _mu(parm.mu()), _sigma(parm.sigma()) {}
  
  // compiler-generated copy ctor and assignment operator are fine
  
  /** Returns the alpha parameter of the distribution. */
  RealType alpha() const { return _alpha; }
  /** Returns the mu parameter of the distribution. */
  RealType mu() const { return _mu; }
  /** Returns the sigma parameter of the distribution. */
  RealType sigma() const { return _sigma; }
  
  /** Return the alpha of the asymptotic stable distribution */
  RealType alpha_stable() const {return std::min(RealType(2),alpha());}
  
  /** Returns the smallest value that the distribution can produce. */
  RealType min () const
  { return _mu; }
  /** Returns the largest value that the distribution can produce. */
  RealType max () const
  { return (std::numeric_limits<RealType>::infinity)(); }
  
  /** Returns the parameters of the distribution. */
  param_type param() const { return param_type(alpha(), mu(), sigma()); }

  /** Sets the parameters of the distribution. */
  void param(const param_type& parm)
  {
    _alpha=parm.alpha();
    _mu=parm.mu();
    _sigma = parm.sigma(); }
  
  /** the cdf of the distribution */
  RealType cdf(RealType x, bool lower_tail = true) const {
    if (x<_mu)
      return RealType(lower_tail ? 0 : 1);
    else {
      RealType uppercdf = pow((1+(x-_mu)/_sigma),-_alpha);
      return lower_tail ? 1 - uppercdf : uppercdf;
    }
  }
  
  /// return the pdf of the distribution
  RealType pdf(RealType x) const {
    if (x < -_mu)
      return RealType(0);
    else
      return _alpha * pow(1+(x-_mu)/_sigma,-_alpha-1)/_sigma;
  }
  
  /** the quantile for a given probability */
  RealType quantile(RealType p) const{
    return _mu + _sigma * (pow(1-p,-1/_alpha) - 1);
  }
  
  /** Return the mean of the distribution */
  RealType mean() const {return _mu + _sigma/(_alpha-1); }
  
  /** Return the MAD of the distribution */
  RealType mad() const {return 2 * pow(_alpha-1,_alpha-2)*pow(_alpha, 1-_alpha);}
  
  /** Return the MAD of the square of the distribution */
  RealType mad2() const {
    Mad2Integrand f(_alpha);
    return 2 * _sigma * gauss_kronrod<RealType,15>::integrate(f, 0, 2/(_alpha - 1));
  }
  
  /// return the characteristic function of the distribution
  complex<RealType> characteristic_function(RealType omega) const{
    if (omega == 0) return complex<RealType>(1,0);
    complex<RealType> arg(0,-_sigma*omega);
    return (_alpha*expint(1+_alpha, arg))/exp(complex<RealType>(0,_sigma*omega));
  }

  /** Return the 95% confidence interval */
  RealType ci(RealType level=RealType(.05)) const {
    return _sigma*(pow(level/2,-1/_alpha)-pow(1-level/2,-1/_alpha));
    
  }
  
  /**
   * Effects: Subsequent uses of the distribution do not depend
   * on values produced by any engine prior to invoking reset.
   */
  void reset() {}
  
  /**
   * Returns a random variate distributed according to the
   * Pareto distribution.
   */
  template<class Engine>
  result_type operator()(Engine& eng) const
  {
    std::uniform_real_distribution<RealType> dist;
    return _mu+_sigma*(pow(dist(eng),-1/_alpha)-1);
  }
  
  /**
   * Returns a random variate distributed according to the
   * Pareto distribution with parameters specified by param.
   */
  template<class Engine>
  result_type operator()(Engine& eng, const param_type& parm)
  { return pareto_distribution(parm)(eng); }
  
  /**
   * Write distribution to std::ostream
   */
  template <class charT, class traits>
  friend std::basic_ostream<charT,traits>& operator<< ( std::basic_ostream<charT,traits>& os,
                                                  const pareto_distribution& dist ) {
    os << "Pareto Type II with alpha = " << dist._alpha
       << ", mu = " << dist._mu << ", sigma =  " << dist._sigma << endl;
    os << "Mean   = " << dist.mean() << endl
    << "MAD    = " << dist.mad() << endl
    << "MAD2   = " << dist.mad2() << endl
    << "kappa2 = " << 2 - log(2)/(log(dist.mad2())-log(dist.mad())) << endl
    << "95% CI = " << dist.ci() << endl;
    return os;
  }

  /** Reads the parameters from a std::istream. */
  template <class charT, class traits>
  friend std::basic_istream<charT,traits>& operator>> ( std::basic_istream<charT,traits>& is,
                                                  pareto_distribution& dist )    {
    is >> dist._alpha >> std::ws >> dist._mu >> std::ws >> dist._sigma;
    return is;
  }
  
  /**
   * Returns true if the two distributions will produce identical
   * sequences of values given equal generators.
   */
  friend bool operator== (const pareto_distribution& lhs, const pareto_distribution& rhs)
  { return lhs._alpha==rhs._alpha && lhs._mu == rhs._mu && lhs._sigma == rhs._sigma; }
  
  /**
   * Returns true if the two distributions may produce different
   * sequences of values given equal generators.
   */
  friend bool operator!= (const pareto_distribution& lhs, const pareto_distribution& rhs) {
    return ! (lhs==rhs); }

private:
  /* return incomplete beta(-alpha,1-alpha) between x1 and x2 */
  static RealType adj_beta(RealType alpha, RealType x1, RealType x2){
      BetaIntegrand f(-alpha, 1-alpha);
      return gauss_kronrod<RealType,15>::integrate(f, x1, x2);
  }
     
  class Mad2Integrand {
    RealType a;
  public:
    Mad2Integrand (RealType alpha) : a(alpha) {}
    RealType operator() (RealType x) const {
      return 2 * a*a * pow(x+2,-2*a - 1) * (2/(a-1) - x) * adj_beta(a, 1/(x+2), (x+1)/(x+2));
    }
  };
  class BetaIntegrand {
    RealType a;
    RealType b;
  public:
    BetaIntegrand(RealType a, RealType b) : a(a), b(b) {}
    RealType operator() (RealType x) const {
      return pow(x,a-1)*pow(1-x, b-1);
    }
  };
  class ExpIntegrand {
    RealType n;
    complex<RealType> z;
  public:
    ExpIntegrand(RealType n, complex<RealType> z) : n(n), z(z) {}
    complex<RealType> operator() (RealType t) const {
      return exp(-z*t)/pow(t,n);
    }
  };
  
  struct ExpIntFraction{
    typedef std::pair<complex<RealType>, complex<RealType> > result_type;
    ExpIntFraction(RealType n_, complex<RealType> z_) : b(z_ + complex<RealType>(n_)), i(-1), n(n_) {}
    std::pair<complex<RealType>, complex<RealType> > operator()()
    {
      std::pair<complex<RealType>, complex<RealType> > result = std::make_pair(-static_cast<complex<RealType> >((i + 1) * (n + i)), b);
      b += 2;
      ++i;
      return result;
    }
  private:
    complex<RealType> b;
    int i;
    RealType n;
  };
  
  //
  // continued_fraction_b
  // Evaluates:
  //
  // b0 +       a1
  //      ---------------
  //      b1 +     a2
  //           ----------
  //           b2 +   a3
  //                -----
  //                b3 + ...
  //
  // Note that the first a0 returned by generator Gen is disarded.
  //
  
  inline static complex<RealType> continued_fraction_b(ExpIntFraction g,
                                                const RealType& factor,
                                                boost::uintmax_t& max_terms)
  {
    
    complex<RealType> tiny = 1.e-100;
    
    std::pair<complex<RealType>, complex<RealType> > v = g();
    
    complex<RealType> f, C, D, delta;
    complex<RealType> zero{0};
    complex<RealType> one{1};
    f = v.second;
    if(f == zero)
      f = tiny;
    C = f;
    D = 0;
    
    boost::uintmax_t counter(max_terms);
    
    do{
      v = g();
      D = v.second + v.first * D;
      if(D == zero)
        D = tiny;
      C = v.second + v.first / C;
      if(C == zero)
        C = tiny;
      D = one/D;
      delta = C*D;
      f = f * delta;
    }while((abs(delta - one) > factor) && --counter);
    
    max_terms = max_terms - counter;
    
    return f;
  }
  
  inline static std::complex<RealType> expint_as_fraction(RealType n, std::complex<RealType> const& z)
  {
    boost::uintmax_t max_iter = 1000;
    ExpIntFraction f(n, z);
    std::complex<RealType> result = continued_fraction_b(f,
                                                          std::numeric_limits<RealType>::epsilon(),
                                                          max_iter);
    result = exp(-z) / result;
    return result;
  }
  
  inline static std::complex<RealType> expint_as_series(RealType n, std::complex<RealType> const& z) {
    int max_iter = 1000;
    RealType s = 1-n;    // the first argument to the incomplete gamma function
    complex<RealType> ret;
    if (s>0 || floor(s) != s) {
      ret = pow(z,-s) * tgamma(s);
      complex<RealType> fac = 1;
      for (int i=0; i<max_iter; ++i){
        complex<RealType> term = fac / (s+i);
        ret -= term;
        if (abs(term) <= abs(ret) * std::numeric_limits<RealType>::epsilon()) break;
        fac *= -z/static_cast<std::complex<RealType> >(i+1);
      }
    } else {
      complex<RealType> gamma0 = -euler<RealType>()-log(z);
      complex<RealType> fac = 1;
      for (int i = 1; i<max_iter; ++i) {
        fac *= (-z)/static_cast<RealType>(i);
        complex<RealType> term = fac/static_cast<RealType>(i);
        gamma0 -= term;
        if (abs(term) <= abs(gamma0) * std::numeric_limits<RealType>::epsilon())
          break;
      }
      ret = pow(-1,-s) * gamma0 * pow(z,-s);
      fac = tgamma(-s);
      complex<double> sum = 0;
      for (int i=0; i< static_cast<int>(-s); ++i) {
        complex<RealType> term = fac;
        sum += term;
        if (i != -s-1) fac *= (-z)/(-s-i-1);
      }
      ret += exp(-z) * sum;
      ret /= tgamma(-s+1);
    }
    if (isnan(ret.real()) || isnan(ret.imag())){
      ostringstream oss;
      oss << "expint: nan encountered w n = " << n
      << " & z = " << z;
      throw std::runtime_error(oss.str());
    }
    return ret;
  }
  
   inline static std::complex<RealType> expint(RealType n, std::complex<RealType> const& z) {
     std::complex<RealType> ret =  (abs(z) < .5) ? expint_as_series(n,z)
                                                 : expint_as_fraction(n,z);
     if (isnan(ret.real()) || isnan(ret.imag())){
       ostringstream oss;
       oss << "expint: nan encountered w n = " << n
       << " & z = " << z;
       throw std::runtime_error(oss.str());
     }
     return ret;
   }
  
  RealType _alpha;
  RealType _mu;
  RealType _sigma;
};

#endif /* pareto_distribution_h */
