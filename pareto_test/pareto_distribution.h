//
//  pareto_distribution.h
//  pareto_test
//
//  Created by Joseph Dunn on 12/31/18.
//  Copyright © 2018 Joseph Dunn. All rights reserved.
//

#ifndef pareto_distribution_h
#define pareto_distribution_h

#include <random>
using std::uniform_real_distribution;

/**
 * Instantiations of class template pareto_distribution model a
 * Pareto Type 2 distribution. Such a distribution produces random numbers
 * with \f$\displaystyle 1-F(x) = (1+\frac{x-\mu}{\sigma})^{-\alpha}\f$
 * for x > 0.
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
    
    /** Constructs the parameters of a lognormal_distribution. */
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
   * Constructs a lognormal_distribution from its parameters.
   */
  explicit pareto_distribution(const param_type& parm)
  : _alpha(parm.alpha()), _mu(parm.mu()), _sigma(parm.sigma()) {}
  
  // compiler-generated copy ctor and assignment operator are fine
  
  /** Returns the alpha parameter of the distribution. */
  RealType alpha() const { return _alpha; }
  /** Returns the m parameter of the distribution. */
  RealType mu() const { return _mu; }
  /** Returns the s parameter of the distribution. */
  RealType sigma() const { return _sigma; }
  
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
  RealType cdf(RealType x) const {
    if (x<_mu)
      return RealType(0);
    else
      return 1 - pow((1+(x-_mu)/_sigma),-_alpha);
  }
  
  /** the quantile for a given probability */
  RealType quantile(RealType p) const{
    return _mu + _sigma * (pow(1-p,-1/_alpha) - 1);
  }
  
  /** Return the mean of the distribution */
  RealType mean() const {return _mu + _sigma/(_alpha-1); }
  
  /** Return the MAD of the distribution */
  RealType mad() const {return 2 * pow(_alpha-1,_alpha-2)*pow(_alpha, 1-_alpha);}
  
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
  RealType _alpha;
  RealType _mu;
  RealType _sigma;
};

#endif /* pareto_distribution_h */
