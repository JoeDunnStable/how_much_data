//
/// \file normal_switch_mean.h
/// \package how_much_data
//
/// \author Joseph Dunn on 3/3/19.
/// \copyright © 2019 Joseph Dunn. All rights reserved.
//

#ifndef normal_switch_mean_h
#define normal_switch_mean_h

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
#include <utility>
using std::pair;
using std::make_pair;

#include <boost/random.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;
using boost::math::constants::one_div_root_pi;
using boost::math::constants::root_half_pi;

#include <boost/math/tools/roots.hpp>
using boost::math::tools::newton_raphson_iterate;

/// class whose instances represent a 50/50 mixture of normal distriubtions
/// with sigma = 1 and means +d and -d
template<typename RealType=double>
class normal_switch_mean {
public:
  normal_switch_mean(RealType d      ///< [in] the means are +- d from 0
                     ) : _d(d), _switch(.5), _nd1(d,1), _nd2(-d,1),
                         _normal1(d,1), _normal2(-d,1) {}
  
  /// return a random number from the normalized distribution
  template <typename Engine>
  RealType operator() (Engine& eng)
  {
    if (_switch(eng))
        return _nd1(eng);
    else
        return _nd2(eng);
  }
  
  RealType min() const {
    return -std::numeric_limits<RealType>::infinity();
  }
  
  RealType max() const {
    return std::numeric_limits<RealType>::infinitye();
  }
  
  /// return the half difference betwee the means
  RealType d() const {return _d;}
  
  /// return the cdf
  RealType cdf(RealType x,             ///< [in] the desired quantile
               bool lower_tail=true   ///< [in] flag for tail
               ) const {
    using namespace boost::math;
    return lower_tail ? .5*(cdf(_normal1, x)+cdf(_normal2,x))
                      : .5*(cdf(complement(_normal1,x))+cdf(complement(_normal2,x)));
  }
  
  /// return the pdf
  RealType pdf(RealType x             ///< [in] the desired quantile
  ) const {
    using namespace boost::math;
    return .5*(pdf(_normal1,x)+pdf(_normal2,x));
  }
  
  /// return the quantile given the probability p
  RealType quantile(RealType p,           ///< [in] the disired probability
                    bool lower_tail=true ///< [in] flag for tail
                   ) const
  {
    using namespace boost::math;
    RealType guess = lower_tail ? .5*(boost::math::quantile(_normal1,p)
                                      +boost::math::quantile(_normal2,p))
                                : .5*(boost::math::quantile(complement(_normal1,p))+
                                      boost::math::quantile(complement(_normal2,p)));
    const int digits = std::numeric_limits<RealType>::digits - 8;
    uintmax_t max_iter = 50;
    QuantileFunctor f(_normal1, _normal2, lower_tail, p);
    const RealType inf = std::numeric_limits<RealType>::infinity();
    RealType ret = newton_raphson_iterate(f, guess,
                                          -inf, inf,
                                          digits, max_iter);
    if (boost::math::isnan(ret) || max_iter==50)
      throw std::runtime_error("quantile: iteration failed");
      return ret;
    
  }
  
  /// return the alpha parameter of the asymptotically equivalent stable distribution
  RealType alpha_stable() const { return static_cast<RealType>(2);}
  
  /// return the char. fun. for real arguments
  complex<RealType> characteristic_function(RealType omega
                                            ) const
  {
    return exp(-pow(omega,2)/2)*cos(_d * omega);
  }
  
  /// return the derivative of the char. fun. for real arguments
  complex<RealType> characteristic_function_prime(RealType omega
                                                  )const
  {
    return -exp(-pow(omega,2)/2)*(omega*cos(_d*omega)+ _d*sin(_d*omega));
  }
  
  RealType mean() const { return RealType(0);}
  
  RealType mad() const {
    return exp(-pow(_d,2)/2)/root_half_pi<RealType>()
           + _d * erf(_d/sqrt(RealType(2)));
  }
  
  RealType mad2() const {
    return (1+exp(-pow(_d,2)))*one_div_root_pi<RealType>()
           + _d * erf(_d);
    
  }
  
  /// Return the confidence interval of the distribution
  RealType ci(RealType level = RealType(.05)) const
  {
    return quantile(level/2, false)-quantile(level/2,true);
  }
  
  /// Write distribution to std::ostream
  friend ostream& operator<< ( ostream& os, normal_switch_mean& dist ) {
    os << "Normal Switching Mean w d = " << dist.d() << endl;
    os << "Mean   = " << dist.mean() << endl
    << "MAD    = " << dist.mad() << endl
    << "MAD2   = " << dist.mad2() << endl
    << "kappa2 = " << 2 - log(2)/(log(dist.mad2())-log(dist.mad())) << endl
    << "95% CI = " << dist.ci() << endl;
    return os;
  }
  
private:
  class QuantileFunctor {
  public:
    
    QuantileFunctor(const boost::math::normal_distribution<RealType>& dist1,
                    const boost::math::normal_distribution<RealType>& dist2,
                    bool lower_tail,
                    RealType p
                    ) : _dist1(dist1), _dist2(dist2) ,
                        _lower_tail(lower_tail), _p(p) {}
    
    pair<RealType, RealType> operator() (RealType x) const{
      if (_lower_tail)
        return make_pair(.5*(boost::math::cdf(_dist1,x)+boost::math::cdf(_dist2,x))-_p,
                         .5*(boost::math::pdf(_dist1,x)+boost::math::pdf(_dist2,x)));
      else
        return make_pair(.5*(boost::math::cdf(boost::math::complement(_dist1,x))
                             +boost::math::cdf(boost::math::complement(_dist2,x)))-_p,
                         -.5*(boost::math::pdf(_dist1,x)+boost::math::pdf(_dist2,x)));
    }

  private:
    const boost::math::normal_distribution<RealType>& _dist1;
    const boost::math::normal_distribution<RealType>& _dist2;
    bool _lower_tail;
    RealType _p;

  };
  RealType _d;
  boost::random::bernoulli_distribution<> _switch;
  boost::random::normal_distribution<RealType> _nd1;
  boost::random::normal_distribution<RealType> _nd2;
  boost::math::normal_distribution<RealType> _normal1;
  boost::math::normal_distribution<RealType> _normal2;
  
};


#endif /* normal_switch_mean_h */
