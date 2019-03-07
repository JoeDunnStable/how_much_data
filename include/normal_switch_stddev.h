//
/// \file normal_switch_stddev.h
/// \package how_much_data
//
/// \author Joseph Dunn on 3/3/19.
/// \copyright © 2019 Joseph Dunn. All rights reserved.
//

#ifndef normal_switch_stddev_h
#define normal_switch_stddev_h

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

/// class whose instances represent a mixture of normal distriubtions
/// with mu=0 and stddev of 1 w prob. (1-p) and a with prob. p
template<typename RealType=double>
class normal_switch_stddev {
public:
  normal_switch_stddev(RealType a,      ///< [in] the means are +- d from 0
                       RealType p       ///< [in] the prob. of sigma=a
                     ) : _a(a), _p(p), _switch(p), _nd1(0,1), _nd2(0,a),
                         _normal1(0,1), _normal2(0,a) {}
  
  /// return a random number from the normalized distribution
  template <typename Engine>
  RealType operator() (Engine& eng)
  {
    if (_switch(eng))
        return _nd2(eng);
    else
        return _nd1(eng);
  }
  
  RealType min() const {
    return -std::numeric_limits<RealType>::infinity();
  }
  
  RealType max() const {
    return std::numeric_limits<RealType>::infinitye();
  }
  
  /// return the second stddev
  RealType a() const {return _a;}
  
  /// return the prob. of a switch to a
  RealType p() const {return _p;}
  
  /// return the cdf
  RealType cdf(RealType x,             ///< [in] the desired quantile
               bool lower_tail=true   ///< [in] flag for tail
               ) const {
    using namespace boost::math;
    return lower_tail ? (1-_p)*cdf(_normal1, x)+ _p*cdf(_normal2,x)
                      : (1-_p)*cdf(complement(_normal1,x))+_p*cdf(complement(_normal2,x));
  }
  
  /// return the pdf
  RealType pdf(RealType x             ///< [in] the desired quantile
  ) const {
    using namespace boost::math;
    return (1-_p)*pdf(_normal1,x)+_p*pdf(_normal2,x);
  }
  
  /// return the quantile given the probability p
  RealType quantile(RealType pp,           ///< [in] the disired probability
                    bool lower_tail=true ///< [in] flag for tail
                   ) const
  {
    using namespace boost::math;
    RealType q1 = lower_tail ? boost::math::quantile(_normal1, pp)
                             : boost::math::quantile(complement(_normal1,pp));
    RealType q2 = lower_tail ? boost::math::quantile(_normal2, pp)
                             : boost::math::quantile(complement(_normal2,pp));
    RealType guess = (1-_p) * q1 + _p * q2;
    const int digits = std::numeric_limits<RealType>::digits - 8;
    uintmax_t max_iter = 50;
    QuantileFunctor f(_normal1, _normal2, _p, lower_tail, pp);
    RealType ret = newton_raphson_iterate(f, guess,
                                          std::min(q1,q2), std::max(q1,q2),
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
    return (1-_p)*exp(-pow(omega,2)/2) + _p * exp(-pow(_a*omega,2)/2);
  }
  
  /// return the derivative of the char. fun. for real arguments
  complex<RealType> characteristic_function_prime(RealType omega
                                                  )const
  {
    return -omega*((1-_p)*exp(-pow(omega,2)/2)+pow(_a,2)*_p*exp(-pow(_a*omega,2)/2));
  }
  
  RealType mean() const { return RealType(0);}
  
  RealType mad() const {
    return (1-_p + _a *_p)/root_half_pi<RealType>();
  }
  
  RealType mad2() const {
    return 2*(1-_p*(2-sqrt(2+2*pow(_a,2))+_p*(-1 - _a + sqrt(2+2*pow(_a,2)))))
            *one_div_root_pi<RealType>();
    
  }
  
  /// Return the confidence interval of the distribution
  RealType ci(RealType level = RealType(.05)) const
  {
    return quantile(level/2, false)-quantile(level/2,true);
  }
  
  /// Write distribution to std::ostream
  friend ostream& operator<< ( ostream& os, normal_switch_stddev& dist ) {
    os << "Normal Switching stddev w a = " << dist.a() << " w prob p ="<< dist.p() << endl;
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
                    const RealType p,
                    bool lower_tail,
                    RealType pp
                    ) : _dist1(dist1), _dist2(dist2) , _p(p),
                        _lower_tail(lower_tail), _pp(pp) {}
    
    pair<RealType, RealType> operator() (RealType x) const{
      if (_lower_tail)
        return make_pair((1-_p)*boost::math::cdf(_dist1,x)+_p*boost::math::cdf(_dist2,x)-_pp,
                         (1-_p)*boost::math::pdf(_dist1,x)+_p*boost::math::pdf(_dist2,x));
      else
        return make_pair((1-_p)*boost::math::cdf(boost::math::complement(_dist1,x))
                             +_p*boost::math::cdf(boost::math::complement(_dist2,x))-_pp,
                         -((1-_p)*boost::math::pdf(_dist1,x)+_p*boost::math::pdf(_dist2,x)));
    }

  private:
    const boost::math::normal_distribution<RealType>& _dist1;
    const boost::math::normal_distribution<RealType>& _dist2;
    const RealType _p;
    bool _lower_tail;
    RealType _pp;

  };
  RealType _a;
  RealType _p;
  boost::random::bernoulli_distribution<> _switch;
  boost::random::normal_distribution<RealType> _nd1;
  boost::random::normal_distribution<RealType> _nd2;
  boost::math::normal_distribution<RealType> _normal1;
  boost::math::normal_distribution<RealType> _normal2;
  
};


#endif /* normal_switch_mean_h */
