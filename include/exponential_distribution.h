//
/// \file exponential_distribution.h
/// \package how_much_data
//
/// \author  Created by Joseph Dunn on 1/3/19.
/// \copyright © 2019 Joseph Dunn. All rights reserved.
//

#ifndef exponential_distribution_h
#define exponential_distribution_h

#include <boost/random.hpp>
#include <boost/math/distributions/exponential.hpp>

/// instances of struct exponential_distribution generate random variates
/// for exponential distribution
template<class RealType = double>
struct exponential_distribution : public boost::random::exponential_distribution<RealType> {
  /// constructor give lambda
  exponential_distribution(RealType lambda) : boost::random::exponential_distribution<RealType>(lambda),
                                              exp_dist(lambda) {}
  
  /// return the cdf or the complement of the cdf
  RealType cdf(RealType x, bool lower_tail = true) const{
    if (x < 0)
      return lower_tail ? 0 : 1;
    else
      return lower_tail ? boost::math::cdf(exp_dist, x)
                        : boost::math::cdf(boost::math::complement(exp_dist, x));
  }
  
  /// return the pdf given x
  RealType pdf(RealType x) const{
    return boost::math::pdf(exp_dist, x);
  }
  
  /// return the quantile give the probability p
  RealType quantile(RealType p) const {
    return boost::math::quantile(exp_dist, p);
  }
  
  /// return the lambda paramter of the distribution
  RealType lambda() const { return boost::random::exponential_distribution<>::lambda();}
  
  /// return the alpha of the asymptotic stable distribution */
  RealType alpha_stable() const { return static_cast<RealType>(2);}
  
  /** Return mean of distribution */
  RealType mean() const {return boost::math::mean(exp_dist);}
  
  /** Return mean absolute deviation of the distribution */
  RealType mad() const {
    return 2 / (exp(1)*lambda());
  }
  
  /** Return the mad of the square of the distribution */
  RealType mad2() const {
    return 8 / (exp(2)*lambda());
  }
  
  /** Return the confidence interval of the distribution */
  RealType ci(RealType level = RealType(.05)) const
  {
    return boost::math::quantile(exp_dist,1-level/2)-boost::math::quantile(exp_dist,level/2);
  }
  
  /// return the characteristic function given real omega
  complex<RealType> characteristic_function(RealType omega) const{
    return lambda()/(lambda()-complex<RealType>(0,omega));
  }
                
  /// return the characteristic function given complex omega
  complex<RealType> characteristic_function(complex<RealType> omega) const{
    const complex<RealType> i{0.,1.};
    return lambda()/(lambda()-i * omega);
  }
  
  /// return the derivative of characteristic function given real omega
  complex<RealType> characteristic_function_prime(RealType omega) const{
    return complex<RealType>{0,lambda()}/pow(lambda()-complex<RealType>(0,omega),2);
  }
  
  /// return the derivative of characteristic function given complex omega
  complex<RealType> characteristic_function_prime(complex<RealType> omega) const{
    const complex<RealType> i{0.,1.};
    return i * lambda()/pow(lambda()-i * omega,2);
  }
  
  /**
   * Write distribution to std::ostream
   */
  template <class charT, class traits>
  friend std::basic_ostream<charT,traits>& operator<< ( std::basic_ostream<charT,traits>& os,
                                                       const exponential_distribution& dist ) {
    os << "Exponential distriubution w lambda = " << dist.lambda() << endl;
    os << "Mean   = " << dist.mean() << endl
    << "MAD    = " << dist.mad() << endl
    << "MAD2   = " << dist.mad2() << endl
    << "kappa2 = " << 2 - log(2)/(log(dist.mad2())-log(dist.mad())) << endl
    << "95% CI = " << dist.ci() << endl;
    return os;
  }
private:
  boost::math::exponential_distribution<RealType> exp_dist;
  
};


#endif /* exponential_distribution_h */
