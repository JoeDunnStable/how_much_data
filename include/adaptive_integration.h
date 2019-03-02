///
/// \file  adaptive_integration.h
/// Routines for adapative integration ala QUADPACK
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef adaptive_integration_h
#define adaptive_integration_h
#include <iostream>
using std::cerr;
#include <vector>
#include <string>
using std::string;
//#include <cstddef>
#include "myFloat.h"
#include "gauss_kronrod.h"

namespace adaptive_integration {
  using namespace gauss_kronrod;
  
  using std::vector;
  using std::ostream;
  
  
  /// a subinterval in the Gauss Kronrod integreation process.
  template<typename myFloat>
  class Subinterval {
  public:
    myFloat a;              ///< the left endpoint
    myFloat b;              ///< the right endpoint
    myFloat r;              ///< the approximation to the integral of f
    myFloat e;              ///< the estimated error of the integration
    myFloat rabs;           ///< the appoximation to the integral of abs(f)
    myFloat defabs;         ///< the approximation of the integral of abs(f-mean(f))
    int level;             ///< the level of the subinterval.  Initial = 0
    int ndin;              ///< ndin = 1 for subintervals that are selected for initial valuation
    
    /// default subinterval with everything set to 0
    Subinterval() : a(0), b(0), r(0), e(0), ndin(0){};
    
    /// test whether the subinterval can be subdivided
    bool is_divisible() const{
      myFloat epmach = std::numeric_limits<myFloat>::epsilon();
      myFloat uflow = std::numeric_limits<myFloat>::min();
      return e!=0 && (std::max<myFloat>(fabs(a), fabs(b)) >
                      (1 + 100 * epmach) * (fabs((a+b)/2) + 1000 * uflow));
    }
    /// calculate the approximations to the integrals over the subinterval
    ///
    /// computes i = integral of f over (a,b), with error estimate
    ///          j = integral of abs(f) over (a,b)
    template<typename F>
    void integrate(F& f,                       ///< [in] the functor to be integrated
                   const int n_gauss,                 ///< [in]the number of Gauss nodes
                   const vector<myFloat>& w_gauss,    ///< [in] the n_gauss weights for the Gauss integration
                   const vector<myFloat>& x_kronrod,  ///< [in] the 2n_gauss+1 abscissae for the Kronrod integration.
                   ///< x_kronrod(1), x_kronrod(3), ...  abscissae of the n_gauss-point
                   ///< gauss rule.
                   ///< x_kronrod(0), x_kronrod(2), ...  abscissae which are optimally
                   ///< added to the n_gauss-point gauss rule.
                   const vector<myFloat>& w_kronrod,  ///< [in] the 2n_gauss+1 weights for the Kronrod integration
                   vector<myFloat>& fv                ///< [in,out] buffer for temporary storage of function values
    );
  };
  
  /// the data and comparison operator used to rank subintervals for subdivision
  template<typename myFloat>
  class SubintervalComp {
  public:
    int level_max;  ///< the cap on the level
    
    /// operator ordering subintervals for subdivision
    bool operator() (const Subinterval<myFloat> lhs, const Subinterval<myFloat> rhs) const {
      if ((!lhs.is_divisible()) && rhs.is_divisible()) return true;
      else if (lhs.is_divisible() && !rhs.is_divisible()) return false;
      else if (lhs.level == level_max && rhs.level < level_max) return true;
      else if (lhs.level < level_max && rhs.level == level_max) return false;
      else if (lhs.ndin < rhs.ndin) return true;
      else if (lhs.ndin > rhs.ndin) return false;
      else if (lhs.e < rhs.e)return true;
      else return false;
    }
    /// constructor for the subinterval comparison
    SubintervalComp(const int level_max ///< [in] cap on the level
    ) : level_max(level_max) {}
  };
  
  /// print out subintervals
  template<typename myFloat>
  void print_subs(ostream & os,                    ///< [in,out] the output stream to use
                  const vector<Subinterval<myFloat> >& subs, ///< [in] a reference to the vector of subintervals
                  const int last,                  ///< [in] the number of subintervals to print
                  const vector<myFloat>& points     ///< [in]the original endpoints before the subdivision process
  );
  
  /// print out a summary of the subintervals
  template<typename myFloat>
  void print_subs_summary(ostream & os,            ///< [in,out] the output stream to use
                          const vector<Subinterval<myFloat> >& subs, ///< [in] a reference to the vector of subintervals
                          const int last,                  ///< [in] the number of subintervals to print
                          const vector<myFloat>& points    ///< [in] the original endpoints before the subdivision process
  );
  
  template<typename myFloat> class IntegrationController;
  template<typename myFloat> ostream& operator<< (ostream& os, const IntegrationController<myFloat>& ctl);
  
  /// Integration controller for adaptive Gauss Kronrod integration over a
  /// finite interval given an initial subdivision.
  ///
  /// The routine calculates an approximation, \f$ result \f$, to a given definite integral
  /// \f$ i = \int_a^b f \f$, hopefully satisfying the following claim for accuracy
  /// \f$ |i-result| < \max(epsabs, epsrel*|i|)) \f$. Break points of the integration
  /// interval, where local difficulties of the integrand may occur(e.g. singularities,
  /// discontinuities), are provided by user.
  ///
  /// \author Joseph Dunn based on the original Fortran code of
  /// \author piessens,robert, appl. math. & progr. div. - k.u.leuven
  /// \author de doncker,elise, appl. math & progr. div. - k.u.leuven
  ///
  template<typename myFloat>
  class IntegrationController {
  private:
    bool  noext;               ///< flag indicating no extrapolation is to be used.
    Kronrod<myFloat> g_k;      ///< the number of Gauss nodes
    int limit;                 ///< the maximum number of subintervals
    int verbose;               ///< flag indicating verbose output
    vector<myFloat> fv;         ///< buffer to temporarily hold function values
  public:
    /// constructor for integration controller.  Nodes and weights are input
    template<typename BigFloat>
    IntegrationController(bool noext,       ///< flag indicating no extrapolation is to used
                          const Kronrod<BigFloat>& g_k_big, ///< A higher precision class of nodes & weights
                          myFloat epsabs,   ///< the target absolute error
                          myFloat epsrel,   ///< the target relative error
                          int limit,        ///< the maximum number of subintervals
                          int verbose       ///< flag indicating verbose output
    ) : noext(noext), g_k(g_k_big), limit(limit),
    verbose(verbose), fv(2*g_k_big.n_gauss+1), epsabs(epsabs), epsrel(epsrel),
    subs(limit)
    {
      myFloat epmach = std::numeric_limits<myFloat>::epsilon();
      if (epsabs==0 && epsrel < 50*epmach) {
        IntegrationController::epsrel = 50*epmach;
        cerr << "IntegrationController: Resetting epsrel to minimum allowable: " << IntegrationController::epsrel << endl;
      }
    }
    
    /// constructor using just myFloat.  Might be slightly inaccurate
    IntegrationController(bool noext,       ///< flag indicating no extrapolation is to used
                          int n_gauss,      ///< the number of Gauss nodes
                          myFloat epsabs,   ///< the target absolute error
                          myFloat epsrel,   ///< the target relative error
                          int limit,        ///< the maximum number of subintervals
                          int verbose       ///< flag indicating verbose output
    ) : IntegrationController(noext, Kronrod<myFloat>(n_gauss), epsabs, epsrel, limit, verbose){}
    
    myFloat epsabs;           ///< the target absolute error
    myFloat epsrel;           ///< the target relative error
    vector<Subinterval<myFloat> > subs; ///< work area for the integrator
    
    /// the termination codes returned by the integrator
    enum TerminationCode {normal, maximum_subdivisions, roundoff_error, bad_integrand,
      extrapolation_roundoff_error, divergent_integral, bad_input};
    
    /// integrate the function f over the subinterval using the adaptive integration
    template<typename F>
    void integrate(F& f,  ///< [in] the function to integrate
                   const vector<myFloat>& points,      ///< [in] the initial subdivsioin of the interval
                   myFloat& result,                    ///< [out] the approximation to the integral of f
                   myFloat& abserr,                    ///< [out] the estimated abs error of the approximation
                   int& neval,                         ///< [out] the number of function evaluations needed
                   TerminationCode& termination_code,  ///< [out] an error indicator.
                   int& last                           ///< [out] the number of subintervals actually used
    );
    
    /// return human readable termination code message
    static string tc_msg(TerminationCode tc) {
      const string msgs[] = {"OK","Maximum subdivisions reached","Roundoff error detected",
        "Bad integrand behavior","Roundoff error in the extrapolation table",
        "Integral probably divergent","Input is invalid"};
      return msgs[tc];
    }


    /// edge spacing
    myFloat edge_spacing() {return (1+g_k.x_kronrod.at(0))/2;}
    
    /// set the extrapolation indicator
    void set_noext(int next) { noext = next; }
    
    /// get the verbose indicator
    int get_verbose() {return verbose;}
    
    /// print out the control information
    friend ostream& operator<< <>(ostream& os, const IntegrationController<myFloat>& ctl);
  };
  
  /// A class with information needed to integrate a partiuclar Integrand
  template<typename myFloat, typename F>
  class Integral {
  private:
    F f;
    const vector<myFloat>& points;      ///< the initial subdivsioin of the interval
    
    IntegrationController<myFloat>& controller;  ///< the controller to use
    int verbose;             ///< the level of diagnostic output
    
  public:
    
    /// construct the functor
    Integral(F& f, ///< pointer to the integrand
             const vector<myFloat>& points, ///< A reference to the initial subdivision
             IntegrationController<myFloat>& ctl, ///< [in] pointer to the integration controller
             int verbose                    ///< the level of diagnostic output
    )
    : f(f), points(points), controller(ctl), verbose(verbose){}
    
    
     Integral(F&& f, ///< pointer to the integrand
             const vector<myFloat>& points, ///< A reference to the initial subdivision
             IntegrationController<myFloat>& ctl, ///< [in] pointer to the integration controller
             int verbose                    ///< the level of diagnostic output
    )
    : f(f), points(points), controller(ctl), verbose(verbose){}
    
    /// return a human readable error message
    string msg() {
      string msgs[] = {"OK","Maximum subdivisions reached","Roundoff error detected",
        "Bad integrand behavior","Roundoff error in the extrapolation table",
        "Integral probably divergent","Input is invalid"};
      return msgs[termination_code];};
    
    /// return a pointer to the instance of StandardStableDistribution
    myFloat result;        ///< the approximation of integral
    myFloat abserr;        ///< the estimated absolute error of the approximation
    int neval;            ///< the number of function evaluations used
    typename IntegrationController<myFloat>::TerminationCode termination_code; ///< the error indicator returned from IntegrationController::integrate
    int last;             ///< the number of subintervals required
    
    /// return the approximation to the integral of f(g(x))
    myFloat operator() ();
  };
  
} //namespace adaptive_integration

#include "adaptive_integration_impl.h"

#endif // adaptive_integration_h
