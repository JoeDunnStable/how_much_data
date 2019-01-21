//
//  main.cpp
//  pareto_test
//
//  Created by Joseph Dunn on 12/31/18.
//  Copyright © 2018 Joseph Dunn. All rights reserved.
//

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
#include <iomanip>
using std::setw;
using std::setprecision;
using std::right;
using std::fixed;
using std::scientific;
using std::defaultfloat;
#include <string>
using std::string;
using std::to_string;
#include <sstream>
using std::istringstream;
using std::ostringstream;
#include <fstream>
using std::ofstream;
#include <random>
using std::mt19937;
#include <vector>
#include <array>
using std::array;
using std::vector;
#include <algorithm>
using std::max;
using std::min;
using std::max_element;
using std::sort;
#include <numeric>
using std::accumulate;
#include <utility>
using std::pair;
#include <mutex>
using std::mutex;
using std::unique_lock;
#include <thread>
using std::thread;

#include <complex>
using dcomplex = std::complex<double>;

#define EIGEN_FFTW_DEFAULT
#include <unsupported/Eigen/FFT>

Eigen::FFT<double> fft_eng;

#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;
using boost::timer::cpu_timer;
#include <boost/filesystem.hpp>
using boost::filesystem::path;
#include <boost/math/tools/roots.hpp>
using boost::math::tools::bracket_and_solve_root;
using boost::math::tools::eps_tolerance;
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;
using boost::math::constants::one_div_root_two_pi;
using boost::math::constants::one_div_two_pi;

#include "pareto_distribution.h"
#include "student_t_distribution.h"
#include "taleb_results.h"

struct Numpunct: public std::numpunct<char>{
protected:
  virtual char do_thousands_sep() const{return ',';}
  virtual std::string do_grouping() const{return "\03";}
};

template<typename RealType>
RealType rel_err(RealType a, RealType b) {
  return fabs(a-b)/std::max(a,b);
}

struct KappaResult {
  KappaResult (vector<int> ns) : ns(ns) {}
  double alpha;
  int nsize;
  double mad_rel_err;
  vector<int> ns;
  vector<double> mad;
  vector<double> kappa_mad;
  double ci_rel_err;
  vector<double> ci;
  vector<double> kappa_ci;
  void calc_kappa() {
    for (size_t i = 1; i<ns.size(); ++i) {
      kappa_mad.push_back(2- (log(ns.at(i))-log(ns.at(0)))/
                          (log(mad.at(i))-log(mad.at(0))));
      kappa_ci.push_back(2- (log(ns.at(i))-log(ns.at(0)))/
                         (log(ci.at(i))-log(ci.at(0))));
    }
  }
};

ostream& operator<< (ostream& os, const KappaResult& k) {
  os << fixed << setw(7) << setprecision(2) << k.alpha
  << setw(14) << k.nsize
  << setw(12) << setprecision(2) << k.mad_rel_err*100 << "%";
  for (size_t i=1; i<k.mad.size(); ++i) {
    os << setw(13) << setprecision(3) << k.kappa_mad.at(i-1);
  }
  os << setw(12) << setprecision(2) << k.ci_rel_err*100 << "%";
  for (size_t i=1; i<k.ci.size(); ++i) {
    os << setw(13) << setprecision(3) << k.kappa_ci.at(i-1);
  }
  os << endl;
  return os;
}

struct KappaResults {
  KappaResults(const vector<int>& ns, size_t n_alphas, size_t taleb_offset)
  : ns(ns), kr(n_alphas, ns), taleb_offset(taleb_offset) {}
  KappaResult& at(size_t i) {
    unique_lock<mutex> lock(kr_mutex);
    return this->kr.at(i);
  }
  vector<int> ns;
  vector<KappaResult> kr;
  size_t taleb_offset;
  mutex kr_mutex;
};

ostream& operator<< (ostream& os, KappaResults& ks) {
  os << setw(7) << right << "alpha"
     << setw(14) << right << "nsize"
     << setw(13) << right << "mad_rel_err";
  for (int i=1; i<ks.ns.size(); ++i)
    os << setw(13) << right << "kappa"+to_string(ks.ns.at(i))+"_mad";
  os << setw(13) << right << "ci_rel_err";
  for (int i=1; i<ks.ns.size(); ++i)
    os << setw(13) << right << "kappa"+to_string(ks.ns.at(i))+"_ci";
  os << endl << endl;
  for (auto kr : ks.kr)
    os << kr;
  os << endl;
  os << setw(68) << right << "Error vs Taleb's Results" << endl << endl;
  for (size_t i=0; i<ks.kr.size(); ++i) {
    os << setw(7) << setprecision(2) << ks.kr.at(i).alpha << setw(27) << " ";
    for (size_t j=0; j<ks.ns.size()-1; ++j) {
      double kappa_mad = ks.kr.at(i).kappa_mad.at(j);
      double err = kappa_mad - taleb_results.at(i).at(j+ks.taleb_offset);
      os << setw(13) << setprecision(3) << err;
    }
    os << endl;
  }
  os << endl;
  return os;
}

int mod (int a, int b)
{
  int ret = a % b;
  if(ret < 0)
    ret+=b;
  return ret;
}

double confidence_interval(int nmin, int nmax, double delta,
                           const vector<double>& pdf, const vector<double>& x,
                           double ci_level) {
  int nsize = static_cast<int>(pdf.size());
  double ci_lower_left{x.at(mod(nmin,nsize))};
  double ci_lower_right=ci_lower_left;
  double cdf_lower_left=.5*pdf.at(mod(nmin,nsize))*delta;
  double cdf_lower_right = cdf_lower_left;
  for (int j=nmin; j<=nmax ; ++j) {
    if (cdf_lower_right > ci_level/2) break;
    ci_lower_left = ci_lower_right;
    cdf_lower_left = cdf_lower_right;
    ci_lower_right= x.at(mod(j+1, nsize));
    cdf_lower_right = cdf_lower_left+.5*(pdf.at(mod(j, nsize))+pdf.at(mod(j+1,nsize)))*delta;
  }
  double ci_lower;
  if (ci_lower_left == ci_lower_right)
    ci_lower = ci_lower_right;
  else {
    double t = (ci_level/2-cdf_lower_left)/(cdf_lower_right-cdf_lower_left);
    ci_lower = (1-t) * ci_lower_left + t * ci_lower_right;
  }
  
  double ci_upper_left{x.at(mod(nmax,nsize))};
  double ci_upper_right=ci_upper_left;
  double cdf_upper_left{.5*pdf.at(mod(nmax,nsize))*delta};
  double cdf_upper_right=cdf_upper_left;
  for (int j=nmax; j>=nmin ; --j) {
    if (cdf_upper_left > ci_level/2) break;
    ci_upper_right = ci_upper_left;
    cdf_upper_right = cdf_upper_left;
    ci_upper_left= x.at(mod(j-1, nsize));
    cdf_upper_left = cdf_upper_right+.5*(pdf.at(mod(j,nsize))+pdf.at(mod(j-1, nsize)))*delta;
  }
  double ci_upper;
  if (ci_upper_right == ci_upper_left)
    ci_upper = ci_upper_left;
  else {
    double t = (ci_level/2-cdf_upper_right)/(cdf_upper_left-cdf_upper_right);
    ci_upper = (1-t) * ci_upper_right + t * ci_upper_left;
  }
  return ci_upper-ci_lower;
}

/// functor for use in finding upper limit w target mad
template<typename Dist>
class Upper {
public:
  Upper (double delta, Dist& dist) : delta(delta), dist(dist) {}
  double operator() (double x) {
    double a = dist.alpha()/(dist.alpha()-1);
    double ret = fabs(x-dist.mean())* a * (dist.cdf(x, false))-delta;
    return ret;
  }
private:
  double delta;
  Dist& dist;
};

template<typename Dist>
class Lower {
public:
  Lower (double delta, Dist&
         dist) : delta(delta), dist(dist) {}
  double operator() (double x) {
    double a = dist.alpha()/(dist.alpha()-1);
    double ret = fabs(x-dist.mean())* a * (dist.cdf(x, true))-delta;
    return ret;
  }
private:
  double delta;
  Dist& dist;
};

/// Check whether the factorization of n contains only primes
bool factor_check(int n,              ///< the number to check
                  vector<int> primes  ///< the array of candidate primes
                  )
{
  int n_current = n;
  for (auto p : primes) {
    while (n_current % p == 0 ) {
      n_current = n_current / p;
    }
  }
  return n_current == 1;
}

struct Job {
  int nmin;
  int nmax;
  int ncurrent;
  int nchunk;
  mutex job_mutex;
  Job(int nmin, int nmax, int nchunk) : nmin(nmin), nmax(nmax), ncurrent(nmin),
                                        nchunk(nchunk) {}
  bool get_next(int& nstart, int& nend) {
    unique_lock<mutex> lock;
    if (ncurrent > nmax) return false;
    nstart = ncurrent;
    ncurrent = min(nmax+1,ncurrent+nchunk);
    nend = ncurrent;
    return true;
  }
};

template <typename Dist>
void calc_characteristic_function(const Dist& dist,
                                const double mean,
                                const int n,
                                const double delta_omega,
                                Job *job,
                                vector<dcomplex> *adj_cf) {
  double alpha_stable = min(2., dist.alpha());
  int nsize = static_cast<int>(adj_cf->size());
  int n_start, n_end;
  while (job->get_next(n_start, n_end)) {
    for (int i=n_start; i<n_end; ++i) {
      double omega = delta_omega*i/pow(n,1/alpha_stable);
      // remove the location parameter from the characteristic function
      dcomplex fac = exp(dcomplex(0,-mean*omega));
      dcomplex cf = pow(fac*dist.characteristic_function(omega),n);
      adj_cf->at(mod(i,nsize)) = cf;
      if (i !=0)
        adj_cf->at(mod(-i,nsize)) = conj(cf);
    }
  }
}

template<typename Dist>
void calculate_kappa(double delta,
                     double delta2,
                     int m,
                     vector<int> ns,
                     Dist dist,
                     double ci_level,
                     KappaResult& k,
                     bool verbose = false) {
  double mean = dist.mean();
  double alpha = dist.alpha();
  double alpha_stable = min(2.,alpha);
  // Calculate xmax and xmin
  // so that the contriubution of the tails to MAD is small
  // for upper the contribution is approximately |xmax-mean| * (1-CDF(xmax))
  if (verbose)
    cout << "alpha = " << alpha << endl
         << "mean  = " << mean << endl;
  Upper<Dist> upper(delta2, dist);
  double guess = dist.quantile(1-delta2/2);
  double factor = 2;
  bool rising = false;
  eps_tolerance<double> tol;
  boost::uintmax_t max_iter{1000};
  pair<double, double> root = bracket_and_solve_root(upper, guess, factor,
                                                     rising, tol, max_iter);
  double xmax0 = root.second - mean;
  
  // for lower the contribution is approximately |xmin-mean| * CDF(xmin)
  Lower<Dist> lower(delta2, dist);
  guess = dist.quantile(delta2/2);
  rising = true;
  max_iter = 1000;
  root = bracket_and_solve_root(lower, guess, factor, rising, tol, max_iter);
  double xmin0 = root.first - mean;
  double xmax = max(fabs(xmax0),fabs(xmin0));
  
  int nmax = static_cast<int>(min(double(m),xmax/delta));
  if (nmax < 1000000) {
    // Split the extra points between increasing the range and the density
    delta = delta/sqrt(1000000./nmax);
    nmax = 1000000;
  }
  int nsize = 2*nmax+1;
  
  // fftw works best for array sizes that are products of small primes
  vector<int> primes{3,5,7,11};
  for (; !factor_check(nsize, primes); nsize+=2) {}
  nmax = ((nsize-1)/2);
  int nmin = -nmax;

  xmax = nmax * delta;
  double uppermad = delta2 + upper(xmax);
  double lowermad = delta2 + lower(-xmax);
  k.nsize = nsize;
  if (verbose) {
    cout << setw(12) << "xmin0 = " << setw(15) << setprecision(2) << xmin0 << endl
         << setw(12) << "xmax0 = " << setw(15) << setprecision(2) << xmax0 << endl
         << setw(12) << "nmin = " << setw(15) << nmin << endl
         << setw(12) << "nmax = " << setw(15) << nmax << endl
         << setw(12) << "xmin = " << setw(15) << setprecision(2) << -xmax << endl
         << setw(12) << "xmax = " << setw(15) << setprecision(2) << xmax << endl
         << setw(12) << "mad = " << setw(15) << setprecision(4) << dist.mad() << endl
         << scientific
         << setw(12) << "uppermad = " << setw(15) << setprecision(2) << uppermad << endl
         << setw(12) << "lowermad = " << setw(15) << setprecision(2) << lowermad << endl
         << defaultfloat
         << setw(10) << "nsize = " << setw(15) << nsize << endl << endl;
  }
  vector<double> x(nsize,0.);
  for (int i=-nmax; i<=nmax; ++i)
    x.at(mod(i,nsize)) = delta * i;
  double omega_max = pi<double>()/delta;
  double delta_omega = (2*omega_max/nsize);
  vector<dcomplex> adj_cf(nsize);
  if (verbose) {
    cout << setw(5) << " "
    << setw(13) << right << "x.at(nmin)"
    << setw(13) << right << "x.at(-1)"
    << setw(13) << right << "x.at(0)"
    << setw(13) << right << "x.at(1)"
    << setw(13) << right << "x.at(nmax)"
    << endl;
  
    cout << setw(5) << " "
    << setw(13) << fixed << setprecision(3) << x.at(mod(nmin,nsize))
    << setw(13) << fixed << setprecision(3) << x.at(mod(-1,nsize))
    << setw(13) << fixed << setprecision(3) << x.at(mod(0,nsize))
    << setw(13) << fixed << setprecision(3) << x.at(mod(1,nsize))
    << setw(13) << fixed << setprecision(3) << x.at(mod(nmax,nsize))
    << endl << endl;
    cout << setw(5) << right << "n"
    << setw(13) << right << "pdf.at(nmin)"
    << setw(13) << right << "pdf.at(-1)"
    << setw(13) << right << "pdf.at(0)"
    << setw(13) << right << "pdf.at(1)"
    << setw(13) << right << "pdf.at(nmax)"
    << setw(13) << right << "p_total"
    << setw(9) << right << "mad"
    << setw(9) << right << "madlower"
    << setw(9) << right << "madupper"
    << endl;
  }
  for (size_t j=0; j<ns.size(); ++j){
    int n = ns.at(j);
    array<thread,8> threads;
    int nchunk = 1000;
    Job job(nmin, nmax, nchunk);
    for (int i=0; i<8; ++i) {
      threads.at(i)=thread{calc_characteristic_function<Dist>,
                           dist, mean, n, delta_omega,
                           &job, &adj_cf};
    }
    for (int i=0; i<8; ++i)
      threads.at(i).join();
    vector<double> pdf(nsize,0.);
    fft_eng.inv(pdf, adj_cf);
    for (int i=-nmax; i<=nmax; ++i)
      pdf.at(mod(i,nsize)) *= delta_omega * nsize * one_div_two_pi<double>();
    double mad = 0;
    for (int i=-nmax; i<=nmax; ++i) {
      mad += pdf.at(mod(i,nsize)) * fabs(x.at(mod(i,nsize)))*delta;
    }
    double mad_upper_tail = mad*uppermad/dist.mad();
    double mad_lower_tail = mad*lowermad/dist.mad();;
    mad += mad_upper_tail+mad_lower_tail;
    mad *= pow(n,1/alpha_stable);
    k.mad.push_back(mad);
    double ci = confidence_interval(nmin, nmax, delta, pdf, x, ci_level);
    ci *= pow(n,1/alpha_stable);
    k.ci.push_back(ci);
    if (verbose) {
      double p_total = accumulate(pdf.begin(),pdf.end(),0.)*delta;
      cout << setw(5) << n
      << setw(13) << fixed << setprecision(8) << pdf.at(mod(nmin,nsize))
      << setw(13) << fixed << setprecision(8) << pdf.at(mod(-1,nsize))
      << setw(13) << fixed << setprecision(8) << pdf.at(mod(0,nsize))
      << setw(13) << fixed << setprecision(8) << pdf.at(mod(1,nsize))
      << setw(13) << fixed << setprecision(8) << pdf.at(mod(nmax,nsize))
      << setw(13) << fixed << setprecision(8) << p_total
      << setw(9) << fixed << setprecision(4) << mad
      << setw(9) << fixed << setprecision(4) << mad_lower_tail*pow(n,1/alpha_stable)
      << setw(9) << fixed << setprecision(4) << mad_upper_tail*pow(n,1/alpha_stable)
      << endl;
    }
  } // j over ns
  k.calc_kappa();
  k.mad_rel_err = rel_err(k.mad.at(0), dist.mad());
  k.ci_rel_err = rel_err(k.ci.at(0), dist.ci(ci_level));
}

void show_usage(path p) {
  cerr << "Usage: " << p.filename().string() << " delta [delta2=.001] [m=1e7]" << endl
       << "  where delta is the interval of integration" << endl
       << "  and delta2 is the tail contribution to mad" << endl;
}


int main(int argc, const char * argv[]) {
  path p(argv[0]);
  if ( (argc < 2) || (argc > 4) ) {
    show_usage(p);
    return 1;
  }
  istringstream iss1{string(argv[1])};
  double delta;
  if ( !(iss1 >> delta) || ( delta <=0 ))
    show_usage(p);
  
  double delta2 = .001;
  if (argc >= 3) {
    istringstream iss2{string(argv[2])};
    if ( !(iss2 >> delta2) || (delta2 <=0 ) )
      show_usage(p);
  }
  int m = 1e7;
  if (argc >= 4) {
    istringstream iss3{string(argv[3])};
    if ( !(iss3 >> m) || (m <=0 ) )
      show_usage(p);
  }
  cout.imbue({std::locale(), new Numpunct});

  vector<double> alphas;
  
  for (size_t i=0; i<taleb_results.size(); ++i)
    alphas.push_back(taleb_results.at(i).at(0));
  
  vector<int> ns{1, 2, 30, 100};
  KappaResults ks_pareto(ns, alphas.size(),1);
  KappaResults ks_student(ns, alphas.size(),4);
  
  double ci_level = .05;
  
  string out_dir{"../output"};
  ostringstream oss;
  oss << out_dir << "/convolution_test_"
      << delta << "_"
      << delta2 << "_"
      << m
      << ".out";
  cout << "Writing to file " << oss.str() << endl << endl;
  cout.flush();
  
  ofstream out{oss.str()};
  out.imbue({std::locale(), new Numpunct});
  auto_cpu_timer t(out);
  for (size_t i=0; i<alphas.size(); ++i) {
    double alpha=alphas.at(i);
    pareto_distribution<> pd(alpha);
    KappaResult kr(ns);
    kr.alpha = alpha;
    calculate_kappa(delta, delta2, m, ns, pd, ci_level, kr, true);
    ks_pareto.at(i) = kr;
  }
  out << "Pareto Distribution" << endl << endl;
  out << ks_pareto;
  out.flush();
  
  for (size_t i=0; i<alphas.size(); ++i) {
    double alpha = alphas.at(i);
    student_t_distribution<> td(alpha);
    KappaResult kr(ns);
    kr.alpha = alpha;
    calculate_kappa(delta, delta2, m, ns, td, ci_level, kr, true);
    ks_student.at(i)=kr;
  }
  out << "Student Distribution" << endl << endl;
  out << ks_student;

    return 0;
}
