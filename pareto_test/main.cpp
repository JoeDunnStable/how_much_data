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
using std::vector;
#include <algorithm>
using std::max;
using std::min;
using std::max_element;
using std::sort;
#include <numeric>
using std::accumulate;
#include <mutex>
using std::mutex;
using std::unique_lock;

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
#include <boost/math/special_functions/beta.hpp>
using boost::math::beta;
#include <boost/math/distributions/students_t.hpp>
using boost::math::students_t;

#include "pareto_distribution.h"
#include "student_t_distribution.h"

template<typename RealType>
RealType rel_err(RealType a, RealType b) {
  return fabs(a-b)/std::max(a,b);
}

template<typename myFloat>
vector<myFloat> quantile(const vector<myFloat> &x, const vector<myFloat>& probs)
{
  myFloat eps = 100 * std::numeric_limits<myFloat>::epsilon();
  long np = probs.size();
  for (auto p : probs) {
    if (p < -eps || p > 1 + eps)
      throw std::range_error("quantile: 'probs' outside [0,1]");
  }
  vector<myFloat> qs(np);
  long n = x.size();
  if (n > 0 && np > 0) {
    vector<myFloat> x_sort = x;
    sort(x_sort.data(), x_sort.data()+n);
    for (int j=0; j<np; ++j) {
      myFloat index = (n - 1) * probs.at(j);
      int lo = static_cast<int>(floor(index));
      int hi = static_cast<int>(ceil(index));
      myFloat h = index - lo;
      qs.at(j) = (1-h) * x_sort.at(lo) + h * x_sort.at(hi);
    }
    return qs;
  } else {
    throw std::range_error("quantile: Both x and prob must be of length greater than 0");
  }
}

struct KappaResult {
  double alpha;
  double mad_rel_err;
  vector<double> kappa_mad;
  double ci_rel_err;
  vector<double> kappa_ci;
};

ostream& operator<< (ostream& os, const KappaResult& k) {
  os << fixed << setw(10) << setprecision(2) << k.alpha
  << setw(12) << setprecision(2) << k.mad_rel_err*100 << "%";
  for (size_t i=0; i<k.kappa_mad.size(); ++i)
    os << setw(13) << setprecision(3) << k.kappa_mad.at(i);
  os << setw(12) << setprecision(2) << k.ci_rel_err*100 << "%";
  for (size_t i=0; i<k.kappa_ci.size(); ++i)
    os << setw(13) << setprecision(3) << k.kappa_ci.at(i);
  os << endl;
  return os;
}

struct KappaResults {
  KappaResults(const vector<int>& ns, size_t n_alphas)
     : ns(ns), kr(n_alphas) {}
  KappaResult& at(size_t i) {
    unique_lock<mutex> lock(kr_mutex);
    return this->kr.at(i);
  }
  vector<int> ns;
  vector<KappaResult> kr;
  mutex kr_mutex;
};

ostream& operator<< (ostream& os, KappaResults& ks) {
  os << setw(10) << right << "alpha"
     << setw(13) << right << "mad_rel_err";
  for (auto n : ks.ns)
    os << setw(13) << right << "kappa"+to_string(n)+"_mad";
  os << setw(13) << right << "ci_rel_err";
  for (auto n : ks.ns)
    os << setw(13) << right << "kappa"+to_string(n)+"_ci";
  os << endl << endl;
  for (auto kr : ks.kr)
    os << kr;
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

double confidence_interval(int nmin, int nmax,
                           const vector<double>& p, const vector<double>& x,
                           double ci_level) {
  int nsize = static_cast<int>(p.size());
  double ci_lower;
  double p_lower=0;
  for (int j=nmin-1; j<=nmax && p_lower<ci_level/2; ++j) {
    ci_lower = x.at(mod(j, nsize));
    p_lower += p.at(mod(j, nsize));
  }
  double ci_upper;
  double p_upper = 0;
  for (int j=nmax; j>=nmin-1 && p_upper<ci_level/2; --j) {
    ci_upper = x.at(mod(j, nsize));
    p_upper += p.at(mod(j, nsize));
  }
  return ci_upper-ci_lower;
}

template<typename Dist>
void calculate_kappa(double delta,
                     double delta2,
                     vector<int> ns, Dist dist,
                     double ci_level,
                     KappaResult& k,
                     bool verbose = false) {
  int ns_max = *max_element(ns.begin(), ns.end());
  double mean = dist.mean();
  double alpha = dist.alpha();
  double dist_low = dist.quantile(delta2)-mean;
  double dist_high = dist.quantile(1-delta2)-mean;
  if (verbose)
    cout << "alpha = " << alpha << endl
         << "mean  = " << mean << endl
         << "dist_low = " << dist_low << endl
         << "dist_high = " << dist_high << endl;
  double xmin = min(0.,double(dist.quantile(.001)));
  xmin = std::min(xmin,mean*ns_max+pow(ns_max,max(.5,1/alpha))*dist_low);
  double xmax = max(0.,double(dist.quantile(.999)));
  xmax = std::max(xmax,mean*ns_max+pow(ns_max,max(.5,1/alpha))*dist_high);
  int nmax = static_cast<int>(xmax/delta);
  int nmin = static_cast<int>(xmin/delta);
  xmax = nmax * delta;
  xmin = nmin * delta;
  int nsize = ns_max*(nmax-nmin+2);
  if (verbose) {
    cout << "xmin = " << xmin << endl
         << "xmax = " << xmax << endl
         << "nmin = " << nmin << endl
         << "nmax = " << nmax << endl
         << "ns_max = " << ns_max << endl
         << "nsize = " << nsize << endl << endl;
  }
  vector<double> p(nsize,0.);
  vector<double> x(nsize,0.);
  double cdf_old = dist.cdf(xmin);
  p.at(mod(nmin-1, nsize)) = cdf_old;
  // The condiction tail expectation of x assuming x ~ -pareto distribution
  x.at(mod(nmin-1, nsize)) = xmin * alpha/((alpha-1)) + ((xmin>0)-(xmin<0))/(alpha-1);
  for (int i = nmin; i<nmax; ++i) {
    x.at(mod(i,nsize)) = i*delta+delta/2;
    double cdf_new = dist.cdf((i+1)*delta);
    p.at(mod(i,nsize)) = cdf_new-cdf_old;
    cdf_old = cdf_new;
  }
  p.at(mod(nmax, nsize)) = 1-cdf_old;
  // The conditional tail expectation of x assuming x has Pareto tail
  x.at(mod(nmax,nsize)) = xmax * alpha/(alpha-1) + ((xmax>0)-(xmax<0))/(alpha-1);
  double mad0 = 0;
  for (int j=nmin-1; j<=nmax; ++j) {
    mad0 += p.at(mod(j,nsize)) * fabs(x.at(mod(j,nsize))-mean);
  }
  k.mad_rel_err = rel_err(mad0, dist.mad());
  k.ci_rel_err = rel_err(confidence_interval(nmin, nmax, p, x, ci_level),
                         dist.ci(ci_level));

  vector<dcomplex> fft_p;
  fft_eng.fwd(fft_p,p);
  
  if (verbose) {
    cout << setw(5) << " "
         << setw(15) << right << "x.at(nmin-1)"
         << setw(15) << right << "x.at(-1)"
         << setw(15) << right << "x.at(0)"
         << setw(15) << right << "x.at(nmax)"
         << endl;
    cout << setw(5) << " "
         << setw(15) << fixed << setprecision(8) << x.at(mod(nmin-1,nsize))
         << setw(15) << fixed << setprecision(8) << x.at(mod(-1,nsize))
         << setw(15) << fixed << setprecision(8) << x.at(mod(0,nsize))
         << setw(15) << fixed << setprecision(8) << x.at(mod(nmax,nsize))
         << endl << endl;

    cout << setw(5) << right << "n"
         << setw(15) << right << "pn.at(nmin-1)"
         << setw(15) << right << "pn.at(-1)"
         << setw(15) << right << "pn.at(0)"
         << setw(15) << right << "pn.at(nmax)"
         << setw(15) << right << "pn_total"
         << setw(15) << right << "madn"
         << endl;
    double p_total = accumulate(p.begin(),p.end(),0.)-accumulate(p.begin()+mod(nmax+1,nsize),p.begin()+mod(nmin-1,nsize),0.);
    cout << setw(5) << 1
         << setw(15) << fixed << setprecision(8) << p.at(mod(nmin-1,nsize))
         << setw(15) << fixed << setprecision(8) << p.at(mod(-1,nsize))
         << setw(15) << fixed << setprecision(8) << p.at(mod(0,nsize))
         << setw(15) << fixed << setprecision(8) << p.at(mod(nmax,nsize))
         << setw(15) << fixed << setprecision(8) << p_total
         << setw(15) << fixed << setprecision(8) << mad0
         << endl;
  }
  for (auto n : ns) {
    vector<dcomplex> fft_pn(nsize);
    for (int i=0; i<nsize; ++i)
      fft_pn.at(i) = pow(fft_p.at(i),n);
    vector<double> pn;
    fft_eng.inv(pn, fft_pn);
    for (int j = nmax+1; j<=n*nmax; ++j) {
      pn.at(mod(nmax,nsize)) += pn.at(mod(j, nsize));
      pn.at(mod(j, nsize)) = 0;
    }
    for (int j = nmin-2; j>= n*(nmin-1); --j) {
      pn.at(mod(nmin-1, nsize)) += pn.at(mod(j,nsize));
      pn.at(mod(j, nsize)) = 0.;
    }
    if (verbose) {
      double p_total = accumulate(pn.begin(),pn.end(),0.)-accumulate(pn.begin()+mod(nmax+1,nsize),pn.begin()+mod(nmin-1,nsize),0.);
      cout << setw(5) << n
           << setw(15) << fixed << setprecision(8) << pn.at(mod(nmin-1, nsize))
           << setw(15) << fixed << setprecision(8) << pn.at(mod(-1, nsize))
           << setw(15) << fixed << setprecision(8) << pn.at(mod(0, nsize))
           << setw(15) << fixed << setprecision(8) << pn.at(mod(nmax,nsize))
           << setw(15) << fixed << setprecision(8) << p_total;
    }
    double madn=0;
    for (int j=nmin-1; j<=nmax; ++j) {
      madn += pn.at(mod(j,nsize)) * fabs(x.at(mod(j,nsize))-n*mean);
    }
    if (verbose)
      cout << setw(15) << fixed << setprecision(8) << madn << endl;
    k.kappa_mad.push_back(2- (log(n) - log(1))/(log(madn) - log(dist.mad())));
    double ci_n = confidence_interval(nmin, nmax, pn, x, ci_level);
    k.kappa_ci.push_back(2-(log(n)-log(1.))/(log(ci_n)-log(dist.ci(ci_level))));
  } // for n
}

template<typename Dist>
void calc_kappa(int m, int n0, vector<int> ns, Dist dist,
                          double ci_level,
                          KappaResult& k,
                          bool verbose = false) {
  k.kappa_mad.clear();
  k.kappa_ci.clear();
  mt19937 eng(5489u);
  double sum_at_n0=0;
  vector<double> p{ci_level/2,1-ci_level/2};
  vector<double> x;
  for (int i=0; i<m; ++i) {
    double sum = 0;
    for (int j=0; j<n0; ++j)
      sum += dist(eng);
    x.push_back(sum);
    sum_at_n0+=sum;
  }
  double mean_at_n0 = sum_at_n0/m;
  double mad_at_n0=0;
  for (int i=0; i<m; ++i) {
    mad_at_n0 += fabs(x.at(i)-mean_at_n0);
  }
  mad_at_n0 /= m;
  k.mad_rel_err = rel_err(mad_at_n0, dist.mad());
  vector<double> tmp = quantile(x,p);
  double ci_at_n0 = tmp.at(1)-tmp.at(0);
  k.ci_rel_err = rel_err(ci_at_n0, dist.ci(ci_level));
  if (verbose)
    cout << "n0 = " << setw(3)<< n0 << ", mad = "<< mad_at_n0 << ", ci = " << ci_at_n0 << endl;
  
  int n_old = n0;
  double sum_at_n = sum_at_n0;
  for (auto n : ns) {
    for (int i=0; i<m; ++i) {
      double sum = 0;
      for (int j=n_old; j<n; ++j)
        sum += dist(eng);
      x.at(i) += sum;
      sum_at_n += sum;
    }
    double mean_at_n = sum_at_n/m;
    double mad_at_n=0;
    for (int i=0; i<m; ++i) {
      mad_at_n += fabs(x.at(i)-mean_at_n);
    }
    mad_at_n /= m;
    tmp = quantile(x,p);
    double ci_at_n = tmp.at(1)-tmp.at(0);
    if (verbose)
      cout << "n  = " << setw(3)<< n << ", mad = "<< mad_at_n << ", ci = " << ci_at_n;
    k.kappa_mad.push_back(2-(log(n)-log(n0))/(log(mad_at_n)-log(dist.mad())));
    k.kappa_ci.push_back(2-(log(n)-log(n0))/(log(ci_at_n)-log(dist.ci(ci_level))));
    if (verbose)
      cout << ", kappa_mad = " << k.kappa_mad.back() << ", kappa_ci = " << k.kappa_ci.back() << endl;
    n_old = n;
  }
  return;
}

void show_usage(path p) {
  cerr << "Usage: " << p.filename().string() << " delta [delta2=.001]" << endl
       << "  where delta is the interval of integration" << endl
       << "  and delta2 is the confidence level for the endpoints" << endl;
}


int main(int argc, const char * argv[]) {
  path p(argv[0]);
  if ( (argc < 2) || (argc > 3) ) {
    show_usage(p);
    return 1;
  }
  istringstream iss1{string(argv[1])};
  double delta;
  if ( !(iss1 >> delta) || ( delta <=0 ))
    show_usage(p);
  
  double delta2 = .001;
  if (argc == 3) {
    istringstream iss2{string(argv[2])};
    if ( !(iss2 >> delta2) || (delta2 <=0 ) )
      show_usage(p);
  }
  
  vector<double> alphas;
  
  for (double alpha=1.25; alpha<=4; alpha+=.25)
    alphas.push_back(alpha);
  
  vector<int> ns{2};
  KappaResults ks_pareto(ns, alphas.size());
  KappaResults ks_student(ns, alphas.size());
  
  double ci_level = .05;
  
  string out_dir{"../output"};
  ostringstream oss;
  oss << out_dir << "/pareto_test_"
      << delta << "_"
      << delta2
      << ".out";
  cout << "Writing to file " << oss.str() << endl << endl;
  cout.flush();
  
  ofstream out{oss.str()};
  auto_cpu_timer t(out);
  for (size_t i=0; i<alphas.size(); ++i) {
    double alpha=alphas.at(i);
    pareto_distribution<> pd(alpha);
    KappaResult kr;
    kr.alpha = alpha;
    calculate_kappa(delta, delta2, ns, pd, ci_level, kr, true);
    ks_pareto.at(i) = kr;
  }
  out << "Pareto Distribution" << endl << endl;
  out << ks_pareto;
  out.flush();
  
  for (size_t i=0; i<alphas.size(); ++i) {
    double alpha = alphas.at(i);
    student_t_distribution<> td(alpha);
    KappaResult kr;
    kr.alpha = alpha;
    calculate_kappa(delta, delta2, ns, td, ci_level, kr, true);
    ks_student.at(i)=kr;
  }
  out << "Student Distribution" << endl << endl;
  out << ks_student;

    return 0;
}
