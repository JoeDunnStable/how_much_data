//
//  monte_carlo_test.cpp
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
  for (size_t j=1; j<ks.ns.size(); ++j)
    os << setw(13) << right << "kappa"+to_string(ks.ns.at(j))+"_mad";
  os << setw(13) << right << "ci_rel_err";
  for (size_t j=1; j<ks.ns.size(); ++j)
    os << setw(13) << right << "kappa"+to_string(ks.ns.at(j))+"_ci";
  os << endl << endl;
  for (auto kr : ks.kr)
    os << kr;
  os << endl;
  return os;
}

template<typename Dist>
void calc_kappa(int m, vector<int> ns, Dist dist,
                double ci_level,
                KappaResult& k,
                bool verbose = false) {
  k.kappa_mad.clear();
  k.kappa_ci.clear();
  mt19937 eng(5489u);
  int n0 = ns.at(0);
  vector<double> p{ci_level/2,1-ci_level/2};
  vector<double> mean_at_n(ns.size(),0.);
  vector<double> mad_at_n(ns.size(),0.);
  vector<double> ci_at_n(ns.size(),0.);
  int m_ci_limit = min(10000000,m);
  double mean = dist.mean();
  // Up to m_ci_limit we'll save the individual run results calculate
  // confidence intervals
  vector<vector<double> > x_ci(ns.size());
  /*
  k.mad_rel_err = rel_err(mad_at_n0, dist.mad());
  vector<double> tmp = quantile(x,p);
  double ci_at_n0 = tmp.at(1)-tmp.at(0);
  k.ci_rel_err = rel_err(ci_at_n0, dist.ci(ci_level));
  if (verbose)
    cout << "n0 = " << setw(3)<< n0 << ", mad = "<< mad_at_n0 << ", ci = " << ci_at_n0 << endl;
  */
  for (int i=0; i<m; ++i) {
    int n_old = 0;
    double sum = 0;
    for (size_t j=0; j<ns.size(); j++) {
      int n = ns.at(j);
      for (int jj=n_old; jj<n; ++jj)
        sum += dist(eng);
      if (i < m_ci_limit)
        x_ci.at(j).push_back(sum);
      mean_at_n.at(j) += sum;
      mad_at_n.at(j) += fabs(sum-n*mean);
      n_old = n;
    }
  } // i
  for (size_t j=0; j<ns.size(); ++j) {
    mean_at_n.at(j) /= m;
    mad_at_n.at(j) /= m;
    vector<double> tmp = quantile(x_ci.at(j),p);
    ci_at_n.at(j) = tmp.at(1)-tmp.at(0);
    if (verbose)
      cout << "n  = " << setw(3)<< ns.at(j)
      << ", mad = "<< mad_at_n.at(j)
      << ", ci = " << ci_at_n.at(j);
    if (j > 0) {
      int n = ns.at(j);
      k.kappa_mad.push_back(2-(log(n)-log(n0))/(log(mad_at_n.at(j))-log(mad_at_n.at(0))));
      k.kappa_ci.push_back(2-(log(n)-log(n0))/(log(ci_at_n.at(j))-log(ci_at_n.at(0))));
      if (verbose)
        cout << ", kappa_mad = " << k.kappa_mad.back() << ", kappa_ci = " << k.kappa_ci.back() << endl;
    }
  } // k
  k.mad_rel_err = rel_err(mad_at_n.at(0), dist.mad());
  k.ci_rel_err = rel_err(ci_at_n.at(0), dist.ci(ci_level));
  return;
}

void show_usage(path p) {
  cerr << "Usage: " << p.filename().string() << " m" << endl
  << "  where m is number of trials" << endl;
}


int main(int argc, const char * argv[]) {
  path p(argv[0]);
  if ( argc != 2) {
    show_usage(p);
    return 1;
  }
  int m;
  istringstream iss1{argv[1]};
  if (!(iss1 >> m) || m < 1) {
    show_usage(p);
    return 1;
  }

  vector<double> alphas;
  
  for (double alpha=1.25; alpha<=4; alpha+=.25)
    alphas.push_back(alpha);
  
  vector<int> ns{1, 2};
  KappaResults ks_pareto(ns, alphas.size());
  KappaResults ks_student(ns, alphas.size());
  
  double ci_level = .05;
  
  string out_dir{"../output"};
  ostringstream oss;
  oss << out_dir << "/monte_carlo_test_"  << m << ".out";
  cout << "Writing to file " << oss.str() << endl << endl;
  cout.flush();
  
  ofstream out{oss.str()};
  auto_cpu_timer t(out);
  out << "m =     " << m << endl
  << "ns =     ";
  for (auto n : ns)
    out << n << " ";
  out << endl<< endl;
  out.flush();
  
#pragma omp parallel for
  for (size_t i=0; i<alphas.size(); ++i) {
    double alpha=alphas.at(i);
    pareto_distribution<> pd(alpha);
    KappaResult kr;
    kr.alpha = alpha;
    calc_kappa(m, ns, pd, ci_level, kr, false);
    ks_pareto.at(i) = kr;
  }
  out << "Pareto Distribution" << endl << endl;
  out << ks_pareto;
  out.flush();
  
#pragma omp parallel for
  for (size_t i=0; i<alphas.size(); ++i) {
    double alpha = alphas.at(i);
    student_t_distribution<> td(alpha);
    KappaResult kr;
    kr.alpha = alpha;
    calc_kappa(m, ns, td, ci_level, kr, false);
    ks_student.at(i)=kr;
  }
  out << "Student Distribution" << endl << endl;
  out << ks_student;
  return 0;
}
