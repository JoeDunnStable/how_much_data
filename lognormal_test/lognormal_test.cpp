//
/// \file lognormal_test.cpp
/// \package how_much_data
//
/// \author Created by Joseph Dunn on 1/25/19.
/// \copyright © 2019 Joseph Dunn. All rights reserved.
//

#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
using std::setw;
using std::setprecision;
using std::fixed;
#include <vector>
using std::vector;
#include <fstream>
using std::ofstream;
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;
using boost::timer::cpu_timer;
#include "lognormal_distribution.h"

int main(int argc, const char * argv[]) {
  Kronrod<double> k_big(10);
  int noext = 0;
  double epsabs_double = std::numeric_limits<double>::epsilon()/128;
  double epsrel_double = boost::math::tools::root_epsilon<double>();
  int limit = 1000;
  int verbose_integration = 0;
  IntegrationController<double> cf_ctl(noext, k_big,
                                       epsabs_double, epsrel_double,
                                       limit, verbose_integration);

  string ofile_name{"../output/lognormal_test_fourier_integrand.out"};
  cout << "Writing to file " << ofile_name << endl;
  ofstream trace(ofile_name);
  lognormal_distribution<>::print_fourier_integrand(trace, 1000,
                                            complex<double>(128,0.), .1);

  vector<double> sigmas = {.01, .1, 1, 5, 10};
  vector<complex<double> > omegas;
  omegas.push_back(0.);
  complex<double> i{0,1};
  for (int j=-512; j<=512; j++)
    omegas.push_back(pow(2,j/8.));
  
  vector<complex<double> > cfs_std;
  vector<complex<double> > cfprimes_std;
  
  vector<string> names = {"fourier_lnx", "fourier", "lambert_w","fourier_mixed", "adj_series"};
  vector<int> types ={2, 5, 3,4,1};
  
  for (auto type : types)
  {
    string ofile_name{"../output/lognormal_test_"+names.at(type-1)+".out"};
    cout << "Writing to file " << ofile_name << endl;
    ofstream out(ofile_name);
    auto_cpu_timer timer;
    out << setw(11) << "sigma,"
    << setw(14) << "omega,"
    << setw(17) << "cf_Re,"
    << setw(17) << "cf_Im,"
    << setw(17) << "Abs(cf-cf_std),"
    << setw(17) << "Integ_err,"
    << setw(11) << "# eval,"
    << setw(17) << "Abs(cf) -1"
    << setw(17) << "cfprime_Re,"
    << setw(17) << "cfprime_Im,"
    << setw(17) << "Abs(cfp-cfp_std),"
    << setw(17) << "Integ_err,"
    << setw(11) << "# eval"
    << endl;

    int j = 0;
    for (auto sigma : sigmas) {
      lognormal_distribution<> lnd(0, sigma, cf_ctl, type);
      for (auto omega : omegas ) {
        double cf_integ_err;
        double cf_l1_norm;
        int cf_neval;
        bool adjusted = (type==5 || type==6);
        double mean = adjusted ? 0 : lnd.mean();
        complex<double> fac = exp(- i * mean * omega);
        complex<double> cf = fac* lnd.characteristic_function(omega,
                                                              &cf_integ_err,
                                                              &cf_l1_norm,
                                                              &cf_neval
                                                              );
        bool std = type == types.at(0);
        complex<double> cf_std = std ? cf : cfs_std.at(j);
        if (std) cfs_std.push_back(cf);
        double cfp_integ_err;
        double cfp_l1_norm;
        int cfp_neval;
        complex<double> cfprime = fac*lnd.characteristic_function_prime(omega,
                                                                        &cfp_integ_err,
                                                                        &cfp_l1_norm,
                                                                        &cfp_neval
                                                                        )-i*mean*cf;
        complex<double> cfprime_std = std ? cfprime : cfprimes_std.at(j++);
        if (std) cfprimes_std.push_back(cfprime);
        out << setw(10) << setprecision(2) << sigma << ","
        << setw(13) << setprecision(4) << omega << ","
        << setw(16) << setprecision(8) << real(cf) << ","
        << setw(16) << setprecision(8) << imag(cf) << ","
        << setw(16) << setprecision(8) << abs(cf-cf_std) << ","
        << setw(16) << setprecision(8) << cf_integ_err << ","
        << setw(10) << cf_neval << ","
        << setw(16) << setprecision(8) << abs(cf) - 1. << ","
        << setw(16) << setprecision(8) << real(cfprime) << ","
        << setw(16) << setprecision(8) << imag(cfprime) << ","
        << setw(16) << setprecision(8) << abs(cfprime-cfprime_std) << ","
        << setw(16) << setprecision(8) << cfp_integ_err << ","
        << setw(10) << cfp_neval
        << endl;
      }
    }
    out << endl;
  }

}
