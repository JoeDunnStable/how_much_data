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
  
  vector<double> sigmas = {.01, .1, 1, 5, 10};
  vector<double> omegas;
  for (int i=64; i>=-64; --i)
    omegas.push_back(-pow(2,i/8.));
  omegas.push_back(0.);
  for (int i=-64; i<=64; i++)
    omegas.push_back(pow(2,i/8.));

 {
    string ofile_name{"../output/lognormal_test_lambert_w.out"};
    cout << "Writing to file " << ofile_name << endl;
    ofstream out(ofile_name);
    auto_cpu_timer timer;
    out << setw(11) << "sigma,"
    << setw(11) << "omega,"
    << setw(17) << "cf_Re,"
    << setw(17) << "cf_Im,"
    << setw(17) << "cf_Abs"
    << endl;
    
    for (auto sigma : sigmas) {
      lognormal_distribution<> lnd(0, sigma, 3);
      for (auto omega : omegas ) {
        complex<double> cf = lnd.characteristic_function(omega);
        out << setw(10) << setprecision(2) << sigma << ","
        << setw(10) << setprecision(4) << omega << ","
        << setw(16) << setprecision(8) << real(cf) << ","
        << setw(16) <<setprecision(8) << imag(cf) << ","
        << setw(16) << setprecision(8) << abs(cf)
        << endl;
      }
    }
    out << endl;
  }
  {
    string ofile_name{"../output/lognormal_test_fourier_x.out"};
    cout << "Writing to file " << ofile_name << endl;
    ofstream out(ofile_name);
    auto_cpu_timer timer;
    out << setw(11) << "sigma,"
    << setw(11) << "omega,"
    << setw(17) << "cf_Re,"
    << setw(17) << "cf_Im,"
    << setw(17) << "cf_Abs"
    << endl;
    
    for (auto sigma : sigmas) {
      lognormal_distribution<> lnd(0, sigma, 2);
      for (auto omega : omegas ) {
        complex<double> cf = lnd.characteristic_function(omega);
        out << setw(10) << setprecision(2) << sigma << ","
        << setw(10) << setprecision(4) << omega << ","
        << setw(16) << setprecision(8) << real(cf) << ","
        << setw(16) <<setprecision(8) << imag(cf) << ","
        << setw(16) << setprecision(8) << abs(cf)
        << endl;
      }
    }
    out << endl;
  }

  {
    string ofile_name{"../output/lognormal_test_fourier_lnx.out"};
    cout << "Writing to file " << ofile_name << endl;
    ofstream out(ofile_name);
    auto_cpu_timer timer;
    out << setw(11) << "sigma,"
    << setw(11) << "omega,"
    << setw(17) << "cf_Re,"
    << setw(17) << "cf_Im,"
    << setw(17) << "cf_Abs"
    << endl;
    
    for (auto sigma : sigmas) {
      lognormal_distribution<> lnd(0, sigma, 1);
      for (auto omega : omegas ) {
        complex<double> cf = lnd.characteristic_function(omega);
        out << setw(10) << setprecision(2) << sigma << ","
        << setw(10) << setprecision(4) << omega << ","
        << setw(16) << setprecision(8) << real(cf) << ","
        << setw(16) <<setprecision(8) << imag(cf) << ","
        << setw(16) << setprecision(8) << abs(cf)
        << endl;
      }
    }
    out << endl;
  }
 
  {
    string ofile_name{"../output/lognormal_test_spline.out"};
    cout << "Writing to file " << ofile_name << endl;
    ofstream out(ofile_name);
    auto_cpu_timer timer;
    out << setw(11) << "sigma,"
    << setw(11) << "omega,"
    << setw(17) << "cf_Re,"
    << setw(17) << "cf_Im,"
    << setw(17) << "cf_Abs"
    << endl;
    
    for (auto sigma : sigmas) {
      cpu_timer timer;
      lognormal_distribution<> lnd(0, sigma, 4);
      cout << format(timer.elapsed()) << endl;
      for (auto omega : omegas ) {
        complex<double> cf = lnd.characteristic_function(omega);
        out << setw(10) << setprecision(2) << sigma << ","
        << setw(10) << setprecision(4) << omega << ","
        << setw(16) << setprecision(8) << real(cf) << ","
        << setw(16) <<setprecision(8) << imag(cf) << ","
        << setw(16) << setprecision(8) << abs(cf)
        << endl;
      }
    }
    out << endl;
  }

}
