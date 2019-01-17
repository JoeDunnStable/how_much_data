//
//  taleb_results.h
//  pareto_test
//
//  Created by Joseph Dunn on 1/12/19.
//  Copyright © 2019 Joseph Dunn. All rights reserved.
//

#ifndef taleb_results_h
#define taleb_results_h

// From
// How Much Data Do You Need?
// A Pre-asymptotic Metric for Fat-tailedness
// Nassim Nicholas Taleb
// Tandon School of Engineering, New York University
// November 2018
// Forthcoming, International Journal of Forecasting

//               Table III
//  COMPARING PARETO TO STUDENT T (SAME TAIL EXPONENT )
//        Pareto  Pareto  Pareto  Student  Student  Student
// alpha     1     1;30   1;100      1       1;30    1;100

vector<vector<double> > taleb_results ={
  {1.25,   0.829,  0.787,  0.771,  0.792,    0.765,  0.756},
  {1.5,    0.724,  0.65,   0.631,  0.647,    0.609,  0.587},
  {1.75,   0.65,   0.556,  0.53,   0.543,    0.483,  0.451},
  {2,      0.594,  0.484,  0.449,  0.465,    0.387,  0.352},
  {2.25,   0.551,  0.431,  0.388,  0.406,    0.316,  0.282},
  {2.5,    0.517,  0.386,  0.341,  0.359,    0.256,  0.227},
  {2.75,   0.488,  0.356,  0.307,  0.321,    0.224,  0.189},
  {3,      0.465,  0.3246, 0.281,  0.29,     0.191,  0.159},
  {3.25,   0.445,  0.305,  0.258,  0.265,    0.167,  0.138},
  {3.5,    0.428,  0.284,  0.235,  0.243,    0.149,  0.121},
  {3.75,   0.413,  0.263,  0.222,  0.225,    0.13,   0.1},
  {4,      0.4,    0.2532, 0.211,  0.209,    0.126,  0.093}};

#endif /* taleb_results_h */
