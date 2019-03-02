src = ['pareto_distribution.h',
       'student_t_distribution.h',
       'exponential_distribution.h',
       'lognormal_distribution.h',
       'lambert_w.h',
       'convolution_test.cpp']

convolution_test = executable('convolution_test', 
                         src,
                         dependencies : alldep,
                         include_directories : incdir)

test('convolution_test', convolution_test, args: '.025 .0001', timeout : 3600)
