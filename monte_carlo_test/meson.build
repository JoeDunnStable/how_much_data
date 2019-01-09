src = ['monte_carlo_test.cpp']

monte_carlo_test = executable('monte_carlo_test', 
                         src,
                         dependencies : alldep,
                         include_directories : incdir)

test('monte_carlo_test', monte_carlo_test, args: '10000000', timeout : 3600)
