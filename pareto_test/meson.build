src = ['pareto_distribution.h',
       'student_t_distribution.h',
       'main.cpp']

pareto_test = executable('pareto_test', 
                         src,
                         dependencies : alldep,
                         include_directories : incdir)

test('pareto_test', pareto_test, args: '.025 .0001', timeout : 3600)
