project('How Much Data', 'cpp', default_options : ['cpp_std=c++11'],
         version : '0.6.0')

compiler=meson.get_compiler('cpp').get_id()

if compiler == 'gcc'
  add_global_arguments(
      [ '-fopenmp'] ,
      language : 'cpp' )

endif

#Boost header and some libraries.
#I had to define BOOST_ROOT to get this to work.
boostdep = dependency('boost', modules : ['filesystem', 'timer', 'system'],  required : true)

gompdep = meson.get_compiler('cpp').find_library('gomp')

fftwdep = dependency('fftw3' , required : true)

#Eigen3 headers
eigendep = dependency('eigen3', version : '>3.0', required : true)

alldep = [boostdep, gompdep, fftwdep, eigendep]

if host_machine.system() == 'linux'
  alldep = [alldep, meson.get_compiler('cpp').find_library('pthread')]
endif


incdir=include_directories([ 'convolution_test'])

subdir('convolution_test')
subdir('monte_carlo_test')

add_test_setup('boost_gcc', env: 'DYLD_LIBRARY_PATH=/usr/local/boost_gcc/lib')
