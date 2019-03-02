doxygen = find_program('doxygen')

if(doxygen.found())
  src_doxygen = [meson.source_root()+'/README.md', \
                 meson.source_root()+'/include/pareto_distribution.h', \
                 meson.source_root()+'/include/student_t_distribution.h', \
                 meson.source_root()+'/include/exponential_distribution.h', \
                 meson.source_root()+'/include/lognormal_distribution.h', \
                 meson.source_root()+'/monte_carlo_test/monte_carlo_test.cpp', \
		 meson.source_root()+'/pinelis_taleb_test/pinelis_taleb_test.cpp']

	doc_config = configuration_data()
	doc_config.set('PACKAGE_NAME', meson.project_name())
	doc_config.set('PACKAGE_VERSION', meson.project_version())
	doc_config.set('PROJECT_SOURCE_DIR', meson.source_root())
	doc_config.set('INPUT', ' '.join(src_doxygen))

	doxyfile = configure_file(input : 'Doxyfile.in',
				  output : 'Doxyfile',
				  configuration : doc_config,
				  install : false)
	custom_target('doxygen',
		      input : doxyfile,
		      output : ['html', 'latex']  ,
		      command : [ doxygen, doxyfile ],
		      depend_files: src_doxygen,
		      install : false,
		      build_by_default : true)

endif
