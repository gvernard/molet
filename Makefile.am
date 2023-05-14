ACLOCAL_AMFLAGS=-I m4

AM_CXXFLAGS = -std=c++11 -fPIC -g -frounding-math

BUILT_SOURCES = libinstruments.la


lib_LTLIBRARIES = libinstruments.la

libinstruments_la_SOURCES= \
	instruments/src/instruments.cpp \
	instruments/src/noise.cpp
libinstruments_la_CXXFLAGS = -std=c++17	-fPIC -g -frounding-math
libinstruments_la_CPPFLAGS = -Iinstruments/include -DINSTRUMENT_PATH='"@INSTRUMENT_PATH@"'
libinstruments_la_LIBADD = -lgmp -lCGAL -lvkl -lfftw3 -ljsoncpp




bin_PROGRAMS = \
	aux_create_instrument \
	angular_diameter_distances \
	check_time_vector \
	check_initialization \
	check_get_map_path \
	check_confirm_convolutions \
	check_confirm_convolutions_vs \
	check_coolest_conversion \
	vkl_fproject \
	vkl_point_source \
	vkl_lens \
	var_moving_fixed_source \
	var_expanding_source \
	var_moving_variable_source \
	var_match_to_gerlumph \
	var_drw \
	final_combine_light


aux_create_instrument_SOURCES = instruments/src/create_new_instrument.cpp
aux_create_instrument_CXXFLAGS = -std=c++17 -fPIC -g -frounding-math
aux_create_instrument_CPPFLAGS = -Iinstruments/include
aux_create_instrument_LDADD = -linstruments -lcfitsio -lCCfits -ljsoncpp

angular_diameter_distances_SOURCES = \
	cosmology/angular_diameter_distances/src/angular_diameter_distances.cpp \
	cosmology/angular_diameter_distances/src/auxiliary_functions.cpp
angular_diameter_distances_CPPFLAGS = -Icosmology/angular_diameter_distances/include
angular_diameter_distances_LDADD = -ljsoncpp

check_time_vector_SOURCES = checks/src/time_vector_checks.cpp
check_time_vector_CXXFLAGS = -std=c++17
check_time_vector_CPPFLAGS = -Iinstruments/include
check_time_vector_LDADD = -linstruments -lcfitsio -lCCfits -lgerlumph -lpng -ljsoncpp

check_initialization_SOURCES = checks/src/initialization_checks.cpp
check_initialization_CPPFLAGS = -Iinstruments/include
check_initialization_LDADD = -linstruments -ljsoncpp

check_get_map_path_SOURCES = checks/src/get_map_path.cpp
check_get_map_path_LDADD = -lgerlumph

check_confirm_convolutions_SOURCES = checks/src/confirm_convolutions.cpp
check_confirm_convolutions_LDADD = -lgerlumph -ljsoncpp

check_confirm_convolutions_vs_SOURCES = checks/src/confirm_convolutions_vs.cpp
check_confirm_convolutions_vs_LDADD = -lgerlumph -ljsoncpp

check_coolest_conversion_SOURCES = checks/src/coolest_conversion.cpp
check_coolest_conversion_CXXFLAGS = -std=c++17 -fPIC -g -frounding-math
check_coolest_conversion_CPPFLAGS = -Iinstruments/include
check_coolest_conversion_LDADD = -linstruments -ljsoncpp

vkl_fproject_SOURCES = \
	lensed_extended_source/vkl_fproject/src/fproject.cpp \
	lensed_extended_source/vkl_fproject/src/caustics.cpp
vkl_fproject_CPPFLAGS = -Iinstruments/include -Ilensed_extended_source/vkl_fproject/include
vkl_fproject_LDADD = -linstruments -lgmp -lCGAL -lvkl -lgfortran -lCCfits -lcfitsio -ljsoncpp

vkl_point_source_SOURCES = \
	lensed_point_source/vkl_point_source/src/point_source.cpp \
	lensed_point_source/vkl_point_source/src/polygons.cpp
vkl_point_source_CPPFLAGS = -Iinstruments/include -Ilensed_point_source/vkl_point_source/include
vkl_point_source_LDADD = -linstruments -lgmp -lCGAL -lvkl -lgfortran -lCCfits -lcfitsio -ljsoncpp

vkl_lens_SOURCES = lens_light_mass/vkl_llm/src/lens_light_mass.cpp
vkl_lens_CPPFLAGS = -Iinstruments/include
vkl_lens_LDADD = -linstruments -lgmp -lCGAL -lvkl -lgfortran -lCCfits -lcfitsio -ljsoncpp

var_moving_fixed_source_SOURCES = \
	variability/extrinsic/moving_fixed_source/src/moving_fixed_source.cpp \
	variability/extrinsic/moving_fixed_source/src/auxiliary_functions.cpp
var_moving_fixed_source_CPPFLAGS = -Iinstruments/include -Ivariability/extrinsic/moving_fixed_source/include
var_moving_fixed_source_LDADD = -linstruments -lgerlumph -lCCfits -lcfitsio -lpng -ljsoncpp

var_expanding_source_SOURCES = variability/extrinsic/expanding_source/src/expanding_source.cpp
var_expanding_source_CPPFLAGS = -Iinstruments/include
var_expanding_source_LDADD = -linstruments -lgerlumph -lCCfits -lcfitsio -lpng -ljsoncpp

var_moving_variable_source_SOURCES = variability/extrinsic/moving_variable_source/src/moving_variable_source.cpp
var_moving_variable_source_CPPFLAGS = -Iinstruments/include
var_moving_variable_source_LDADD = -linstruments -lgerlumph -lCCfits -lcfitsio -lpng -ljsoncpp

var_match_to_gerlumph_SOURCES = variability/extrinsic/match_to_gerlumph/src/match_to_gerlumph.cpp
var_match_to_gerlumph_CPPFLAGS = -Ivariability/extrinsic/match_to_gerlumph/include
var_match_to_gerlumph_LDADD = -lsqlite3 -ljsoncpp

var_drw_SOURCES = \
	variability/intrinsic/DRW/src/drw.cpp \
	variability/intrinsic/DRW/src/auxiliary_functions.cpp
var_drw_CPPFLAGS = -Iinstruments/include -Ivariability/intrinsic/DRW/include
var_drw_LDADD = -linstruments -lgerlumph -lfftw3 -lCCfits -lcfitsio -lpng -ljsoncpp

final_combine_light_SOURCES = \
	combine_light/src/combine_light.cpp \
	combine_light/src/auxiliary_functions.cpp
final_combine_light_CPPFLAGS = -Iinstruments/include -Icombine_light/include
final_combine_light_LDADD = -linstruments -lgmp -lCGAL -lvkl -lfftw3 -lCCfits -lcfitsio -ljsoncpp


bin_SCRIPTS = molet_driver
CLEANFILES = $(bin_SCRIPTS)


install-exec-hook:
	$(abs_top_builddir)/data/initialize_instruments.sh $(abs_top_builddir)
	rm $(bin_PROGRAMS)

uninstall-hook:
	rm -r $(abs_top_builddir)/instrument_modules/*
