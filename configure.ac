AC_INIT([molet], [0.1], [gvernardos@amnh.org])
AC_LANG(C++)
AC_CONFIG_AUX_DIR(aux-dist)
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_MACRO_DIRS([m4])
AC_PREFIX_DEFAULT($PWD)

AC_SUBST(MOLET_HOME,$PWD)
AC_SUBST(INSTRUMENT_PATH,$PWD/instrument_modules/)

AC_PROG_CXX
AC_PROG_MKDIR_P # MKDIR_P
AC_PROG_SED # SED

LT_INIT([disable-static])



AC_DEFUN([AC_INSTALLATION_INSTRUCTIONS],
	[
		name=$1
		echo ""
		echo "Program/library $(tput bold)$(tput setaf 1)${name} is missing$(tput sgr0), you need to install it using your package manager."
		echo "Alternatively, you can install ${name} from source using the following instructions:"
		echo $(tput setaf 6)
		cat ${MOLET_HOME}/third_party/install_${name}.txt
		echo $(tput sgr0)
		echo "If you use the 'prefix' option while running ./configure above, then you need to pass:"
		echo "   $(tput setaf 2)--with-${name}=prefix$(tput sgr0)"
		echo "when running the MOLET ./configure script."
		AC_MSG_ERROR("Program/library ${name} not found")
	]
)




AC_ARG_WITH([jq],
    [AS_HELP_STRING([--with-jq=prefix],[give the path to the jq program binary file.])],
    [
	AC_SUBST(AM_JQ,${withval}/bin/jq)
    ],
    [
	AC_PATH_PROG(AM_JQ,jq,[no],[path=$PATH])
	if test "$AM_JQ" = no ; then
	   AC_INSTALLATION_INSTRUCTIONS(jq)
	fi
    ])



ac_new_CPPFLAGS=""
ac_new_LDFLAGS=""
ac_new_LIBS=""

AC_ARG_WITH([fftw3],
    [AS_HELP_STRING([--with-fftw3=prefix],[give the path to the fftw3 library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lfftw3"
    ],
    [
	AC_CHECK_LIB(fftw3,main,[],[AC_INSTALLATION_INSTRUCTIONS(fftw3)])
	AC_CHECK_HEADERS(fftw3.h,[],[AC_MSG_ERROR("you need to have the fftw3 header files installed.")])
    ])

AC_ARG_WITH([cfitsio],
    [AS_HELP_STRING([--with-cfitsio=prefix],[give the path to the cfitsio library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lcfitsio"
    ],
    [
	AC_CHECK_LIB(cfitsio,main,[],[AC_INSTALLATION_INSTRUCTIONS(cfitsio)])
	AC_CHECK_HEADERS(fitsio.h,[],[AC_MSG_ERROR("you need to have the cfitsio header files installed.")])
    ])

AC_ARG_WITH([CCfits],
    [AS_HELP_STRING([--with-CCfits=prefix],[give the path to the CCfits library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lCCfits"
    ],
    [
	AC_CHECK_LIB(CCfits,main,[],[AC_INSTALLATION_INSTRUCTIONS(CCfits)])
	AC_CHECK_HEADERS(CCfits/CCfits.h,[],[AC_MSG_ERROR("you need to have the CCfits header files installed.")])
    ])

AC_ARG_WITH([jsoncpp],
    [AS_HELP_STRING([--with-jsoncpp=prefix],[give the path to the jsoncpp library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -ljsoncpp"
    ],
    [
	AC_CHECK_LIB(jsoncpp,main,[],[AC_INSTALLATION_INSTRUCTIONS(jsoncpp)])
	AC_CHECK_HEADERS(json/json.h,[],[AC_MSG_ERROR("you need to have the jsoncpp header files installed.")])
    ])

AC_ARG_WITH([png],
    [AS_HELP_STRING([--with-png=prefix],[give the path to the png library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lpng"
    ],
    [
	AC_CHECK_LIB(png,main,[],[AC_INSTALLATION_INSTRUCTIONS(png)])
	AC_CHECK_HEADERS(png.h,[],[AC_MSG_ERROR("you need to have the png header files installed.")])
    ])

AC_ARG_WITH([sqlite3],
    [AS_HELP_STRING([--with-sqlite3=prefix],[give the path to the sqlite3 library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lsqlite3"
    ],
    [
	AC_CHECK_LIB(sqlite3,main,[],[AC_INSTALLATION_INSTRUCTIONS(sqlite3)])
	AC_CHECK_HEADERS(sqlite3.h,[],[AC_MSG_ERROR("you need to have the sqlite3 header files installed.")])
    ])

AC_ARG_WITH([vkl],
    [AS_HELP_STRING([--with-vkl=prefix],[give the path to the vkl library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lvkl"
    ],
    [
	AC_CHECK_LIB(vkl,main,[],[AC_INSTALLATION_INSTRUCTIONS(vkl)])
	AC_CHECK_HEADERS(vkllib.hpp,[],[AC_MSG_ERROR("you need to have the vkl header files installed.")])
    ])

AC_ARG_WITH([gmp],
    [AS_HELP_STRING([--with-gmp=prefix],[give the path to the gmp library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lgmp"
    ],
    [
	AC_CHECK_LIB(gmp,main,[],[AC_MSG_ERROR("you need to have the gmp library installed.")])
	AC_CHECK_HEADERS(gmp.h,[],[AC_MSG_ERROR("you need to have the gmp header files installed.")])
    ])

AC_ARG_WITH([CGAL],
    [AS_HELP_STRING([--with-CGAL=prefix],[give the path to the CGAL library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lCGAL"
    ],
    [
	AC_CHECK_LIB(CGAL,main,[],[AC_MSG_ERROR("you need to have the CGAL library installed.")])
	#AC_CHECK_HEADERS(CGAL/CGAL.h,[],[AC_MSG_ERROR("you need to have the CGAL header files installed.")])
    ])

AC_ARG_WITH([gerlumph],
    [AS_HELP_STRING([--with-gerlumph=prefix],[give the path to the gerlumph library and header files.])],
    [
	ac_new_CPPFLAGS+=" -I${withval}/include"
	ac_new_LDFLAGS+=" -L${withval}/lib -Wl,-rpath -Wl,${withval}/lib"
	ac_new_LIBS+=" -lgerlumph"
    ],
    [
	AC_CHECK_LIB(gerlumph,main,[],[AC_INSTALLATION_INSTRUCTIONS(gerlumph)])
	AC_CHECK_HEADERS(gerlumph.hpp,[],[AC_MSG_ERROR("you need to have the gerlumph header files installed.")])
    ])
    
dnl echo $ac_new_CPPFLAGS
dnl echo $ac_new_LDFLAGS
dnl echo $ac_new_LIBS

CPPFLAGS+="$ac_new_CPPFLAGS"
LDFLAGS+="$ac_new_LDFLAGS"
LIBS+="$ac_new_LIBS"

dnl echo $CPPFLAGS
dnl echo $LDFLAGS
dnl echo $LIBS

AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])
AC_SUBST([LIBS])


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([molet_driver],[chmod +x molet_driver])
AC_CONFIG_FILES([checks/check_map_files.sh],[chmod +x checks/check_map_files.sh])
AC_CONFIG_FILES([combine_variability/setup_dirs.sh],[chmod +x combine_variability/setup_dirs.sh])

AC_OUTPUT
