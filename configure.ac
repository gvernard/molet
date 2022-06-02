AC_INIT([molet], [0.1], [gvernard@astro.rug.nl])
AC_CONFIG_AUX_DIR(aux-dist)
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_MACRO_DIRS([m4])
AC_PREFIX_DEFAULT($PWD)

AC_PROG_CXX

LT_INIT([disable-static])


dnl AC_ARG_WITH([instrument-path],
dnl     [AS_HELP_STRING([--with-instrument-path=prefix],[give the path to the instrument modules])],
dnl     [INSTRUMENT_PATHSET=1],
dnl     [INSTRUMENT_PATHSET=0])
	
dnl # if optional parameter used, extend path flags for compliler and linker
dnl if test $INSTRUMENT_PATHSET = 1 ; then
dnl     INSTRUMENT_PATH="$with_instrument_path"
dnl else
dnl     AC_MSG_ERROR([You need to provide the INSTRUMENT_PATH to comile the instruments library!])
dnl fi


AC_CHECK_LIB(fftw3,main,[],[
	echo "Error! You need to have fftw3 installed."
	exit -1
	])

AC_CHECK_LIB(jsoncpp,main,[],[
	echo "Error! You need to have jsoncpp installed."
	exit -1
	])




AC_SUBST(INSTRUMENT_PATH,$PWD/instrument_modules/)



AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([molet_driver],[chmod +x molet_driver])

AC_OUTPUT
