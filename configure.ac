
AC_INIT([miraculix], 0.9)
CXX=`"${R_HOME}/bin/R" CMD config CXX`
AC_PROG_CXX
AC_LANG(C++)

AC_MSG_CHECKING([which SIMD flags are available])

##
CPPFLAGS="" ## to be deleted if the above does not cause any problems
if test "$CPPFLAGS" != "" ; then
  AC_MSG_RESULT([$CPPFLAGS])
else
  AC_MSG_RESULT([none])
fi

MIRACULIX_CXXFLAGS="$CXX_FLAGS $CPPFLAGS"
AC_SUBST(MIRACULIX_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
