##-- TOPHAT CONFIGURE --##

# initialize
AC_PREREQ(2.61)
AC_INIT([TopHat],[0.8.2],[cole@cs.umd.edu],[tophat])
AC_PREFIX_DEFAULT([`pwd`])
AC_CONFIG_SRCDIR([src/reads.h])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall gnits])

AC_ARG_ENABLE(debug,
	AC_HELP_STRING([--enable-debug],[build debug library (default=no)]),,enable_debug=no)

# declare some precious variables
AC_ARG_VAR(PYTHON, [python program])
AC_ARG_VAR(PERL, [perl program])
AC_ARG_VAR(TOPHAT_SHELL, [bash program])

CFLAGS=" $CFLAGS"
CXXFLAGS=" $CXXFLAGS"
# Checks for programs.
AC_LANG(C++)
AC_PROG_CC
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_RANLIB
AC_CHECK_TOOL(AR,ar,false)
AC_PATH_PROG(PERL, [perl], [:])
AC_PATH_PROG(TOPHAT_SHELL, bash, no)
AM_PATH_PYTHON([2.4])

# ---------------------------------------------------------------------
# Debug and profile
# ---------------------------------------------------------------------

# these are first because the SET AM_C.*FLAGS, not just append
if test x$enable_debug = xyes; then
  GCJFLAGS="-g -O0"
  if test x$ac_cv_prog_cc_g = xyes; then
    AM_CFLAGS="-g3 -O0 -Wall"
  else
    AM_CFLAGS="-O0 -Wall"
  fi
  if test x$ac_cv_prog_cxx_g = xyes; then
    AM_CXXFLAGS="-g3 -O0 -Wall"
  else
    AM_CXXFLAGS="-O0 -Wall"
  fi
else
  AM_CFLAGS="-O3"
  AM_CXXFLAGS="-O3"
fi


AM_CXXFLAGS="$AM_CXXFLAGS -I./SeqAn-1.1"

# Checks for libraries.
AC_CHECK_LIB([z], [gzread])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_INT32_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_SYS_LARGEFILE

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset strdup strrchr strtol strsep])

# check the platform
AC_CANONICAL_HOST

AC_SUBST(CC)
AC_SUBST(AM_CXXFLAGS)

# makefiles to configure
AC_CONFIG_FILES([Makefile src/Makefile])

# make it happen
AC_OUTPUT

# dump some configuration confirmations
echo \
"
-- ${PACKAGE_STRING} Configuration Results --
  C compiler:          ${CC} ${CFLAGS}
  C++ compiler:        ${CXX} ${CXXFLAGS}"

if test x"${GCC}" = x"yes" ; then
   gcc_version=`${CC} --version | head -n 1`
   echo "  GCC version:         ${gcc_version}"
else
   gcc_version=''
fi

echo \
"  Host System type:    ${host}
  Install prefix:      ${prefix}
  Install eprefix:     ${exec_prefix}

  See config.h for further configuration information.
  Email <${PACKAGE_BUGREPORT}> with questions and bug reports.
"

if test x"${PYTHON}" = x":"  ||  ! test -x "${PYTHON}"; then
  echo "WARNING! python was not found and is required to run tophat"
  echo "  Please install python and point configure to the installed location"
fi
