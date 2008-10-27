##-- TOPHAT CONFIGURE --##

# initialize
AC_PREREQ(2.61)
AC_INIT([TopHat],[0.7.0],[cole@cs.umd.edu],[tophat])
AC_PREFIX_DEFAULT([`pwd`])
AC_CONFIG_SRCDIR([src/reads.h])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall foreign])

# declare some precious variables
AC_ARG_VAR(PYTHON, [python program])

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_RANLIB
AC_PATH_PROG(PYTHON, [python], [:])

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
AC_CHECK_FUNCS([memset strdup strrchr strtol])

# check the platform
AC_CANONICAL_HOST

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
