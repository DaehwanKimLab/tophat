define([gitrev], esyscmd([./git_revision.sh]))dnl
AC_INIT([tophat],[2.1.2],[tophat.cufflinks@gmail.com])
AC_DEFINE(SVN_REVISION, "gitrev", [SVN Revision])

AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
# AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE

AC_ARG_VAR(PYTHON, [python program])

# Make sure CXXFLAGS is defined so that AC_PROG_CXX doesn't set it.
CXXFLAGS="$CXXFLAGS"
CFLAGS="$CFLAGS"
AC_PROG_CXX
AC_LANG([C++])

# AC_LANG([C])

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL
AM_PATH_PYTHON([2.4])

m4_include([ax_boost_base.m4])
m4_include([ax_boost_thread.m4])
# CXXFLAGS="$CXXFLAGS $threadLib"
AX_BOOST_BASE([1.38.0])
AX_BOOST_THREAD
if test -z "$BOOST_THREAD_LIBS"; then
  AC_MSG_ERROR([boost.thread not found. Aborting.])
fi


# BAM related:
 ac_bam_path=samtools-0.1.18
 BAM_LIB="-lbam"
 BAM_LDFLAGS="-L./$ac_bam_path"
 BAM_CPPFLAGS="-I./$ac_bam_path"
 AC_SUBST(BAM_CPPFLAGS)
 AC_SUBST(BAM_LDFLAGS)
 AC_SUBST(BAM_LIB)


# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for header files.
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])
AC_C_CONST
AC_SYS_LARGEFILE

# Checks for libraries.
AC_CHECK_LIB([z], [gzread])

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset strdup strrchr strtol strsep])

# check the platform
AC_CANONICAL_HOST

# ---------------------------------------------------------------------
# Debug and profile
# ---------------------------------------------------------------------

# set CFLAGS and CXXFLAGS
#user_CFLAGS="${CXXFLAGS}"
user_CFLAGS=${CFLAGS}
generic_CFLAGS="-Wall -Wno-strict-aliasing -g -gdwarf-2 -Wuninitialized"
ext_CFLAGS=""
debug_CFLAGS=""
user_LDFLAGS="$LDFLAGS"

AC_ARG_ENABLE(intel64,      [  --enable-intel64        optimize for Intel64 CPU such as Xeon and Core2],
			  [ext_CFLAGS="${ext_CFLAGS} -mtune=nocona"], [])

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
	                      [enable debugging info (default is no)])],
              [], [enable_debug=no])
AC_ARG_ENABLE([optim],
              [AS_HELP_STRING([--enable-optim@<:@=0|1|2|3@:>@],
	                      [set optimization level (default is 3)])],
              [if test "x$enable_optim" = xyes; then enable_optim=3; fi],
              [enable_optim=3])

AS_IF([test "x$enable_optim" != xno], [ext_CFLAGS="$ext_CFLAGS -O$enable_optim"])	
AS_IF([test "x$enable_debug" = xyes],
      [debug_CFLAGS="-DDEBUG"],
      [debug_CFLAGS="-DNDEBUG"])
	  
CFLAGS="${generic_CFLAGS} ${ext_CFLAGS} ${user_CFLAGS} ${debug_CFLAGS}"
CXXFLAGS="$CFLAGS"
CXXFLAGS="$CXXFLAGS $BAM_CPPFLAGS $BOOST_CPPFLAGS -I./SeqAn-1.4.2"
LDFLAGS="$BAM_LDFLAGS $BOOST_LDFLAGS $user_LDFLAGS"

AM_INIT_AUTOMAKE([-Wall foreign tar-pax foreign])

# makefiles to configure
AC_CONFIG_FILES([Makefile src/Makefile])

# make it happen
AC_OUTPUT

# dump some configuration confirmations
echo \
"
-- ${PACKAGE_STRING} Configuration Results --
  C++ compiler:        ${CXX} ${CXXFLAGS}
  Linker flags:        ${LDFLAGS}
  BOOST libraries:     ${BOOST_THREAD_LIBS}"

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
  Email bug reports to <${PACKAGE_BUGREPORT}>.
"

if test x"${PYTHON}" = x":"  ||  ! test -x "${PYTHON}"; then
  echo "WARNING! python was not found and is required to run tophat"
  echo "  Please install python and point configure to the installed location"
fi
