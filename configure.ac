AC_PREREQ(2.59)

AC_INIT([libsequence], [1.8.4], [https://github.com/molpopgen/libsequence])
AC_CONFIG_SRCDIR([src/Seq/Seq.cc])
AC_CONFIG_SRCDIR([test/FastaIO.cc])
AC_CONFIG_SRCDIR([examples/msstats.cc])
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AM_MAINTAINER_MODE([disable])

LT_INIT
AC_PROG_LIBTOOL
AC_LANG(C++)
AC_CONFIG_FILES([Makefile src/Makefile Sequence/Makefile Sequence/bits/Makefile Sequence/Coalescent/Makefile Sequence/Coalescent/bits/Makefile Sequence/SummStats/Makefile test/Makefile examples/Makefile])

AC_ARG_ENABLE(debug,
[  --enable-debug    Turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=false])
AM_CONDITIONAL(DEBUG, test x$debug = xtrue)

AC_ARG_ENABLE(profiling,
[  --enable-profiling    Turn on profiling],
[case "${enableval}" in
  yes) profiling=true ;;
  no)  profiling=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-profiling) ;;
esac],[profiling=false])
AM_CONDITIONAL(PROFILING, test x$profiling = xtrue)


dnl zlib header
AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib headers missing - cannot continue])])

dnl boost unit test library
AC_CHECK_HEADER(boost/test/unit_test.hpp, BUNITTEST=1,[echo "boost/test/unit_test.hpp not found. Unit tests will not be compiled."])
AM_CONDITIONAL([BUNIT_TEST_PRESENT], test x$BUNITTEST = x1)

dnl htslib, http://www.htslib.org/
AC_CHECK_HEADER(htslib/bgzf.h,FOUND_HTSLIB_HEADERS=1,[AC_MSG_WARN([htslib/bgzf.h not found.  We will not be compiling features depending on this library])])
AM_CONDITIONAL([HAVE_HTSLIB],test x$FOUND_HTSLIB_HEADERS = x1 )
dnl zlib runtime
AC_CHECK_LIB([z],gzungetc,,[echo "zlib run time library not found";exit 1])

dnl htslib runtime
AC_CHECK_LIB([hts],bgzf_read,FOUND_HTSLIB_RUNTIME=1,[echo "htslib runtime not found.  Will not compile features requiring htslib"])
AM_CONDITIONAL([HAVE_HTSLIB],test x$FOUND_HTSLIB_RUNTIME = x1 )

AC_OUTPUT
