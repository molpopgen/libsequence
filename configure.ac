AC_PREREQ(2.59)

AC_INIT([libsequence], [1.9.8], [https://github.com/molpopgen/libsequence])
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
AC_CONFIG_FILES([Makefile src/Makefile Sequence/Makefile Sequence/bits/Makefile Sequence/SummStatsDeprecated/Makefile
				 Sequence/variant_matrix/Makefile Sequence/summstats/Makefile test/Makefile examples/Makefile doc/libsequence.doxygen])

dnl AC_ARG_ENABLE(debug,
dnl [  --enable-debug    Turn on debugging],
dnl [case "${enableval}" in
dnl   yes) debug=true ;;
dnl   no)  debug=false ;;
dnl   *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
dnl esac],[debug=false])
dnl AM_CONDITIONAL(DEBUG, test x$debug = xtrue)
dnl 
dnl AC_ARG_ENABLE(profiling,
dnl [  --enable-profiling    Turn on profiling],
dnl [case "${enableval}" in
dnl   yes) profiling=true ;;
dnl   no)  profiling=false ;;
dnl   *) AC_MSG_ERROR(bad value ${enableval} for --enable-profiling) ;;
dnl esac],[profiling=false])
dnl AM_CONDITIONAL(PROFILING, test x$profiling = xtrue)


dnl zlib header
dnl AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib headers missing - cannot continue])])

dnl zlib runtime
dnl AC_CHECK_LIB([z],gzungetc,,[echo "zlib run time library not found";exit 1])

dnl boost unit test library
AC_CHECK_HEADER(boost/test/unit_test.hpp, BUNITTEST=1,[echo "boost/test/unit_test.hpp not found. Unit tests will not be compiled."])
AM_CONDITIONAL([BUNIT_TEST_PRESENT], test x$BUNITTEST = x1)

dnl check for Intel TBB headers that we need
dnl AC_CHECK_HEADER(tbb/parallel_for.h,,[AC_MSG_ERROR([tbb/parallel_for.h not found.  Please either install Intel's TBB library or make sure your CXXFLAGS are set correctly])])
dnl AC_CHECK_HEADER(tbb/task_scheduler_init.h,,[AC_MSG_ERROR([tbb/task_scheduler_init.h not found.  Please either install Intel's TBB library or make sure your CXXFLAGS are set correctly])])

dnl check for Intel TBB library
dnl AC_CHECK_LIB([tbb],main,FOUND_TBB_RUNTIME=1,[AC_MSG_ERROR([tbb runtime library not found.  Please install Intel's TBB library.])])
AC_OUTPUT
