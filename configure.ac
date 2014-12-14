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

dnl check for things that we need

dnl boost (only needed prior to 1.8.3)
dnl AC_CHECK_HEADER(boost/type_traits.hpp, LIBSEQUENCE_USE_BOOST=1 , [AC_MSG_ERROR([boost::type_traits required - either install boost (see http://www.boost.org) or write a replacement])])
dnl AC_CHECK_HEADER(boost/utility.hpp, LIBSEQUENCE_USE_BOOST=1 , [AC_MSG_ERROR([boost::noncopyable required (boost/utility.hpp) - either install boost (see http://www.boost.org) or write a replacement])])
dnl AC_CHECK_HEADER(boost/static_assert.hpp, LIBSEQUENCE_USE_BOOST=1 , [AC_MSG_ERROR([BOOST_STATIC_ASSERT required - either install boost (see http://www.boost.org) or write a replacement])])
dnl AC_CHECK_HEADER(boost/tuple/tuple.hpp,,[AC_MSG_ERROR([boot::tuple required - either install boost (see http://www.boost.org) or write a replacement])])

dnl gsl.  dependency removed in 1.8.3
dnl AC_CHECK_HEADER(gsl/gsl_rng.h, LIBSEQUENCE_USE_GSL=1 , [AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install GSL from http://www.gnu.org/software/gsl])])

dnl zlib header
AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib headers missing - cannot continue])])

dnl htslib, http://www.htslib.org/
AC_CHECK_HEADER(htslib/bgzf.h,FOUND_HTSLIB_HEADERS=1,[AC_MSG_WARN([htslib/bgzf.h not found.  We will not be compiling features depending on this library])])
AM_CONDITIONAL([HAVE_HTSLIB],test x$FOUND_HTSLIB_HEADERS = x1 )
dnl zlib runtime
AC_CHECK_LIB([z],gzungetc,,[echo "zlib run time library not found";exit 1])

dnl htslib runtime
AC_CHECK_LIB([hts],bgzf_read,FOUND_HTSLIB_RUNTIME=1,[echo "htslib runtime not found.  Will not compile features requiring htslib"])
AM_CONDITIONAL([HAVE_HTSLIB],test x$FOUND_HTSLIB_RUNTIME = x1 )
dnl AM_CONFIG_HEADER(config.h)

AC_OUTPUT
