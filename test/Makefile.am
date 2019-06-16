if BUNIT_TEST_PRESENT

check_PROGRAMS=libseq_unit_tests

TESTS=$(check_PROGRAMS)

AM_CXXFLAGS=-g
AM_LDFLAGS=-L../src/.libs -Wl,-rpath,../src/.libs
AM_LIBS=-lsequence

#if DEBUG
#AM_CXXFLAGS+=-g
#else
#AM_CXXFLAGS+=-DNDEBUG
#endif
#
#if PROFILING
#PROFILE= -pg
#else
#PROFILE=
#endif

LIBS+=$(AM_LIBS)

libseq_unit_tests_SOURCES=libseq_unit_tests.cc \
FastaConstructors.cc \
FastaIO.cc \
FastaOperations.cc \
AlignStreamTest.cc \
CountingOperators.cc \
PolyTableConversions.cc \
PolyTableTweaking.cc \
PolyTableBadBehavior.cc \
PolySitesIO.cc \
SimpleSNPIO.cc \
PolySIMtest.cc \
PolySNPtest.cc \
ComparisonsTest.cc \
AlignmentTest.cc \
fastqIO.cc \
fastqConstructors.cc \
SeqConversions.cc \
RedundancyCom95test.cc \
alphabets.cc \
polySiteVectorTest.cc \
PolyTableSliceTest.cc \
stateCounterTest.cc \
VariantMatrixTest.cc \
testAlleleCountMatrix.cc \
testClassicSummstats.cc \
testClassicSummstatsEmptyVariantMatrix.cc \
testLD.cc \
testGarudStatistics.cc \
msformatdata.cc

endif #if BUNIT_TEST_PRESENT
