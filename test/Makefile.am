if BUNIT_TEST_PRESENT

check_PROGRAMS=FastaConstructors FastaIO FastaExplicitIO FastaOperations \
	AlignStreamTest \
	testSimDataIO \
	CountingOperators \
	PolyTableConversions \
	PolyTableTweaking \
	PolyTableBadBehavior \
	PolySitesIO \
	SimpleSNPIO \
	PolySIMtest \
	PolySNPtest \
	ComparisonsTest \
	AlignmentTest \
	fastqIO \
	fastqConstructors \
	SeqConversions \
	RedundancyCom95test \
	alphabets \
	polySiteVectorTest \
	PolyTableSliceTest

TESTS=$(check_PROGRAMS)

#BOOST_LIBS=-lboost_unit_test_framework
#	SEQUENCE=-lsequence

AM_CXXFLAGS=
AM_LDFLAGS=-L../src/.libs -static
AM_LIBS=-lboost_unit_test_framework -lsequence

if HAVE_HTSLIB
AM_CXXFLAGS+=-DHAVE_HTSLIB
AM_LIBS+=-lhts
endif

if DEBUG
AM_CXXFLAGS+=-g
else
AM_CXXFLAGS+=-DNDEBUG
endif

if PROFILING
PROFILE= -pg
else
PROFILE=
endif

LIBS+=$(AM_LIBS)

FastaConstructors_SOURCES=FastaConstructors.cc
FastaIO_SOURCES=FastaIO.cc
FastaExplicitIO_SOURCES=FastaIO.cc
FastaOperations_SOURCES=FastaOperations.cc

AlignStreamTest_SOURCES=AlignStreamTest.cc

testSimDataIO_SOURCES=testSimDataIO.cc

CountingOperators_SOURCES=CountingOperators.cc

PolyTableConversions_SOURCES=PolyTableConversions.cc

PolyTableTweaking_SOURCES=PolyTableTweaking.cc

PolyTableBadBehavior_SOURCES=PolyTableBadBehavior.cc

PolySitesIO_SOURCES=PolySitesIO.cc

SimpleSNPIO_SOURCES=SimpleSNPIO.cc

PolySIMtest_SOURCES=PolySIMtest.cc

PolySNPtest_SOURCES=PolySNPtest.cc

ComparisonsTest_SOURCES=ComparisonsTest.cc

AlignmentTest_SOURCES=AlignmentTest.cc

fastqIO_SOURCES=fastqIO.cc

fastqConstructors_SOURCES=fastqConstructors.cc

SeqConversions_SOURCES=SeqConversions.cc

RedundancyCom95test_SOURCES=RedundancyCom95test.cc

alphabets_SOURCES=alphabets.cc

polySiteVectorTest_SOURCES=polySiteVectorTest.cc

PolyTableSliceTest_SOURCES=PolyTableSliceTest.cc

endif #if BUNIT_TEST_PRESENT
