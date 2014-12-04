check_PROGRAMS=FastaConstructors FastaIO FastaOperations \
	AlignStreamTest \
	testSimDataIO \
	CountingOperators \
	PolyTableConversions \
	PolyTableTweaking \
	PolySitesIO \
	SimpleSNPIO \
	PolySIMtest \
	PolySNPtest \
	ComparisonsTest \
	AlignmentTest

BOOST_LIBS=-lboost_unit_test_framework
SEQUENCE=-lsequence

AM_LDFLAGS=-L../src/.libs -static
LDADD=$(BOOST_LIBS) $(SEQUENCE)

FastaConstructors_SOURCES=FastaConstructors.cc
FastaIO_SOURCES=FastaIO.cc
FastaOperations_SOURCES=FastaOperations.cc

AlignStreamTest_SOURCES=AlignStreamTest.cc

testSimDataIO_SOURCES=testSimDataIO.cc

CountingOperators_SOURCES=CountingOperators.cc

PolyTableConversions_SOURCES=PolyTableConversions.cc

PolyTableTweaking_SOURCES=PolyTableTweaking.cc

PolySitesIO_SOURCES=PolySitesIO.cc

SimpleSNPIO_SOURCES=SimpleSNPIO.cc

PolySIMtest_SOURCES=PolySIMtest.cc

PolySNPtest_SOURCES=PolySNPtest.cc

ComparisonsTest_SOURCES=ComparisonsTest.cc

AlignmentTest_SOURCES=AlignmentTest.cc
