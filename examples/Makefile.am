check_PROGRAMS= ms_to_VariantMatrix \
	nSL_from_ms \
	baseComp valid_dna translateTest \
	slidingWindow slidingWindow2 PolyTableIterators \
	ufs \
	msstats polySiteVector_test

ms_to_VariantMatrix_SOURCES=ms_to_VariantMatrix.cc
nSL_from_ms_SOURCES=nSL_from_ms.cc
baseComp_SOURCES=baseComp.cc
valid_dna_SOURCES=valid_dna.cc
translateTest_SOURCES=translateTest.cc
slidingWindow_SOURCES=slidingWindow.cc
slidingWindow2_SOURCES=slidingWindow2.cc
PolyTableIterators_SOURCES=PolyTableIterators.cc
ufs_SOURCES=ufs.cc
msstats_SOURCES=msstats.cc
polySiteVector_test_SOURCES=polySiteVector_test.cc

AM_CXXFLAGS=
AM_LIBS=-lsequence
AM_LDFLAGS=-L../src/.libs -Wl,-rpath,../src/.libs

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
