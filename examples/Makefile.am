check_PROGRAMS= baseComp valid_dna translateTest \
	slidingWindow slidingWindow2 PolyTableIterators  correlations \
	critical_values ufs msmm freerec bottleneck fragments \
	msstats Ptable_test

baseComp_SOURCES=baseComp.cc
valid_dna_SOURCES=valid_dna.cc
translateTest_SOURCES=translateTest.cc
slidingWindow_SOURCES=slidingWindow.cc
slidingWindow2_SOURCES=slidingWindow2.cc
PolyTableIterators_SOURCES=PolyTableIterators.cc
correlations_SOURCES=correlations.cc
critical_values_SOURCES=critical_values.cc
ufs_SOURCES=ufs.cc
msmm_SOURCES=msmm.cc
freerec_SOURCES=freerec.cc
bottleneck_SOURCES=bottleneck.cc
fragments_SOURCES=fragments.cc
msstats_SOURCES=msstats.cc
Ptable_test_SOURCES=Ptable_test.cc

AM_CXXFLAGS=
AM_LDFLAGS=-L../src/.libs -static

AM_LIBS=-lsequence
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
