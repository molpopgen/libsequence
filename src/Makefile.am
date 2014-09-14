lib_LTLIBRARIES=libsequence.la

libsequence_la_SOURCES=  Grantham.cc\
	PathwayHelper.cc\
	SummStats/Recombination.cc\
	SimParams.cc\
	Translate.cc\
	Comeron95.cc\
	GranthamWeights.cc\
	SummStats/PolySNP.cc\
	SummStats/PolySIM.cc\
	RedundancyCom95.cc\
	SingleSub.cc\
	TwoSubs.cc\
	stateCounter.cc\
	SummStats/FST.cc\
	Comparisons.cc\
	SimpleSNP.cc\
	PolyTable.cc\
	PolyTableFunctions.cc\
	PolyTableManip.cc\
	Seq/Seq.cc\
	ComplementBase.cc\
	Sites.cc\
	Unweighted.cc\
	Seq/Fasta.cc\
	Kimura80.cc\
	PolySites.cc\
	SimData.cc\
	SimDataIO.cc\
	ThreeSubs.cc\
	CodonTable.cc\
	Specializations.cc\
	Seq/FastaExplicit.cc\
	SeqExceptions.cc\
	SeqConstants.cc\
	shortestPath.cc\
	SummStats/HKA.cc\
	SummStats/Snn.cc\
	Coalescent/CoalescentMutation.cc\
	Coalescent/CoalescentRecombination.cc\
	Coalescent/CoalescentSimTypes.cc\
	Coalescent/CoalescentTreeOperations.cc\
	Coalescent/CoalescentInitialize.cc\
	Coalescent/CoalescentCoalesce.cc\
	Coalescent/CoalescentFragmentsRescaling.cc\
	sam/samrecord.cc\
	sam/samflag.cc\
	sam/samfunctions.cc\
	Ptable.cc \
	SummStats/SummStats.cc \
	SummStats/nSL.cc \
	PolyData.cc

AM_LDFLAGS = -version-info 19:0:0

if DEBUG
CXXOTHERFLAGS = -g
COTHERFLAGS = -g
else
CXXOTHERFLAGS = -DNDEBUG -O2
COTHERFLAGS = -DNDEBUG -O2
endif

if PROFILING
PROFILE= -pg
else
PROFILE=
endif

#if GFIVEOPT
#GFIVEFLAGS= -mcpu=G5 -mtune=G5 -mpowerpc64 -mpowerpc-gpopt
#else
#GFIVEFLAGS=
#endif

#if GFOUROPT
#GFOURFLAGS= -mcpu=G4 -mtune=G4 -mpowerpc -mpowerpc-gpopt
#else
#GFOURFLAGS=
#endif

#if HAVE_STRSTREAM
#STRSTREAM = -DHAVE_STRSTREAM
#else
#STRSTREAM =
#endif

#if HAVE_SSTREAM
#SSTREAM = -DHAVE_SSTREAM
#else
#SSTREAM =
#endif

#if HAVE_LIMITS
#LIMITS= -DHAVE_LIMITS
#else
#LIMITS =
#endif

#if HAVE_CLIMITS
#CLIMITS= -DHAVE_CLIMITS
#else
#CLIMITS =
#endif

#if HAVE_IEEEFP
#IEEEFPH= -DHAVE_IEEEFP_H
#else
#IEEEFPH=
#endif

#some bits for the compile
#CXXLINK = 
AM_CXXFLAGS = $(CXXOTHERFLAGS) $(GFIVEFLAGS) $(GFOURFLAGS) $(PROFILE) $(STRSTREAM) $(SSTREAM) $(LIMITS) $(CLIMITS) $(IEEEFPH)  -Wall -W -Woverloaded-virtual  -Wnon-virtual-dtor -Wcast-qual -Wconversion -Wsign-conversion -Wsign-promo -Wsynth -ffor-scope -pedantic -funroll-loops -ansi
CFLAGS =  $(COTHERFLAGS) $(PROFILE) -Wall $(INCLUDES) -pedantic -ansi
