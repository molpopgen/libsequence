lib_LTLIBRARIES=libsequence.la
bin_PROGRAMS=libsequenceConfig

libsequenceConfig_SOURCES=libsequenceConfig.cc

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
	Seq/Seq.cc\
	ComplementBase.cc\
	Sites.cc\
	Unweighted.cc\
	Seq/Fasta.cc\
	Seq/fastq.cc\
	Kimura80.cc\
	PolySites.cc\
	SimData.cc\
	SimDataIO.cc\
	ThreeSubs.cc\
	CodonTable.cc\
	Specializations.cc\
	Seq/FastaExplicit.cc\
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
	hts/samrecord.cc\
	hts/samflag.cc\
	hts/samfunctions.cc\
	polySiteVector.cc \
	SummStats/SummStats.cc \
	SummStats/nSL.cc \
	SummStats/Garud.cc\
	IOhelp.cc \
	hts/bamrecord.cc \
	hts/bamreader.cc \
	SeqAlphabets.cc \
	SummStats/lHaf.cc \
	SummStats/variantCounts.cc \
	SummStats/classic.cc \
	threading.cc \
	variant_matrix/VariantMatrix.cc \
	variant_matrix/VariantMatrixViews.cc \
	variant_matrix/StateCounter.cc


AM_LDFLAGS=-version-info 20:0:0

AM_CXXFLAGS= -Wall -W -Woverloaded-virtual  -Wnon-virtual-dtor -Wcast-qual -Wconversion -Wsign-conversion -Wsign-promo -Wsynth

if HAVE_HTSLIB
AM_CXXFLAGS+=-DHAVE_HTSLIB
AM_LDFLAGS+=-lhts
endif

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

#Link to Intel's TBB runtime library
AM_LDFLAGS+=-ltbb
