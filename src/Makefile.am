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
	SeqConstants.cc\
	shortestPath.cc\
	SummStats/HKA.cc\
	SummStats/Snn.cc\
	polySiteVector.cc \
	SummStats/SummStats.cc \
	SummStats/nSL.cc \
	SummStats/Garud.cc\
	IOhelp.cc \
	SeqAlphabets.cc \
	SummStats/lHaf.cc \
	threading.cc \
	variant_matrix/VariantMatrix.cc \
	variant_matrix/VariantMatrixViews.cc \
	variant_matrix/StateCounts.cc \
	variant_matrix/filtering.cc \
	summstats/thetapi.cc \
	summstats/thetaw.cc \
	summstats/tajd.cc \
	summstats/thetah_thetal.cc \
	summstats/faywuh.cc \
	summstats/hprime.cc \
	summstats/nvariablesites.cc \
	summstats/allele_counts.cc \
	summstats/haplotype_statistics.cc \
	summstats/ld.cc \
	summstats/rmin.cc \
	summstats/nsl.cc \
	summstats/garud.cc \
	summstats/generic.cc \
	summstats/auxillary.cc


AM_LDFLAGS=-version-info 20:0:0

AM_CXXFLAGS= -Wall -W -Woverloaded-virtual  -Wnon-virtual-dtor -Wcast-qual -Wconversion -Wsign-conversion -Wsign-promo -Wsynth

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

