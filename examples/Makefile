CC=cc
CXX=c++
CFLAGS = -O3 -Wall -pedantic -ansi -Wtraditional -I..
CXXFLAGS = -O3 -Wall -W -ansi -I..
LIBS = -L/usr/local/lib -lsequence -lz
REGEX= -lboost_regex
#all: codons.o gestimator.o msstats.o int_handler.o getopt.o getopt1.o baseComp.o valid_dna.o FastaExplicit.o translateTest.o slidingWindow.o slidingWindow2.o PolyTableIterators.o correlations.o critical_values.o ufs.o ms--.o msbeta.o freerec.o bottleneck.o fragments.o test_SimDataIO.o
all: codons.o msstats.o int_handler.o baseComp.o valid_dna.o FastaExplicit.o translateTest.o slidingWindow.o slidingWindow2.o PolyTableIterators.o correlations.o critical_values.o ufs.o ms--.o msbeta.o freerec.o bottleneck.o fragments.o test_SimDataIO.o
	$(CXX) $(CXXFLAGS)   ufs.o -o ufs $(LDFLAGS) $(LIBS)
	$(CXX) $(CXXFLAGS)   codons.o -o codons $(LDFLAGS) $(LIBS)
	$(CXX) $(CXXFLAGS)   msstats.o -o msstats $(LDFLAGS) $(LIBS)
#	$(CXX) $(CXXFLAGS)   gestimator.o int_handler.o getopt.o getopt1.o -o gestimator $(LDFLAGS) $(LIBS)
	$(CXX) $(CXXFLAGS)   baseComp.o -o baseComp $(LDFLAGS) $(LIBS)
	$(CXX) $(CXXFLAGS)   FastaExplicit.o -o FastaExplicit $(LDFLAGS) $(LIBS)
	$(CXX) $(CXXFLAGS)   translateTest.o -o translateTest $(LDFLAGS) $(LIBS)
	$(CXX) $(CXXFLAGS)   slidingWindow.o -o slidingWindow $(LDFLAGS) $(LIBS)
	$(CXX) $(CXXFLAGS)   slidingWindow2.o -o slidingWindow2 $(LDFLAGS) $(LIBS)
	$(CXX) $(CXXFLAGS)   PolyTableIterators.o -o PolyTableIterators $(LDFLAGS) $(LIBS)
	$(CXX) $(CXXFLAGS) correlations.o -o correlations -lgsl -lgslcblas $(LDFLAGS) $(LIBS)
	$(CXX) $(CXXFLAGS) critical_values.o -o critical_values $(LDFLAGS) $(LIBS)
	$(CXX) $(CXXFLAGS)   ms--.o -o ms-- $(LDFLAGS) $(LIBS) -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS)   msbeta.o -o msbeta $(LDFLAGS) $(LIBS) -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS)   freerec.o -o freerec $(LDFLAGS) $(LIBS) -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS)   bottleneck.o -o bottleneck $(LDFLAGS) $(LIBS) -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS)   fragments.o -o fragments $(LDFLAGS) $(LIBS) -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS)   valid_dna.o -o valid_dna $(LDFLAGS) $(LIBS) $(REGEX)
	$(CXX) $(CXXFLAGS)   test_SimDataIO.o -o test_SimDataIO $(LDFLAGS) $(LIBS) -lboost_iostreams

clean:
	rm -f *.o codons msstats baseComp valid_dna FastaExplicit translateTest \
	slidingWindow slidingWindow2 PolyTableIterators  correlations \
	critical_values ufs ms-- msbeta freerec bottleneck fragments

