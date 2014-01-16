CC=cc
CXX=c++
CFLAGS = -O3 -Wall -pedantic -ansi -Wtraditional -I..
CXXFLAGS = -O3 -Wall -W -ansi -I..
LIBS = -L/usr/local/lib -lsequence -lz
REGEX= -lboost_regex
all: codons.o gestimator.o msstats.o int_handler.o getopt.o getopt1.o baseComp.o valid_dna.o FastaExplicit.o translateTest.o slidingWindow.o slidingWindow2.o PolyTableIterators.o correlations.o critical_values.o ufs.o ms--.o msbeta.o freerec.o bottleneck.o fragments.o test_SimDataIO.o
	$(CXX) $(CXXFLAGS)   ufs.o -o ufs $(LIBS)
	$(CXX) $(CXXFLAGS)   codons.o -o codons $(LIBS)
	$(CXX) $(CXXFLAGS)   msstats.o -o msstats $(LIBS)
	$(CXX) $(CXXFLAGS)   gestimator.o int_handler.o getopt.o getopt1.o -o gestimator $(LIBS)
	$(CXX) $(CXXFLAGS)   baseComp.o -o baseComp $(LIBS)
	$(CXX) $(CXXFLAGS)   FastaExplicit.o -o FastaExplicit $(LIBS)
	$(CXX) $(CXXFLAGS)   translateTest.o -o translateTest $(LIBS)
	$(CXX) $(CXXFLAGS)   slidingWindow.o -o slidingWindow $(LIBS)
	$(CXX) $(CXXFLAGS)   slidingWindow2.o -o slidingWindow2 $(LIBS)
	$(CXX) $(CXXFLAGS)   PolyTableIterators.o -o PolyTableIterators $(LIBS)
	$(CXX) $(CXXFLAGS) correlations.o -o correlations -lgsl -lgslcblas $(LIBS)
	$(CXX) $(CXXFLAGS) critical_values.o -o critical_values $(LIBS)
	$(CXX) $(CXXFLAGS)   ms--.o -o ms-- $(LIBS) -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS)   msbeta.o -o msbeta $(LIBS) -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS)   freerec.o -o freerec $(LIBS) -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS)   bottleneck.o -o bottleneck $(LIBS) -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS)   fragments.o -o fragments $(LIBS) -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS)   valid_dna.o -o valid_dna $(LIBS) $(REGEX)
	$(CXX) $(CXXFLAGS)   test_SimDataIO.o -o test_SimDataIO $(LIBS) -lboost_iostreams

clean:
	rm -f *.o codons msstats gestimator baseComp valid_dna FastaExplicit translateTest \
	slidingWindow slidingWindow2 PolyTableIterators  correlations \
	critical_values ufs ms-- msbeta freerec bottleneck fragments

