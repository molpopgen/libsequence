
#include <iostream>
#ifndef _UNISTD_H_
#include "getopt.h"
#endif
#include <vector>
#include <glob.h>
#include <time.h>
#include <cstdio>
#include <memory>
#include <Sequence/PolySites.hpp>

//if using GNU C++, use explicit instantiations
//of routines for FASTA I/O.  
#ifdef __GNUG__
#include <Sequence/FastaExplicit.hpp>
#else
//On other systems, rely on compile-time (implicit)
//template instantiation
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#endif

#include <Sequence/Comeron95.hpp>
#include "int_handler.hpp"

/*!   \include gestimator.cc
 */

/*
  Calculate Ka and Ks between aligned sequences using Comeron's (1995)
  method.  Basically makes heavy use of Sequence::Comeron95
*/
using namespace std;
using namespace Sequence;
using namespace Alignment;

//struct params keeps all the command-line arguments
//together
struct params {
  char *infileglob;
  char *outfile;
  int maxhits;
  bool verbose;
  bool remove_all_gaps;
};

//prototypes
void parseargs(int argc, char *argv[],params *args);
void process(glob_t *files, params *args, ostream &ofstr);
void usage(void);


int main(int argc, char *argv[]) {
  glob_t files;
  ofstream ofstr;
  params args;

  parseargs(argc,argv,&args);

  if(args.infileglob == NULL) {
    cout << "fatal error: no infile(s) specified\n";
    usage();
    exit(10);
  }

  //filenames are obtained via a pattern passed from
  //the shell, i.e. compute -i '*.fasta', where the 
  //pattern in single quotes is what gets processed
  //see man glob(3)
#if defined (__SVR4) && defined (__sun)
  glob(args.infileglob,GLOB_ERR,NULL,&files);
#else
  glob(args.infileglob,GLOB_TILDE|GLOB_ERR,NULL,&files);
#endif
  if (args.outfile != NULL)
    ofstr.open(args.outfile);
  if (args.verbose == 1) {
    if (args.outfile != NULL){
      ofstr <<"seq1\tseq2\tKa\tKs\tomega\tp0\tp2S\tp2V\tp4\tq0\tq2S\tq2V\tq4\t";
      ofstr <<"Aa\tAs\tBa\tBs\tL0\tL2S\tL2V\tL4"<<endl;
    } else {
      cout <<"seq1\tseq2\tKa\tKs\tomega\tp0\tp2S\tp2V\tp4\tq0\tq2S\tq2V\tq4\t";
      cout <<"Aa\tAs\tBa\tBs\tL0\tL2S\tL2V\tL4"<<endl;
    }
  }

  signal(SIGINT,cntrl_c_handler);
  process(&files,&args,ofstr);

  globfree(&files);
  exit(0);
}


void parseargs(int argc, char *argv[],params *args)
{
  if(argc==1){
    usage();
    exit(1);
  }

  //assign some defaults
  args->infileglob = NULL;
  args->outfile = NULL;
  args->maxhits = 3;
  args->verbose = 0;
  args->remove_all_gaps =0;

  int c;

  while ((c = getopt (argc, argv, "i:o:m:vg")) != -1)
    {
      switch (c)
	{
	case 'i':
	  args->infileglob = optarg;
	  break;
	case 'o':
	  args->outfile = optarg;
	  break;
	case 'm':
	  args->maxhits = atoi(optarg);
	  break;
	case 'v':
	  args->verbose = 1;
	  break;
	case 'g':
	  args->remove_all_gaps = 1;
	  break;
	default:
	  usage();
	  exit(1);
	  break;
	}
    }
}

void process(glob_t *files, params *args, ostream &ofstr)
{
  for(unsigned int i=0;i<files->gl_pathc;++i)  {
    char *infile = files->gl_pathv[i];
    vector<Fasta > data;//store the sequence objects
    bool fileValid = 1;
    try {
      GetData(data,infile);
      //GetData(...) is from namespace Alignment
      //there are 2 forms of the function:
      //GetData(vector<T *> &data, const char *infile), which reads data
      //from a file identified by a C-style string
      //GetData(vector<T *> &data, istream &in), which can read
      //from any open istream
    } catch (badFormat &b)
      {
	//badFormat is thrown if
	//the input is not in the expected format,
	//such as a FASTA file that does not begin
	//with '>'
	cerr <<"error: processing of file " << infile << " threw the following\n";
	cerr<<"exception: ";
	b.print(cerr);
	cerr << endl;
	fileValid = 0;
      }
    if (fileValid) {
      if (args->remove_all_gaps == 1) {
	if (Gapped(data))
	  RemoveGaps(data);
      }
      for (unsigned int i = 0 ; i < data.size()-1 ; ++i) {
	for (unsigned int j = i+1 ; j < data.size() ; ++j) {
	  try {
	    auto_ptr<Comeron95> C(new Comeron95(&data[i],&data[j],args->maxhits));
	    if (args->outfile == NULL) {
	      cout.precision(4);
	      cout << data[i].GetName()<< '\t';
	      cout << data[j].GetName() << '\t';
	      cout << C->ka() << '\t';
	      cout << C->ks() << '\t';
	      cout << C->ratio();
	      if (args->verbose) {
		cout <<'\t';
		cout <<C->P0() << '\t' <<C->P2S()<<'\t';
		cout <<C->P2V() <<'\t' <<C->P4()<<'\t';
		cout <<C->Q0() << '\t' <<C->Q2S()<<'\t';
		cout <<C->Q2V() << '\t' <<C->Q4()<<'\t';
		cout <<C->aa()<<'\t'<<C->as()<<'\t'<<C->ba()<<'\t'<<C->bs()<<'\t';
		cout.precision(6);
		cout <<C->L0() << '\t' <<C->L2S()<<'\t';
		cout <<C->L2V() << '\t' <<C->L4()<<endl;
	      } else
		cout << endl;
	    }
	    else {
	      ofstr.precision(4);
	      ofstr << data[i].GetName()<< '\t';
	      ofstr << data[j].GetName() << '\t';
	      ofstr << C->ka() << '\t';
	      ofstr << C->ks() << '\t';
	      ofstr << C->ratio();
	      if (args->verbose) {
		ofstr<<'\t';
		ofstr <<C->P0() << '\t' <<C->P2S()<<'\t';
		ofstr <<C->P2V() <<'\t' <<C->P4()<<'\t';
		ofstr <<C->Q0() << '\t' <<C->Q2S()<<'\t';
		ofstr <<C->Q2V() << '\t' <<C->Q4()<<'\t';
		ofstr <<C->aa()<<'\t'<<C->as()<<'\t'<<C->ba()<<'\t'<<C->bs()<<'\t';
		ofstr.precision(6);
		ofstr <<C->L0() << '\t' <<C->L2S()<<'\t';
		ofstr <<C->L2V() << '\t' <<C->L4()<<endl;
	      } else
		ofstr<<endl;

	    }
	  } catch(Sequence::SeqException &e) {
	    cout << infile <<": ";
	    e.print(cout);
	    cout << endl;
	  }
	}
      }
    }
    data.clear();
  }

}


void usage(void){
  cerr<<"usage: ";
  cerr<<"gestimator -i infile"<<endl;
  cerr<<"\tother options\n";
  cerr<<"\t\t-o outfile : write results to outfile\n";
  cerr<<"\t\t-v : get verbose output\n";
  cerr<<"\t\t-m : max # of hits allowed per codon (default = 3)\n";
  cerr <<"\t\t-g : remove gaps from the whole aligment before calculating (default = FALSE)\n";
  cerr << endl;
}
