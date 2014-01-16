#include <Sequence/Coalescent/Coalescent.hpp>
#include <Sequence/RNG/gsl_rng_wrappers.hpp>
#include <boost/bind.hpp>
#include <boost/static_assert.hpp>
#include <algorithm>
#include <numeric>
#include <functional>

int main(int argc, char **argv)
{
  std::vector< std::pair<int,int> > fragments;
  //a list of fragments is a vector of pairs
  //each pair contains the distance to the next fragment (in bp),
  //and then the length of the fragment (in bp)
  fragments.push_back(std::make_pair(0,500));
  fragments.push_back(std::make_pair(1000,500));
  fragments.push_back(std::make_pair(1000,500));
  fragments.push_back(std::make_pair(1000,500));
  std::vector< std::pair<double,double> > sample_scale,mutation_scale;

  //convert bp to units on the interval 0,1
  Sequence::calculate_scales(fragments,&sample_scale,&mutation_scale);


  const int slength = Sequence::sample_length(fragments);
  const int tlength = Sequence::total_length(fragments);

  std::vector<Sequence::chromosome> isample = 
    Sequence::init_sample( std::vector<int>(1,10),slength );
  Sequence::marginal imarg = Sequence::init_marginal(10);

  //make a genetic map for a region
  std::vector<double> genetic_map;
  double rho_bw_loci = 100.;
  double rho_locus = 50.;
  for(unsigned locus = 0 ; locus < fragments.size() ; ++locus)
    {
      int length = fragments[locus].second;
      for(int site = 0 ; site < length ; ++site)
	{
	  double trho = (site < length-1 ) ? rho_locus/double(length-1) : 
	    ( (locus < fragments.size()-1) ? rho_bw_loci : 0. ) ;
	  genetic_map.push_back( trho );
	}
    }

  const double total_rho = std::accumulate(genetic_map.begin(),genetic_map.end(),0.);
  std::vector<double> genetic_map_pdf(genetic_map);
  std::transform(genetic_map_pdf.begin(),genetic_map_pdf.end(),genetic_map_pdf.begin(),
		 std::bind2nd(std::divides<double>(),total_rho));
  const int total_size = genetic_map.size()-1;

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,time(0));

  Sequence::gsl_uniform01 uni01(r); 
  Sequence::gsl_uniform uni(r);     
  Sequence::gsl_exponential expo(r);
  Sequence::gsl_poisson poiss(r);   

  const int nsam = 10;
  const int NRUNS = 1000;
  int run = 0;
  std::vector<double> reclens;
  while(run++ < NRUNS)
    {
      std::vector<Sequence::chromosome> sample(isample);
      Sequence::arg sample_history(1,imarg);
      int NSAM = nsam;
      int nlinks = NSAM*total_size;
      double t = 0.;
      while( NSAM > 1 )
	{
	  double rc = double(NSAM*(NSAM-1));
	  double rrec = Sequence::integrate_genetic_map(sample,NSAM,genetic_map,&reclens);
	  double tc = expo(1./rc),tr = expo(1./rrec);
	 
	  t += std::min(tc,tr);
	  if( tc < tr ) //CA
	    {
	      //cerr << "CA" << '\n';
	      std::pair<int,int> two = Sequence::pick2(uni,NSAM);
	      NSAM -= Sequence::coalesce(t,nsam,NSAM,two.first,two.second,slength,
			       &nlinks,&sample,&sample_history);
	    }
	  else //RE
	    {
	      //cerr << "REC" << '\n';
	      std::pair<int,int> pos_rec = Sequence::pick_spot(uni01,
						rrec,
						reclens,
						sample.begin(),
						NSAM,&genetic_map_pdf[0]);
	      nlinks -= Sequence::crossover(NSAM,pos_rec.first,pos_rec.second,
					    &sample,&sample_history);
	      NSAM++;
	    }
	  if( unsigned(NSAM) < sample.size()/5 )
	    sample.erase(sample.begin()+NSAM+1,sample.end());
	}
      Sequence::minimize_arg(&sample_history);      
      Sequence::SimData d = infinite_sites_sim_data(poiss,uni,
					  slength,sample_history,10.);
      Sequence::rescale_mutation_positions(&d,sample_scale,mutation_scale);
      std::cout  << d << std::endl;
    }
}
