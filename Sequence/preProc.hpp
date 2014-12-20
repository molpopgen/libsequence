#ifndef __SEQUENCE_PREPROC_HPP__
#define __SEQUENCE_PREPROC_HPP__

#include <memory>
#include <Sequence/Ptable8.hpp>
#include <Sequence/PolyTable.hpp>
#include <Sequence/SummStats/Uhaps.hpp>
#include <Sequence/stateCounter.hpp>

namespace Sequence
{
  //! fwd declaration
  struct preProcImpl;
  class Ptable;
  class preProc
  {
  private:
    std::unique_ptr<preProcImpl> __impl;
  public:
    mutable Ptable8 ptable;
    mutable Uhaps uhaps;
    using state_t = std::vector<stateCounter>;
    using dstate_t = std::vector<std::pair<bool,stateCounter> >;
    //! Unpolarized states per site.
    mutable state_t states;
    /*! 
      Polarized states per site.		       
      \note Will be empty if no info ancestral states are 
      provided
    */
    mutable dstate_t dstates;
  
    preProc() = default;
    preProc( const Ptable &, const bool & haveAncStates = 0, const std::string::size_type & anc = 0 );
    preProc( Ptable &, const bool & haveAncStates = 0, const std::string::size_type & anc = 0 );
    preProc( const PolyTable &, const bool & haveAncStates = 0, const PolyTable::size_type & anc = 0 );
    preProc( PolyTable &&, const bool & haveAncStates = 0, const PolyTable::size_type & anc = 0 );
    ~preProc(void);// = default;
  };
}


#endif
