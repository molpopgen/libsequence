#ifdef HAVE_HTSLIB //Will only compile if ./configure detects htslib

#include <Sequence/bamrecord.hpp>
#include <algorithm>
#include <cstring>
#include <cassert>

using namespace std;

namespace {

  //The valid tag types and the number
  static size_t nValTypes = 12; //C++ 1 extra tradition so we have a test-for-end with std::find
  static const char ValTypes[12] = {'A','c','C','s','S','i','I','f','Z','H','B',char(255)}; //char(255) is not an allowed value, and terminates the array

  size_t auxTagSize(const char & c)
  {
    switch(c)
      {
      case 'c':
	return sizeof(int8_t);
	break;
      case 'C':
	return sizeof(uint8_t);
	break;
      case 's':
	return sizeof(int16_t);
	break;
      case 'S':
	return sizeof(uint16_t);
	break;
      case 'i':
	return sizeof(int32_t);
	break;
      case 'I':
	return sizeof(uint32_t);
	break;
      default:
	return 0;
	break;
      }
    return 0;
  }

  size_t sizeofTag( const char & valtype,
		    const char * auxpos,
		    const char * auxposend)
  {
    if( std::find(&ValTypes[0],&ValTypes[0]+nValTypes,valtype) ==
	&ValTypes[0]+nValTypes ) return 0; //error

    if( valtype == 'A' ) return sizeof(char);
    else if (valtype == 'Z')
      {
	return size_t(std::find(auxpos,auxposend,'\0') - auxpos)+1;
      }
    else if (valtype == 'B' ) //Experimental/untested
      {
	char Btype = *auxpos;
	std::int32_t Bsize = *(const std::int32_t*)(auxpos+1);
	return size_t(sizeof(char) + sizeof(std::int32_t) + size_t(Bsize)*sizeof(auxTagSize(Btype)));
      }
    else if (valtype == 'H')
      {
	return 0;
      }
    else if (valtype == 'f') return sizeof(float);

    return auxTagSize(valtype);
  }
}

using U32 = std::uint32_t;
using I32 = std::int32_t;
using U16 = std::uint16_t;
using I16 = std::int16_t;
using U8 = std::uint8_t;
using I8 = std::int8_t;

namespace Sequence 
{
  namespace bamutil {
    const char int2seq[16] = {'=','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};
    const char bamCig[9] = {'M','I','D','N','S','H','P','=','X'};
  }
  using bamutil::int2seq;
  using bamutil::bamCig;
  bamaux::bamaux( ) : size(0),
		      value_type(char()),
		      //tag(),
		      value(nullptr)
  {
  }

  bamaux::bamaux( size_t __size,
		  //std::unique_ptr<char[]> & __tag,
		  char __tag[3],
		  char __value_type,
		  std::unique_ptr<char[]> & __value) : size(std::move(__size)),
						       value_type(std::move(__value_type)),
						       value(std::move(__value))
  {
    tag[0]=__tag[0];
    tag[2]=__tag[1];
    tag[1]=__tag[2];
  }

  bamaux::bamaux( bamaux && ba ) : size(std::move(ba.size)),
				   value_type(std::move(ba.value_type)),
				   value(std::move(ba.value))
  {
    tag[0]=ba.tag[0];
    tag[2]=ba.tag[1];
    tag[1]=ba.tag[2];
  }

  //Implementation class details
  class bamrecordImpl
  {
  private:
    template<typename T> T toItype(const char * beg)
    {
      return *(const T*)(beg);
    }
    void parse_block();
  public:
    bool __empty;
    I32 __block_size;
    std::unique_ptr<char[]> __block;
    I32 __refID,__pos;
    U32 __bin_mq_nl,__bin,__mq,__nl,__flag_nc,__flag,__nc;
    I32 __l_seq,__next_refID,__next_pos,__tlen;
    const char * __rname_beg, * __rname_end;
    const U32 * __cig_beg, * __cig_end;
    const U8 * __seq_beg, * __seq_end;
    const char * __qual_beg, * __qual_end, * __aux_beg, *__aux_end;
    bamrecordImpl(  );
    bamrecordImpl( std::int32_t blocksize,
		   std::unique_ptr<char[]> && block);
    bamrecordImpl( const bamrecordImpl & );
    bamrecordImpl& operator=( const bamrecordImpl &);
  };

  bamrecordImpl::bamrecordImpl() : __empty(true),
				   __block_size(0),
				   __block(nullptr),
				   __refID(-1),__pos(-1),
				   __bin_mq_nl(0),__flag_nc(0),
				   __l_seq(0),__next_refID(0),
				   __next_pos(0),
				   __tlen(0)
  {
  }
  
  bamrecordImpl::bamrecordImpl( const bamrecordImpl & bI) {
    *this = bI;
  }

  bamrecordImpl& bamrecordImpl::operator=( const bamrecordImpl & bI )  
  {
    this->__empty = bI.__empty;
    this->__block_size = bI.__block_size;
    if(this->__block_size)
      {
	this->__block = std::unique_ptr<char[]>(new char[__block_size]);
	std::copy(bI.__block.get(),bI.__block.get()+__block_size,
		  __block.get());
      }
    else this->__block == nullptr;
    this->parse_block();
    return *this;
  }

  void bamrecordImpl::parse_block() 
  {
    __refID= std::move(toItype<I32>(__block.get()));
    __pos = std::move(toItype<I32>(__block.get()+sizeof(I32)));
    __bin_mq_nl = std::move(toItype<U32>(__block.get()+2*sizeof(I32)));
    __bin = __bin_mq_nl>>16;
    __mq = ((__bin_mq_nl ^ __bin<<16)>>8 );
    __nl= __bin_mq_nl ^ __bin<<16 ^ __mq << 8 ;
    __flag_nc = std::move(toItype<U32>(__block.get()+2*sizeof(I32)+sizeof(U32)));
    __flag =  __flag_nc >> 16;
    __nc = __flag_nc ^ __flag << 16;
    __l_seq = std::move(toItype<I32>(__block.get()+2*(sizeof(I32)+sizeof(U32))));
    __next_refID = std::move(toItype<I32>(__block.get()+3*sizeof(I32)+2*sizeof(U32)));
    __next_pos = std::move(toItype<I32>(__block.get()+4*sizeof(I32)+2*sizeof(U32)));
    __tlen = std::move(toItype<I32>(__block.get()+5*sizeof(I32)+2*sizeof(U32)));
    __rname_beg = __block.get()+6*sizeof(I32)+2*sizeof(U32);
    __rname_end  = __block.get()+6*sizeof(I32)+2*sizeof(U32)+__nl;
    __cig_beg = (U32*)(__block.get()+6*sizeof(I32)+2*sizeof(U32)+__nl);
    __cig_end = (U32*)(__block.get()+6*sizeof(I32)+2*sizeof(U32)+__nl + __nc*sizeof(U32));
    __seq_beg = (U8*)(__block.get()+6*sizeof(I32)+2*sizeof(U32)+__nl + __nc*sizeof(U32));
    __seq_end = (U8*)(__block.get()+6*sizeof(I32)+2*sizeof(U32)+__nl + __nc*sizeof(U32) + (__l_seq+1)/2);
    __qual_beg =__block.get()+6*sizeof(I32)+2*sizeof(U32)+__nl + __nc*sizeof(U32) + (__l_seq+1)/2;
    __qual_end = __block.get()+6*sizeof(I32)+2*sizeof(U32)+__nl + __nc*sizeof(U32) + (__l_seq+1)/2 + __l_seq;
    __aux_beg = __block.get()+6*sizeof(I32)+2*sizeof(U32)+__nl + __nc*sizeof(U32) + (__l_seq+1)/2 + __l_seq;
    __aux_end = __block.get()+__block_size;
  }

  bamrecordImpl::bamrecordImpl( std::int32_t blocksize,
				std::unique_ptr<char[]> && block) : __empty(blocksize>0 ? false : true),
								    __block_size(std::move(blocksize)),
								    __block(std::move(block))															  
  {
    parse_block();
  }

  bamrecord::bamrecord(  ) : __impl( std::unique_ptr<bamrecordImpl>(new bamrecordImpl()) ) {}

  bamrecord::bamrecord( std::int32_t blocksize,
			std::unique_ptr<char[]> && block) :
    __impl(std::unique_ptr<bamrecordImpl>(new bamrecordImpl(blocksize,std::move(block)))) 
  {
  }

  bamrecord::bamrecord( const bamrecord & rhs) : __impl(std::unique_ptr<bamrecordImpl>(new bamrecordImpl(*rhs.__impl))) {}

  bool bamrecord::empty() const 
  {
    return __impl->__empty;
  }

  bamrecord::bamrecord(bamrecord&&rhs) : __impl(nullptr)
  {
    std::swap(this->__impl,rhs.__impl);
  }

  bamrecord & bamrecord::operator=(bamrecord&&rhs) 
  {
    bamrecordImpl  * t = this->__impl.release();
    delete t;
    std::swap(this->__impl,rhs.__impl);
    return *this;
  }

  bamrecord & bamrecord::operator=(const bamrecord &rhs )
  {
    bamrecordImpl * t = this->__impl.release();
    delete t;
    this->__impl = std::unique_ptr<bamrecordImpl>(new bamrecordImpl(*rhs.__impl));
    return *this;
  }
  
  bamrecord::~bamrecord() {}
  
  std::string bamrecord::read_name() const
  {
    return std::string(__impl->__rname_beg,__impl->__rname_end);
  }

  std::string bamrecord::seq() const 
  {
    std::string rv;
    std::for_each(this->seq_cbeg(),
     		  this->seq_cend(),
     		  [&](const U8 & i) {
     		    rv += int2seq[ (((i) >> 4) & 0x0F) ];
     		    rv += int2seq[ (((i)) & 0x0F) ];
     		  });
    return rv;
  }

  const std::uint8_t * bamrecord::seq_cbeg() const
  {
    return __impl->__seq_beg;
  }

  const std::uint8_t * bamrecord::seq_cend() const
  {
    return __impl->__seq_end;
  }

  std::string bamrecord::cigar() const 
  {
    std::string rv;
    std::for_each( __impl->__cig_beg,
     		   __impl->__cig_end,
     		   [&](const U32 & i )
     		   {
     		     U32 __op = (i>>4);
     		     rv += to_string(__op) + bamCig[i^(__op<<4)];
     		   } );
    return rv;
  }

  std::string bamrecord::qual() const 
  {
    return std::string(this->qual_cbeg(),
     		       this->qual_cend());
  }

  const char * bamrecord::qual_cbeg() const 
  {
    return __impl->__qual_beg;
  }

  const char * bamrecord::qual_cend() const 
  {
    return __impl->__qual_end;
  }

  std::pair< std::int32_t, const char * >
  bamrecord::raw() const
  {
    return std::make_pair(__impl->__block_size,__impl->__block.get());
  }

  samflag bamrecord::flag() const
  {
    return samflag(I32(__impl->__flag));
  }
 
  std::uint32_t bamrecord::mapq() const
  {
    return __impl->__mq;
  }

  std::int32_t bamrecord::pos() const {
    return __impl->__pos;
  }

  std::int32_t bamrecord::next_pos() const {
    return __impl->__next_pos;
  }

  std::int32_t bamrecord::refid() const {
    return __impl->__refID;
  }

  std::int32_t bamrecord::next_refid() const {
    return __impl->__next_refID;
  }

  std::int32_t bamrecord::tlen() const {
    return __impl->__tlen;
  }

  std::int32_t bamrecord::l_seq() const {
    return __impl->__l_seq;
  }
  
  const char * bamrecord::hasTag(const char * tag) const {
    if(tag==NULL||tag==nullptr || (__impl->__aux_beg == __impl->__aux_end) ) {
      return nullptr;
    }
    const char * rv = __impl->__aux_beg;
    
    for( ; rv < __impl->__aux_end-1 ; rv++)
      {
	if( *rv == tag[0] &&
	    *(rv+1) == tag[1]) 
	  {
	    if( std::find(&ValTypes[0],&ValTypes[0]+nValTypes,*(rv+2)) !=
		&ValTypes[0]+nValTypes )
	      return rv;
	  }
      }
    return nullptr;
  }
  
  bamaux bamrecord::aux(const char * tag) const {
    const char * tagspot = this->hasTag(tag);
    if( tagspot == nullptr ) return bamaux();
    char __tag[3];
    __tag[0]=*(tagspot++);
    __tag[1]=*(tagspot++);
    __tag[2] = '\0';
    char __val_type = char(*tagspot++);
    size_t valsize = sizeofTag(__val_type,tagspot,__impl->__aux_end);
    if(!valsize) return bamaux();
    std::unique_ptr<char[]> value(new char[valsize]);
    std::copy( tagspot,tagspot+valsize,value.get() );
    return bamaux(valsize,__tag,__val_type,value);
  }
    
  std::string bamrecord::allaux() const
  {
    if(__impl->__aux_beg == __impl->__aux_end) return std::string();
    string rv;
    const char * start = __impl->__aux_beg, * end = __impl->__aux_end;
    char tag[3];
    tag[2]='\0';
    while( start < end )
      {
	tag[0] = (*start++);
	tag[1] = (*start++);
	char val_type = (*start++);
	rv += string(tag) + ":" + val_type + ":";
	if ( val_type == 'A' )
	  {
	    rv += (*start++);
	  }
	else if ( val_type == 'c' )
	  {
	    I32 v = *(const I8*)(start++);
	    rv += to_string(v);
	  }
	else if (val_type == 'C')
	  {
	    I32 v = *(const U8*)(start++);
	    rv += to_string(v);
	  }
	else if (val_type == 's')
	  {
	    I32 v = *(const I16*)(start++);
	    rv += to_string(v);
	  }
	else if (val_type == 'S')
	  {
	    I32 v = *(const U16*)(start++);
	    rv += to_string(v);
	  }
	else if (val_type == 'i')
	  {
	    I32 v = *(const I32*)(start++);
	    rv += to_string(v);
	  }
	else if (val_type == 'I')
	  {
	    U32 v = *(const U32*)(start++);
	    rv += to_string(v);
	  }
	else if (val_type == 'B')//EXPERIMENTAL
	  {
	    char Btype = (*start++);
	    I32 Bsize = *(const I32*)(start++);
	    for( auto i = 0 ; i < Bsize ; ++i )
	      {
		if ( Btype == 'c' )
		  {
		    I32 v = *(const I8*)(start++);
		    rv += to_string(v);
		  }
		else if (Btype == 'C')
		  {
		    I32 v = *(const U8*)(start++);
		    rv += to_string(v);
		  }
		else if (Btype == 's')
		  {
		    I32 v = *(const I16*)(start++);
		    rv += to_string(v);
		  }
		else if (Btype == 'S')
		  {
		    I32 v = *(const U16*)(start++);
		    rv += to_string(v);
		  }
		else if (Btype == 'i')
		  {
		    I32 v = *(const I32*)(start++);
		    rv += to_string(v);
		  }
		else if (Btype == 'I')
		  {
		    U32 v = *(const U32*)(start++);
		    rv += to_string(v);
		  }
		if(i<Bsize-1)
		  rv += ',';
	      }
	  }
	else if (val_type == 'Z')
	  {
	    const char * __x = std::find(start,__impl->__aux_end,'\0');
	    rv += string( start,__x );
	    start += (__x - start)+1;
	  }
	rv += '\t';
      }
    return rv;
  }
}

#endif
