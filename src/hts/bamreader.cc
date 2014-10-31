#ifdef HAVE_HTSLIB //Will only compile if ./configure detects htslib

#include <Sequence/bamreader.hpp>
#include <htslib/bgzf.h>

using std::string;

namespace Sequence
{
  using I32 = std::int32_t;
  //!Impl class for bamreader
  class bamreaderImpl
  {
  public:
    using I32 = std::int32_t;
    BGZF * in;
    bool __EOF,__errorstate;
    char __magic[4];
    I32 __l_text,__n_ref;
    std::unique_ptr<char[]> __htext;
    std::vector< std::pair<std::string,I32> > __refdata;

    bamreaderImpl(const char * bamfilename);
    ~bamreaderImpl();
  }; 

  bamreaderImpl::~bamreaderImpl() 
  {
    if(in != NULL) bgzf_close(in);
  }

  bamreaderImpl::bamreaderImpl(const char * bamfilename) :
    in((bamfilename != nullptr) ? bgzf_open(bamfilename,"rb") : NULL),
    __EOF(false),
    __errorstate(false),
    __htext(nullptr),
    __refdata(std::vector< std::pair<std::string,I32> >())
  {
    if(gzopen != NULL)
      {
	auto rv = bgzf_read( in, &__magic[0], 4*sizeof(char) );
	if (!rv) __EOF = true;
	if(rv==-1) __errorstate = true;
	if(string({__magic[0],__magic[1],__magic[2]}) != string("BAM")) __errorstate = 1;
	if(!__errorstate && !__EOF)
	  {
	    rv = bgzf_read( in, &__l_text, sizeof(I32) );
	    if (!rv) {__EOF = true; return; }
	    if(rv==-1){ __errorstate = true; return; }
	    __htext = std::unique_ptr<char[]>( new char[__l_text] );
	    rv = bgzf_read(in,__htext.get(),size_t(__l_text)*sizeof(char));
	    if (!rv) {__EOF = true; return; }
	    if(rv==-1){ __errorstate = true; return; }
	    rv = bgzf_read(in,&__n_ref,sizeof(I32));
	    if (!rv) {__EOF = true; return; }
	    if(rv==-1){ __errorstate = true; return; }
	    for(decltype(__n_ref) i = 0 ; i < __n_ref ; ++i )
	      {
		I32 l_name,l_ref;
		rv = bgzf_read( in,&l_name,sizeof(I32) );
		if (!rv) {__EOF = true; return; }
		if(rv==-1){ __errorstate = true; return; }
		char name[l_name];
		rv = bgzf_read( in,&name[0],size_t(l_name)*sizeof(char));
		if (!rv) {__EOF = true; return; }
		if(rv==-1){ __errorstate = true; return; }
		rv = bgzf_read( in,&l_ref,sizeof(I32) );
		if (!rv) {__EOF = true; return; }
		if(rv==-1){ __errorstate = true; return; }
		__refdata.push_back(std::make_pair(std::string(name),l_ref));
	      }
	  }
      }
    else __errorstate = 1;
  }

  bamreader::bamreader( const char * bamfilename) :
    __impl( new bamreaderImpl(bamfilename) )
  {
  }

  bamreader::~bamreader(){}

  bamrecord bamreader::next_record() const
  {
    I32 bsize;
    auto rv = bgzf_read(__impl->in,&bsize,sizeof(I32));
    if(!rv) { __impl->__EOF=1; return bamrecord(); }
    if(rv==-1) { __impl->__errorstate = 1; return bamrecord(); }
    std::unique_ptr<char[]> block(new char[bsize]);
    rv = bgzf_read(__impl->in,&block[0],size_t(bsize));
    if(!rv) { __impl->__EOF=1; return bamrecord(); }
    if(rv==-1) { __impl->__errorstate = 1; return bamrecord(); }
    return bamrecord(bsize,std::move(block));
  }

  bamrecord bamreader::record_at_pos( std::int64_t offset ) const 
  {
    auto current = bgzf_tell(__impl->in);
    auto rv = bgzf_seek(__impl->in, offset, SEEK_SET );
    if(rv==-1) {__impl->__errorstate=true;return bamrecord();}

    I32 bsize;
    rv = bgzf_read(__impl->in,&bsize,sizeof(I32));
    if(rv==0) {__impl->__EOF=true;return bamrecord();}
    if(rv==-1) {__impl->__errorstate=true;return bamrecord();}
    std::unique_ptr<char[]> block(new char[bsize]);

    rv = bgzf_read(__impl->in,block.get(),size_t(bsize)*sizeof(char));
    if(rv==0) {__impl->__EOF=true;return bamrecord();}
    if(rv==-1) {__impl->__errorstate=true;return bamrecord();}
    bamrecord b(bsize,std::move(block));

    //restore offset
    rv = bgzf_seek(__impl->in, current, SEEK_SET );
    if(rv==-1) {__impl->__errorstate=true;return bamrecord();}
    return b;
  }

  bool bamreader::eof() const
  {
    return __impl->__EOF;
  }

  bool bamreader::has_eof() const 
  {
    return bgzf_check_EOF(__impl->in);
  }

  bool bamreader::error() const
  {
    return __impl->__errorstate;
  }

  std::int64_t bamreader::rewind() 
  {
    return bgzf_seek(__impl->in,0L,SEEK_SET);
  }

  std::int64_t bamreader::seek( std::int64_t offset, int whence )
  {
    return bgzf_seek(__impl->in,std::move(offset),std::move(whence));
  }

  int bamreader::close()
  {
    return bgzf_close(__impl->in);
  }

  std::int64_t bamreader::tell() 
  {
    return bgzf_tell(__impl->in);
  }

  bamreader::refdataObj bamreader::operator[](const size_type & i) {
    return __impl->__refdata[i];
  }

  bamreader::refdata_citr bamreader::ref_cbegin() const
  {
    return __impl->__refdata.cbegin();
  }

  bamreader::refdata_citr bamreader::ref_cend() const
  {
    return __impl->__refdata.cend();
  }

  std::string bamreader::header() const
  {
    return std::string(__impl->__htext.get());
  }

  std::int32_t bamreader::n_ref() const {
    return __impl->__n_ref;
  }

}

#endif
