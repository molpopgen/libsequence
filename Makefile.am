AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
LIBTOOL_DEPS = @LIBTOOL_DEPS@	
SUBDIRS=src Sequence
includedir=$(oldincludedir) -ISequence
EXTRA_DIST=doc/libsequence.doxygen doc/Makefile doc/images/* examples/Makefile examples/*.cc examples/*.h examples/*.hpp examples/*.c README.md config-h.in

libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool
