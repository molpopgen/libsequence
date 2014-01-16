DOX=doxygen
all: 
	$(DOX) libsequence.doxygen
	cd latex && make
	cd ..
clean:
	rm -rf html latex
install:
	install -d html /usr/share/doc/libsequence/html
	install html/*  /usr/share/doc/libsequence/html
	install  latex/refman.pdf /usr/share/doc/libsequence/libsequence-ref.pdf