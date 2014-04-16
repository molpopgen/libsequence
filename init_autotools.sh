#!sh
libtoolize --force
autoreconf -fis
autoheader
automake --add-missing