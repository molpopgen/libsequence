#!sh
libtoolize --force --copy
autoreconf -fis
autoheader
automake --add-missing --copy