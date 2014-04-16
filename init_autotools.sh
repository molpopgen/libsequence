#!sh
libtoolize --force --copy
autoreconf -fi
autoheader
automake --add-missing --copy