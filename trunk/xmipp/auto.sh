# Bootstrap script

aclocal -I m4
autoheader
libtoolize --force --copy
automake --gnu --add-missing
autoconf

./configure
#./configure --disable-static --enable-debug --enable-profiling
make
make dist

