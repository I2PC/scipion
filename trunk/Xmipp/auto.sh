gettextize --force --copy
libtoolize --force --copy
aclocal -I m4
automake --gnu --add-missing
autoconf
