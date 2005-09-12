# As developer
# Packages needed: gettext, gettext-devel, automake, autoconf

gettextize --force --copy
aclocal -I m4
autoheader
libtoolize --force --copy
automake --gnu --add-missing
autoconf

./configure
./configure --disable-static
make
make dist

# As user:
# 
# ./configure
# make
