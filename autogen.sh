#!/bin/sh

echo "Bootstrapping libbpm"
echo "Generate configure scripts and Makefiles..."

echo " -- running aclocal..."
aclocal  

echo " -- running autoconf, creating configure scripts..."
autoconf --force

# libtoolize is called glibtoolize on osx.. check this !
echo " -- library has libtool support, libtoolizing..."
if [ `uname` == "Darwin" ]; then
echo "  + You seem to be using Darwin, using glibtoolize !"
  glibtoolize --copy --force 
else
  libtoolize --copy --force
fi

echo " -- running automake, creating Makefile template..."
automake --gnu --copy --add-missing

echo "All done :)"

