#!/usr/bin/env bash

#simple script to pack up a precompiled binary package
if [[ -z "$1" ]]; then
 echo -e "Usage:\n./make_bin.sh <binary_package_name> [<boost_prefix>]"
 echo "  <binary_package_name> cannot be a full path, but only a file name"
 echo "  e.g.: tophat-2.1.0.Linux_x86_64"
 exit 1
fi
#echo "packing up $1.tar.gz, using BOOST installation in $2"
build="$1"
boostpre="$2"
if [[ -z "$boostpre" ]]; then
  echo "Preparing binary package $build.tar.gz using the system BOOST installation"
else
  echo "Preparing binary package $build.tar.gz using BOOST from $boostpre"
fi
/bin/rm -rf $build $build.tar.gz
mkdir $build
if [[ -f Makefile ]]; then
  make clean
  make distclean
fi
if [[ $(uname -m) = "x86_64" ]]; then
 echo "Linking statically on x86_64.."
 export LDFLAGS="-static-libgcc -static-libstdc++"
fi
if [[ $(uname) = "Darwin" ]]; then
 export CFLAGS="-mmacosx-version-min=10.7"
fi

./configure --prefix=`pwd`/$build --with-boost=$boostpre
make install
mv $build/bin/* $build

cp README $build
cp LICENSE $build
cp AUTHORS $build

/bin/rm -rf $build/bin

tar cvfz $build.tar.gz $build
