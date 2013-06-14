#!/bin/bash

#simple script to pack up a precompiled binary package
if [[ -z "$1" ]]; then
 echo -e "Usage:\n./make_bin.sh <binary_package_name> [<bam_prefix> <boost_prefix>]"
 exit 1
fi
echo "packing up $1.tar.gz, using BAM installation in $2, BOOST installation in $3"
/bin/rm -rf $1 $1.tar.gz
mkdir $1
make clean
make distclean
if [[ $(uname -m) = "x86_64" ]]; then
 echo "Linking statically on x86_64.."
 export LDFLAGS="-static-libgcc -static-libstdc++"
fi
if [[ $(uname) = "Darwin" ]]; then
 export CFLAGS="-mmacosx-version-min=10.6"
fi


l2="$2"
l3="$3"
if [[ -z "$l3" ]]; then
  l3="$l2"
fi

./configure --prefix=`pwd`/$1 --with-bam=$l2 --with-boost=$l3
sed -e 's|__PREFIX__||' src/tophat2.in > src/tophat2
make
make install
cp $1/bin/* $1

cp README $1
cp COPYING $1
cp AUTHORS $1

/bin/rm -rf $1/bin

tar cvfz $1.tar.gz $1
