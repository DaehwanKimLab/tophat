#!/bin/bash

#simple script to pack up a precompiled binary package
if [[ -z "$1" ]]; then
 echo -e "Usage:\n./make_bin.sh <binary_package_name> [<bam_prefix> <boost_prefix>]"
 exit 1
fi
echo "packing up $1.tar.gz, using BAM installation in $2, BOOST installation in $3"
mkdir $1
make clean
./configure --prefix=`pwd`/$1 --with-bam=$2 --with-boost=$3
sed -e 's|__PREFIX__||' src/tophat2.in > src/tophat2
make
make install
cp $1/bin/* $1

cp README $1
cp COPYING $1
cp AUTHORS $1
cp -r annotation $1

rm -r $1/bin

tar cvfz $1.tar.gz $1
