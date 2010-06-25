#!/bin/bash

#simple script to pack up a precompiled binary package

echo "packing up $1.tar.gz"
mkdir $1
make clean
./configure --enable-intel64 --prefix=`pwd`/$1
make
make install
cp $1/bin/* $1

cp README $1
cp COPYING $1
cp AUTHORS $1

rm -r $1/bin

tar cvfz $1.tar.gz $1