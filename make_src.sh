#!/usr/bin/env bash

# simple script to pack up source

echo "packing up $1.tar.gz"
mkdir $1

make clean
cp -r --preserve=all *.m4 src config.h.in configure.ac Makefile.am Makefile.in configure install-sh config.sub config.guess depcomp ChangeLog NEWS README COPYING AUTHORS INSTALL THANKS $1
(cd $1; find -name .svn -exec rm -rf {} \;)
tar cvzf $1.tar.gz $1
rm -rf $1
