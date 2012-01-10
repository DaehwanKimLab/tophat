#!/bin/bash

# simple script to pack up source

echo "packing up $1.tar.gz"
mkdir $1

cp -r --preserve=all *.m4  build-aux src config.h.in configure.ac Makefile.am Makefile.in configure ChangeLog NEWS README COPYING AUTHORS INSTALL THANKS annotation $1
(cd $1; find -name .svn -exec rm -rf {} \;)
tar cvzf $1.tar.gz $1
rm -rf $1
