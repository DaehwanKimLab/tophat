include $(top_srcdir)/config/tophat.mk

#-- progs to be installed in $prefix/bin
bin_PROGRAMS = \
	cvg_islands \
	extract_reads \
	polyA_reads \
	read_ids \
	spanning_reads

#-- scripts to be installed in $prefix/bin
dist_bin_SCRIPTS = \
	contig_to_chr_coords.py \
	junctions.py \
	tophat.sh

#-- tophat library for linking convienence
noinst_LIBRARIES = libtophat.a
noinst_HEADERS = \
	const.h \
	islands.h \
	reads.h \
	timer.h \
	tokenize.h
libtophat_a_SOURCES = \
	reads.cpp \
	tokenize.cpp

#-- program sources
cvg_islands_SOURCES = cvg_islands.cpp
cvg_islands_LDADD = -lz

extract_reads_SOURCES = extract_reads.cpp
extract_reads_LDADD = $(top_builddir)/src/libtophat.a

polyA_reads_SOURCES = polyA_reads.cpp
polyA_reads_LDADD = $(top_builddir)/src/libtophat.a

read_ids_SOURCES = read_ids.cpp
read_ids_LDADD = $(top_builddir)/src/libtophat.a

spanning_reads_SOURCES = spanning_reads.cpp
spanning_reads_LDADD = -lz $(top_builddir)/src/libtophat.a

