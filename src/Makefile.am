include $(top_srcdir)/config/tophat.mk

AM_CPPFLAGS = -I$(top_srcdir)/SeqAn-1.1

#-- progs to be installed in $prefix/bin
bin_PROGRAMS = \
	cvg_islands \
	extract_reads \
	polyA_reads \
	read_ids \
	spanning_reads \
	tophat_reports \
	juncs_db 

#-- scripts to be installed in $prefix/bin
dist_bin_SCRIPTS = \
	contig_to_chr_coords.py \
	junctions.py \
	fq_all2std.pl \
	tophat.sh

#-- tophat library for linking convienence
noinst_LIBRARIES = libtophat.a
noinst_HEADERS = \
	const.h \
	islands.h \
	reads.h \
	common.h \
	bwt_map.h \
	wiggles.h \
	junctions.h \
	assert_helpers.h \
	alphabet.h \
	timer.h 
libtophat_a_SOURCES = \
	reads.cpp \
	alphabet.c \
	bwt_map.cpp \
	common.cpp \
	wiggles.cpp \
	junctions.cpp

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

juncs_db_SOURCES = juncs_db.cpp
juncs_db_LDADD = $(top_builddir)/src/libtophat.a

tophat_reports_SOURCES = tophat_reports.cpp
tophat_reports_LDADD = $(top_builddir)/src/libtophat.a