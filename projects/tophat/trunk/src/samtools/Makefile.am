#include $(top_srcdir)/config/tophat.mk

#-- tophat library for linking convienence
noinst_LIBRARIES = libbam.a
noinst_HEADERS = \
	 bam.h \
	 bam_endian.h \
	 bam_maqcns.h \
	 bgzf.h \
	 faidx.h \
	 glf.h \
	 khash.h \
	 kseq.h \
	 ksort.h \
	 razf.h \
	 zutil.h 

AM_CPPFLAGS = -D_IOLIB=2 -D_FILE_OFFSET_BITS=64


libbam_a_SOURCES = \
	 bam.c \
	 bam_aux.c \
	 bam_import.c \
	 bam_lpileup.c \
	 bam_maqcns.c \
	 bam_mate.c \
	 bam_pileup.c \
	 bam_plcmd.c \
	 bam_sort.c \
	 bamtk.c \
	 bgzf.c \
	 bgzip.c \
	 faidx.c \
	 glf.c \
	 razf.c \
	 razip.c 
	
