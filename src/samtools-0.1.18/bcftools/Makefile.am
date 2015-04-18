AM_CFLAGS = -I$(top_srcdir)/src/samtools-0.1.18

lib_LIBRARIES = libbcf.a
#bcf.o vcf.o bcfutils.o prob1.o em.o kfunc.o kmin.o index.o fet.o mut.o bcf2qcall.o
libbcf_a_SOURCES = bcf.h kmin.h prob1.h \
		bcf.c kmin.c prob1.c \
		vcf.c bcfutils.c em.c kfunc.c index.c fet.c mut.c bcf2qcall.c

bin_PROGRAMS = bcftools
#call1.o main.o $(OMISC)/kstring.o $(OMISC)/bgzf.o $(OMISC)/knetfile.o $(OMISC)/bedidx.o
bcftools_SOURCES = main.c call1.c \
		../kstring.h ../bgzf.h ../knetfile.h \
		../kstring.c ../bgzf.c ../knetfile.c ../bedidx.c
bcftools_LDADD = libbcf.a -lm -lz
