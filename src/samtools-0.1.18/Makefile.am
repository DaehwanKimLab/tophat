
AM_CFLAGS = -g -Wall -O2 #-m64 #-arch ppc
DFLAGS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -D_CURSES_LIB=0
ARFLAGS = -csru

#ldflags for samtools -lbcf -lm -lz -lcurses
#ldflags for razip & bgzip : -lz
#faidx_main.o:faidx.h razf.h <--- might need this

#libbam uses Lobjs
#samtools uses aobjs and links against libbam
#razip and bgzip use knetfile

#actual libraries
noinst_LIBRARIES = libbam.a

libbam_a_SOURCES = 	bgzf.h bgzf.c kstring.h kstring.c bam_aux.h bam_aux.c bam.c bam.h razf.h bam_endian.h \
			sam_header.h bam_import.c kseq.h khash.h sam.c sam.h bam_index.c ksort.h bam_pileup.c \
			bam_lpileup.c bam_md.c faidx.h razf.c faidx.c bedidx.h bedidx.c knetfile.h knetfile.c \
			bam_sort.h bam_sort.c sam_header.c bam_reheader.h bam_reheader.c kprobaln.h kprobaln.c \
			bam_cat.h bam_cat.c

bin_PROGRAMS = samtools_0.1.18
EXTRA_PROGRAMS = razip bgzip

samtools_0_1_18_SOURCES = bam_tview.h bam_tview.c bam_plcmd.c bam.h faidx.h bcftools/bcf.h bam2bcf.h sam_view.h \
			sam_view.c bam_rmdup.h bam_rmdup.c bam_rmdupse.h bam_rmdupse.c bam_mate.h bam_mate.c \
			bam_stat.h bam_stat.c bam_color.h bam_color.c bamtk.c kaln.h kaln.c bam2bcf.c errmod.h \
			bam2bcf_indel.c errmod.c sample.h sample.c cut_target.h cut_target.c phase.c khash.h \
			ksort.h bam2depth.h bam2depth.c

samtools_0_1_18_LDADD = libbam.a bcftools/libbcf.a -lm -lz -lcurses

razip_SOURCES =  razip.c razf.c razf.h knetfile.h knetfile.c
razip_LDADD = -lz
bgzip_SOURCES =  bgzip.c bgzf.h bgzf.c knetfile.h knetfile.c
bgzip_LDADD = -lz


SUBDIRS= bcftools .
