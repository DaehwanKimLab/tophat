#!/bin/sh

# TopHat
#
# Created by Cole Trapnell on 9/3/08.
# Copyright 2008 Cole Trapnell. All rights reserved.
BINDIR=.

ARGS=("$@")

FASTA_FORMAT=0
FASTQ_FORMAT=0
ANCHOR_LEN=5
SPAN_MISMATCHES=0
MAX_INTRON_LENGTH=20000
MIN_INTRON_LENGTH=70
SEED_SIZE=28
SELF_ISLAND_JUNCS=300
SOLEXA_SCALE=""
PHRED_CHAR="-R"

LC_ALL=C

USAGE_MSG="Usage: tophat [options] <index prefix> <reads>\n
Options:\n
\t-h\t\tPrint this message\n
\t-q\t\tExpect reads in FASTQ format [default: on]\n
\t-f\t\tExpect reads in FASTA format [default: off]\n
\t-a\t<3-6>\tanchor length [default: 5]\n
\t-m\t<0-2>\tmismatches allowed in extension [default: 0]\n
\t-s\t<int>\tseed length [default: 28]\n
\t-I\t<int>\tMaximum intron length [default: 20000]\n
\t-i\t<int>\tMinimum intron length [default: 70]\n
\t-X\t\tUse the solexa scale for qualities [default: off, use Phred scale]\n\n
Experimental features (USE WITH CAUTION):\n
\t-D\t<0-1000>\tMinimum normalized depth to look for junctions within single islands [default: 300]\n
\t-Q\t<';'-'I'>\t Allow Maq to call SNPs in islands sequences at positions >= this Phred qual [default: infinity, always use reference sequence]\n"

while getopts  "a:s:m:fqI:i:D:hX" flag
do
  #echo "$flag" $OPTIND $OPTARG
  case "$flag" in
		h) echo -e $USAGE_MSG; exit;;
        a) ANCHOR_LEN=$OPTARG; echo "anchor length set to $ANCHOR_LEN";;
        m) SPAN_MISMATCHES=$OPTARG; echo "span mismatches set to $SPAN_MISMATCHES";;
		s) SEED_SIZE=$OPTARG; echo "seed size set to $SEED_SIZE";;
		I) MAX_INTRON_LENGTH=$OPTARG; echo "Maximum intron length set to $MAX_INTRON_LENGTH";;
		i) MIN_INTRON_LENGTH=$OPTARG; echo "Minimum intron length set to $MIN_INTRON_LENGTH";;
		q) FASTQ_FORMAT=1; echo "Expecting reads in FASTQ format";;
        f) FASTA_FORMAT=1; echo "Expecting reads in FASTA format";;
		D) SELF_ISLAND_JUNCS=$OPTARG; echo "Will look for junctions within islands with D >= $SELF_ISLAND_JUNCS";;
		X) SOLEXA_SCALE="--solexa-quals"; echo "Qualities expected to be on the Solexa scale";;
		Q) PHRED_CHAR="-Q$OPTARG"; echo "Allowing SNPs at island consensus positions with $OPTARG Phred score or better";;
  esac
done

#ECHO=/bin/echo
ECHO="echo -e"

if [ -z `which bowtie` ]; then
    $ECHO "Error: Bowtie was not detected on this system.  Please verify that Bowtie is installed and that the Bowtie executables are in your path"
    exit
fi

# bowtie-convert was renamed in 0.9.8, this is for backwards-compatibility
if [ -z `which bowtie-maqconvert` ]; then
	BOWTIE_CONVERT=bowtie-maqconvert
else
	BOWTIE_CONVERT=bowtie-convert
fi

if [ -z `which maq` ]; then
    $ECHO "Error: Maq was not detected on this system.  Please verify that Maq is installed and that the Maq executables (maq, fq_all2std.pl) are in your path"
    exit
fi

if [ $ANCHOR_LEN -lt 3 ] || [ $ANCHOR_LEN -gt 6 ]; then
	$ECHO "Error: anchor length must be greater than 2 and less than 8."
	exit
fi

if [ $SEED_SIZE -lt 20 ]; then
	$ECHO "Error: seed length must be at least 20"
	exit
fi

if [ $SPAN_MISMATCHES -gt 2 ]; then
	$ECHO "Error: no more than two spanning mismatches currently supported."
	exit
fi

if  (($FASTA_FORMAT && $FASTQ_FORMAT)) ; then
	$ECHO "Error: please supply only one of -q, -f"
	exit
fi

if !(($FASTA_FORMAT || $FASTQ_FORMAT)); then
	$ECHO "Defaulting to FASTQ"
	FASTQ_FORMAT=1
fi

if (($FASTA_FORMAT)); then
	FORMAT="-f"
else
	FORMAT="-q"
fi

OPTIND=($OPTIND-1)

INDEX_FWD_1="${ARGS[($OPTIND)]}.1.ebwt"
INDEX_FWD_2="${ARGS[($OPTIND)]}.2.ebwt"
INDEX_REV_1="${ARGS[($OPTIND)]}.rev.1.ebwt"
INDEX_REV_2="${ARGS[($OPTIND)]}.rev.2.ebwt"
if ! [ -f $INDEX_FWD_1  ]; then
	$ECHO "Error: Bowtie index file $INDEX_FWD_1 not found"
	exit
fi 

if ! [ -f $INDEX_FWD_2  ]; then
	$ECHO "Error: Bowtie index file $INDEX_FWD_2 not found"
	exit
fi

if ! [ -f $INDEX_REV_1  ]; then
	$ECHO "Error: Bowtie index file $INDEX_REV_1 not found"
	exit
fi

if ! [ -f $INDEX_REV_2  ]; then
	$ECHO "Error: Bowtie index file $INDEX_REV_2 not found"
	exit
fi

INDEX_BFA="${ARGS[($OPTIND)]}.bfa"

if ! [ -f $INDEX_BFA  ]; then
	$ECHO "Error: Binary fasta $INDEX_BFA not found"
	exit
fi

if ! [ -f "${ARGS[($OPTIND)+1]}" ]; then
	echo "Error: no reads provided"
	exit
fi 

EBWT=${ARGS[$OPTIND]}
EBWT_SHORT=`echo $EBWT | awk -F/ '{print $NF}'`
#echo $EBWT_SHORT  

BWTMAP=reads_to_$EBWT_SHORT.bwtout
$ECHO  "Filtering out garbage reads in ${ARGS[($OPTIND+1)]} : \c"
date

cat ${ARGS[($OPTIND+1)]} | $BINDIR/polyA_reads $FORMAT > kept_reads


$ECHO  "Mapping reads in ${ARGS[($OPTIND+1)]} to $EBWT : \c"
date
echo "   bowtie -l $SEED_SIZE $FORMAT $SOLEXA_SCALE $EBWT kept_reads > $BWTMAP"

bowtie -l $SEED_SIZE $FORMAT $SOLEXA_SCALE $EBWT kept_reads > $BWTMAP

MAQMAP=reads_to_$EBWT_SHORT.map

$ECHO  "Converting $BWTMAP to Maq format : \c"
date
MAQ_VERSION_STR=`maq 2>&1 | grep Version | awk '{print $2}'`
MAQ_MAJOR_VERSION=`echo $MAQ_VERSION_STR | awk -F"." '{print $1}'`
MAQ_MINOR_VERSION=`echo $MAQ_VERSION_STR | awk -F"." '{print $2}'`
MAQ_FIX_VERSION=`echo $MAQ_VERSION_STR | awk -F"." '{print $3}'`
#echo $MAQ_MAJOR_VERSION
#echo $MAQ_MINOR_VERSION
#echo $MAQ_FIX_VERSION
if [ $MAQ_MINOR_VERSION -lt 7 ] ; then 
	$ECHO "Detected Maq older than 0.7.0, using 64bp read mapping format"
	bowtie-convert -o  $BWTMAP $MAQMAP $EBWT.bfa 2>map_convert.log
else
	bowtie-convert $BWTMAP $MAQMAP $EBWT.bfa 2>map_convert.log
fi

$ECHO  "Collecting initially unmapped reads : \c"
date

LC_ALL=C

if  [[ "$FORMAT" == "-q" ]] ; then
	# Get all read ids
	$BINDIR/read_ids -q < kept_reads | sort > kept_read_ids

	# Get all mapped read ids
	awk '{print $1}' $BWTMAP | sort > mapped_read_ids
		
	# Record the ids of the unmapped reads
	comm -3 mapped_read_ids kept_read_ids > unmapped_read_ids
	
	# Grab the FASTQ files for the unmapped reads
	$BINDIR/extract_reads -q unmapped_read_ids < kept_reads > unmapped_reads.fq
	$ECHO "Converting initially unmapped reads to FASTA format"
	$BINDIR/fq_all2std fq2fa unmapped_reads.fq > unmapped_reads.fa
	rm unmapped_reads.fq
else
	# Get all read ids
	$BINDIR/read_ids -f < kept_reads | sort > kept_read_ids

	# Get all mapped read ids
	awk '{print $1}' $BWTMAP | sort > mapped_read_ids

	comm -3 mapped_read_ids kept_read_ids > unmapped_read_ids

	$BINDIR/extract_reads -f unmapped_read_ids < kept_reads > unmapped_reads.fa
fi

MAQCNS=reads_to_$EBWT_SHORT.cns

$ECHO  "Assembling the consensus : \c"
date
$ECHO "   maq assemble -s $MAQCNS $EBWT.bfa $MAQMAP 2>/dev/null"
maq assemble -s $MAQCNS $EBWT.bfa $MAQMAP 2>/dev/null

$ECHO  "Extracting coverage islands from assembly : \c"
date
$ECHO  "   cvg_islands -d 0.0 -b 6 -e 45 $PHRED_CHAR $MAQCNS islands.fa islands.gff"

$BINDIR/cvg_islands -d 0.0 -b 6 -e 45 $PHRED_CHAR $MAQCNS islands.fa islands.gff

$ECHO  "Mapping initially unmapped reads against possible exon junctions : \c"
date
$ECHO  "   spanning_reads -v -a $ANCHOR_LEN -s $SEED_SIZE -m $SPAN_MISMATCHES -I $MAX_INTRON_LENGTH -i $MIN_INTRON_LENGTH -S 300 islands.fa islands.gff unmapped_reads.fa > reads_to_$EBWT_SHORT.splices"
$BINDIR/spanning_reads  -a $ANCHOR_LEN -s $SEED_SIZE -m $SPAN_MISMATCHES -I $MAX_INTRON_LENGTH -i $MIN_INTRON_LENGTH -S $SELF_ISLAND_JUNCS islands.fa islands.gff unmapped_reads.fa > reads_to_$EBWT_SHORT.splices
$ECHO  "Collecting junctions from spliced reads : \c"
date

# Convert .splices file into a UCSC genome browser BED file of unique junctions
$BINDIR/junctions reads_to_$EBWT_SHORT.splices > junctions.bed
$ECHO "TopHat run complete"
date
