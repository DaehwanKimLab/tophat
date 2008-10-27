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

LC_ALL=C

while getopts  "a:s:m:fq" flag
do
  #echo "$flag" $OPTIND $OPTARG
  case "$flag" in
        a) ANCHOR_LEN=$OPTARG; echo "anchor length set to $ANCHOR_LEN";;
        m) SPAN_MISMATCHES=$OPTARG; echo "span mismatches set to $SPAN_MISMATCHES";;
		s) SEED_SIZE=$OPTARG; echo "seed size set to $SEED_SIZE";;
		q) FASTQ_FORMAT=1; echo "Expecting reads in FASTQ format";;
        f) FASTA_FORMAT=1; echo "Expecting reads in FASTA format";;
  esac
done

if [ -z `which bowtie` ]; then
    echo "Error: Bowtie was not detected on this system.  Please verify that Bowtie is installed and that the Bowtie executables are in your path"
    exit
fi

if [ -z `which maq` ]; then
    echo "Error: Maq was not detected on this system.  Please verify that Maq is installed and that the Maq executables (maq, fq_all2std.pl) are in your path"
    exit
fi

if [ $ANCHOR_LEN -lt 3 ]; then
	echo "Error: anchor length must be at least 3."
	exit
fi

if [ $SPAN_MISMATCHES -gt 1 ]; then
	echo "Error: no more than one spanning mismatch currently supported."
	exit
fi

if  (($FASTA_FORMAT && $FASTQ_FORMAT)) ; then
	echo "Error: please supply only one of -q, -f"
	exit
fi

if !(($FASTA_FORMAT || $FASTQ_FORMAT)); then
	echo "Defaulting to FASTQ"
	FASTQ_FORMAT=1
fi

if (($FASTA_FORMAT)); then
	FORMAT="-f"
else
	FORMAT="-q"
fi

OPTIND=($OPTIND-1)

if [ -z "${ARGS[($OPTIND)]}" ]
then
	echo "Error: no Bowtie index specified"
	exit
fi 

if [ -z "${ARGS[($OPTIND+1)]}" ]
then
	echo "Error: no reads provided"
	exit
fi 

EBWT=${ARGS[$OPTIND]}
EBWT_SHORT=`echo $EBWT | awk -F/ '{print $NF}'`
#echo $EBWT_SHORT  

BWTMAP=reads_to_$EBWT_SHORT.bwtout
echo -n "Filtering out garbage reads in ${ARGS[($OPTIND+1)]} : "
date

cat ${ARGS[($OPTIND+1)]} | $BINDIR/polyA_reads $FORMAT > kept_reads


echo -n "Mapping reads in ${ARGS[($OPTIND+1)]} to $EBWT : "
date
echo "bowtie -l $SEED_SIZE $FORMAT $EBWT kept_reads > $BWTMAP"

bowtie -l $SEED_SIZE $FORMAT $EBWT kept_reads > $BWTMAP

MAQMAP=reads_to_$EBWT_SHORT.map

echo -n "Converting $BWTMAP to Maq format : "
date
bowtie-convert $BWTMAP $MAQMAP $EBWT.bfa 2>map_convert.log

echo -n "Collecting initially unmapped reads : "
date

LC_ALL=C

if  [[ "$FORMAT" -eq "-q" ]] ; then
	# Get all read ids
	$BINDIR/read_ids -q < kept_reads | sort > kept_read_ids

	# Get all mapped read ids
	awk '{print $1}' $BWTMAP | sort > mapped_read_ids
		
	# Record the ids of the unmapped reads
	comm -3 mapped_read_ids kept_read_ids > unmapped_read_ids
	
	# Grab the FASTQ files for the unmapped reads
	$BINDIR/extract_reads -q unmapped_read_ids < kept_reads > unmapped_reads.fq
	echo "Converting initially unmapped reads to FASTA format"
	fq_all2std.pl fq2fa unmapped_reads.fq > unmapped_reads.fa
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

echo -n "Assembling the consensus : "
date
echo "maq assemble -s $MAQCNS $EBWT.bfa $MAQMAP 2>/dev/null"
maq assemble -s $MAQCNS $EBWT.bfa $MAQMAP 2>/dev/null

echo -n "Extracting coverage islands from assembly : "
date
$BINDIR/cvg_islands -d 0.0 -b 6 -e 45 -R $MAQCNS islands.fa islands.gff

echo -n "Mapping initially unmapped reads against possible exon junctions : "
date
echo  "spanning_reads -v -a $ANCHOR_LEN -s $SEED_SIZE -m $SPAN_MISMATCHES -I $MAX_INTRON_LENGTH -i $MIN_INTRON_LENGTH -S 300 islands.fa islands.gff unmapped_reads.fa > reads_to_$EBWT_SHORT.splices"
$BINDIR/spanning_reads  -a $ANCHOR_LEN -s $SEED_SIZE -m $SPAN_MISMATCHES -I $MAX_INTRON_LENGTH -i $MIN_INTRON_LENGTH -S 300 islands.fa islands.gff unmapped_reads.fa > reads_to_$EBWT_SHORT.splices

# Convert .splices file into a UCSC genome browser BED file of unique junctions
$BINDIR/junctions reads_to_$EBWT_SHORT.splices > junctions.bed
