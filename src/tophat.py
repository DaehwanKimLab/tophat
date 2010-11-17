#!/usr/bin/env python

# encoding: utf-8
"""
tophat.py

Created by Cole Trapnell on 2008-12-25.
Copyright (c) 2008 Cole Trapnell. All rights reserved.
"""

import sys
try:
    import psyco
    psyco.full()
except ImportError:
    pass

import sys
import getopt
import subprocess
import errno
import os
import tempfile
import warnings
import shutil
import copy
from datetime import datetime, date, time

use_message = '''
TopHat maps short sequences from spliced transcripts to whole genomes.

Usage:
    tophat [options] <bowtie_index> <reads1[,reads2,...,readsN]> [reads1[,reads2,...,readsN]] [quals1,[quals2,...,qualsN]] [quals1[,quals2,...,qualsN]]
    
Options:
    -v/--version
    -o/--output-dir                <string>    [ default: ./tophat_out ]
    -a/--min-anchor                <int>       [ default: 8            ]
    -m/--splice-mismatches         <0-2>       [ default: 0            ]
    -i/--min-intron                <int>       [ default: 50           ]
    -I/--max-intron                <int>       [ default: 500000       ]
    -g/--max-multihits             <int>       [ default: 40           ]
    -F/--min-isoform-fraction      <float>     [ default: 0.15         ]
    --solexa-quals                          
    --solexa1.3-quals                          (same as phred64-quals)
    --phred64-quals                            (same as solexa1.3-quals)
    -Q/--quals
    --integer-quals
    -C/--color                                 (Solid - color space)
    --color-out
    --library-type                             (--fr-unstranded, --fr-firststrand, --fr-secondstrand, --ff-unstranded, --ff-firststrand, --ff-secondstrand)
    -p/--num-threads               <int>       [ default: 1            ]
    -G/--GTF                       <filename>
    -j/--raw-juncs                 <filename>
    -r/--mate-inner-dist           <int>       
    --mate-std-dev                 <int>       [ default: 20           ]
    --no-novel-juncs                           
    --no-gtf-juncs                             
    --no-coverage-search
    --coverage-search                                              
    --no-closure-search
    --closure-search
    --fill-gaps        
    --microexon-search
    --butterfly-search
    --no-butterfly-search
    --keep-tmp
    --tmp-dir                      <dirname>
    
Advanced Options:

    --segment-mismatches           <int>       [ default: 2            ]
    --segment-length               <int>       [ default: 25           ]
    --min-closure-exon             <int>       [ default: 100          ]
    --min-closure-intron           <int>       [ default: 50           ]
    --max-closure-intron           <int>       [ default: 5000         ]
    --min-coverage-intron          <int>       [ default: 50           ]
    --max-coverage-intron          <int>       [ default: 20000        ]
    --min-segment-intron           <int>       [ default: 50           ]
    --max-segment-intron           <int>       [ default: 500000       ] 

SAM Header Options (for embedding sequencing run metadata in output):
    --rg-id                        <string>    (read group ID)
    --rg-sample                    <string>    (sample ID)
    --rg-library                   <string>    (library ID)
    --rg-description               <string>    (descriptive string, no tabs allowed)
    --rg-platform-unit             <string>    (e.g Illumina lane ID)
    --rg-center                    <string>    (sequencing center name)
    --rg-date                      <string>    (ISO 8601 date of the sequencing run)
    --rg-platform                  <string>    (Sequencing platform descriptor)  
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

output_dir = "./tophat_out/"
logging_dir = output_dir + "logs/"
run_log = None
run_cmd = None
fasta_linebuf  = "" # buffer for fa_ungetline(), checked by fa_getline() before reading from file
fasta_lastline = "" # last line read by fa_getline()
tmp_dir = output_dir + "tmp/"
bin_dir = sys.path[0] + "/"
#ok_str = "\t\t\t\t[OK]\n"
fail_str = "\t[FAILED]\n"

sam_header = tmp_dir + "stub_header.sam"

# TopHatParams captures all of the runtime paramaters used by TopHat, and many
# of these are passed as command line options to exectubles run by the pipeline

# This class and its nested classes also do options parsing through parse_options()
# and option validation via the member function check()

class TopHatParams:

    # SpliceConstraints is a group of runtime parameters that specify what 
    # constraints to put on junctions discovered by the program.  These constraints
    # are used to filter out spurious/false positive junctions.
    
    class SpliceConstraints:
        def __init__(self, 
                     min_anchor_length,
                     min_intron_length,
                     max_intron_length, 
                     splice_mismatches,
                     min_isoform_fraction):
            self.min_anchor_length = min_anchor_length
            self.min_intron_length = min_intron_length
            self.max_intron_length = max_intron_length
            self.splice_mismatches = splice_mismatches
            self.min_isoform_fraction = min_isoform_fraction
        
        def parse_options(self, opts):
            for option, value in opts:
                if option in ("-m", "--splice-mismatches"):
                    self.splice_mismatches = int(value)
                elif option in ("-a", "--min-anchor"):
                    self.min_anchor_length = int(value)
                elif option in ("-F", "--min-isoform-fraction"):
                    self.min_isoform_fraction = float(value)
                elif option in ("-i", "--min-intron-length"):
                    self.min_intron_length = int(value)
                elif option in ("-I", "--max-intron-length"):
                    self.max_intron_length = int(value)
        
        def check(self):
            if self.splice_mismatches not in [0,1,2]:
                print >> sys.stderr, "Error: arg to --splice-mismatches must be 0, 1, or 2"
                sys.exit(1)
            if self.min_anchor_length < 4:
                print >> sys.stderr, "Error: arg to --min-anchor-len must be greater than 4"
                sys.exit(1)
            if self.min_isoform_fraction < 0.0 or self.min_isoform_fraction > 1.0:
                print >> sys.stderr, "Error: arg to --min-isoform-fraction must be between 0.0 and 1.0"
                sys.exit(1)
            if self.min_intron_length <= 0:
                print >> sys.stderr, "Error: arg to --min-intron-length must be greater than 0"
                sys.exit(1)                    
            if self.max_intron_length <= 0:
                print >> sys.stderr, "Error: arg to --max-intron-length must be greater than 0"
                sys.exit(1)
    
    # SystemParams is a group of runtime parameters that determine how to handle
    # temporary files produced during a run and how many threads to use for threaded
    # stages of the pipeline (e.g. Bowtie)
    
    class SystemParams:
        def __init__(self,
                     bowtie_threads,
                     keep_tmp):
            self.bowtie_threads = bowtie_threads
            self.keep_tmp = keep_tmp
            
        def parse_options(self, opts):
            for option, value in opts:
                if option in ("-p", "--num-threads"):
                    self.bowtie_threads = int(value)
                elif option == "--keep-tmp":
                    self.keep_tmp = True
        
        def check(self):
            pass
        
    
    # ReadParams is a group of runtime parameters that specify various properties
    # of the user's reads (e.g. which quality scale their are on, how long the 
    # fragments are, etc).
    
    class ReadParams:
        def __init__(self,
                     solexa_quals,
                     phred64_quals,
                     quals,
                     integer_quals,
                     color,
                     color_out,
                     library_type,
                     seed_length,
                     reads_format,
                     mate_inner_dist,
                     mate_inner_dist_std_dev,
                     read_group_id,
                     sample_id,
                     library_id,
                     description,
                     seq_platform_unit,
                     seq_center,
                     seq_run_date,
                     seq_platform):
            self.solexa_quals = solexa_quals
            self.phred64_quals = phred64_quals
            self.quals = quals
            self.integer_quals = integer_quals
            self.color = color
            self.color_out = color_out
            self.library_type = library_type
            self.seed_length = seed_length
            self.reads_format = reads_format
            self.mate_inner_dist = mate_inner_dist
            self.mate_inner_dist_std_dev = mate_inner_dist_std_dev
            self.read_group_id = read_group_id 
            self.sample_id = sample_id
            self.library_id = library_id
            self.description = description
            self.seq_platform_unit = seq_platform_unit
            self.seq_center = seq_center
            self.seq_run_date = seq_run_date
            self.seq_platform = seq_platform
            
        def parse_options(self, opts):
            for option, value in opts:
                if option == "--solexa-quals":
                    self.solexa_quals = True
                elif option in ("--solexa1.3-quals", "--phred64-quals"):
                    self.phred64_quals = True
                elif option in ("-Q", "--quals"):
                    self.quals = True
                elif option == "--integer-quals":
                    self.integer_quals = True
                elif option in ("-C", "--color"):
                    self.color = True
                elif option == "--color-out":
                    self.color_out = True
                elif option == "--library-type":
                    self.library_type = value
                elif option in ("-s", "--seed-length"):
                    self.seed_length = int(value)
                elif option in ("-r", "--mate-inner-dist"):
                    self.mate_inner_dist = int(value)
                elif option == "--mate-std-dev":
                    self.mate_inner_dist_std_dev = int(value)
                elif option == "--rg-id":
                    self.read_group_id = value
                elif option == "--rg-sample":
                    self.sample_id = value
                elif option == "--rg-library":
                    self.library_id = value
                elif option == "--rg-description":
                    self.description = value
                elif option == "--rg-platform-unit":
                    self.seq_platform_unit = value
                elif option == "--rg-center":
                    self.seq_center = value
                elif option == "--rg-date":
                    self.seq_run_date = value    
                elif option == "--rg-platform":
                    self.seq_platform = value            

        def check(self):
            if self.seed_length != None and self.seed_length < 20:
                print >> sys.stderr, "Error: arg to --seed-length must be at least 20"
                sys.exit(1)
            if self.mate_inner_dist_std_dev != None and self.mate_inner_dist_std_dev < 0:
                print >> sys.stderr, "Error: arg to --mate-std-dev must at least 0"
                sys.exit(1)
            if (not self.read_group_id and self.sample_id) or (self.read_group_id and not self.sample_id):
                print >> sys.stderr, "Error: --rg-id and --rg-sample must be specified or omitted together"
                sys.exit(1)
    
    # SearchParams is a group of runtime parameters that specify how TopHat will
    # search for splice junctions
    
    class SearchParams:
        def __init__(self,
                     min_closure_exon,
                     min_closure_intron,
                     max_closure_intron,
                     min_coverage_intron,
                     max_coverage_intron,
                     min_segment_intron,
                     max_segment_intron):
                     
             self.min_closure_exon_length = min_closure_exon
             self.min_closure_intron_length = min_closure_intron
             self.max_closure_intron_length = max_closure_intron
             self.min_coverage_intron_length = min_coverage_intron
             self.max_coverage_intron_length = max_coverage_intron
             self.min_segment_intron_length = min_segment_intron
             self.max_segment_intron_length = max_segment_intron

        def parse_options(self, opts):
            for option, value in opts:
                if option == "--min-closure-exon":
                    self.min_closure_exon_length = int(value)
                if option == "--min-closure-intron":
                    self.min_closure_intron_length = int(value)
                if option == "--max-closure-intron":
                    self.max_closure_intron_length = int(value)
                if option == "--min-coverage-intron":
                    self.min_coverage_intron_length = int(value)
                if option == "--max-coverage-intron":
                    self.max_coverage_intron_length = int(value)
                if option == "--min-segment-intron":
                    self.min_segment_intron_length = int(value)
                if option == "--max-segment-intron":
                    self.max_segment_intron_length = int(value)

        def check(self):
            if self.min_closure_exon_length < 0:
                print >> sys.stderr, "Error: arg to --min-closure-exon must be at least 20"
                sys.exit(1)
            if self.min_closure_intron_length < 0:
                print >> sys.stderr, "Error: arg to --min-closure-intron must be at least 20"
                sys.exit(1)
            if self.max_closure_intron_length < 0:
                print >> sys.stderr, "Error: arg to --max-closure-intron must be at least 20"
                sys.exit(1)
            if self.min_coverage_intron_length < 0:
                print >> sys.stderr, "Error: arg to --min-coverage-intron must be at least 20"
                sys.exit(1)
            if self.max_coverage_intron_length < 0:
                print >> sys.stderr, "Error: arg to --max-coverage-intron must be at least 20"
                sys.exit(1)
            if self.min_segment_intron_length < 0:
                print >> sys.stderr, "Error: arg to --min-segment-intron must be at least 20"
                sys.exit(1)
            if self.max_segment_intron_length < 0:
                print >> sys.stderr, "Error: arg to --max-segment-intron must be at least 20"
                sys.exit(1)
                    
    def __init__(self):        
        self.splice_constraints = self.SpliceConstraints(8,     # min_anchor 
                                                         50,    # min_intron
                                                         500000, # max_intron
                                                         0,     # splice_mismatches
                                                         0.15)  # min_isoform_frac
        
        self.read_params = self.ReadParams(False,               # solexa_scale
                                           False,
                                           False,               # quals
                                           None,                # integer quals
                                           False,               # SOLiD - color space
                                           False,               # SOLiD - color out instead of base pair,
                                           "",                  # library type (e.g. "illumina-stranded-pair-end")
                                           None,                # seed_length
                                           "fastq",             # quality_format
                                           None,                # mate inner distance
                                           20,                  # mate inner dist std dev
                                           None,                # read group id
                                           None,                # sample id
                                           None,                # library id
                                           None,                # description
                                           None,                # platform unit (i.e. lane)
                                           None,                # sequencing center
                                           None,                # run date
                                           None)                # sequencing platform
        
        self.system_params = self.SystemParams(1,               # bowtie_threads
                                               False)           # keep_tmp   
                                               
        self.search_params = self.SearchParams(100,             # min_closure_exon_length
                                               50,              # min_closure_intron_length
                                               5000,            # max_closure_intron_length
                                               50,              # min_coverage_intron_length
                                               20000,           # max_coverage_intron_length
                                               50,              # min_segment_intron_length
                                               500000)          # max_segment_intron_length                          
        
        self.gff_annotation = None
        self.raw_junctions = None
        self.find_novel_juncs = True
        self.find_GFF_juncs = True
        self.skip_check_reads = False
        self.max_hits = 40
        self.segment_length = 25
        self.segment_mismatches = 2
        self.closure_search = None
        self.coverage_search = None
        self.microexon_search = False
        self.butterfly_search = None
        
    def check(self):
        self.splice_constraints.check()
        self.read_params.check()
        self.system_params.check()
       
        if self.segment_length <= 4:
            print >> sys.stderr, "Error: arg to --segment-length must at least 4"
            sys.exit(1)
        if self.segment_mismatches < 0 or self.segment_mismatches > 3:
            print >> sys.stderr, "Error: arg to --segment-mismatches must in [0, 3]"
            sys.exit(1)

        if self.read_params.color == True and self.butterfly_search == True:
            print >> sys.stderr, "Error: butterfly-search in colorspace is not yet supported"
            sys.exit(1)

        library_types = ["fr-unstranded", "fr-firststrand", "fr-secondstrand"]

        if self.read_params.library_type != "" and self.read_params.library_type not in library_types:
            print >> sys.stderr, "Error: libary-type should be one of", library_types
            sys.exit(1)
        
        self.search_params.max_closure_intron_length = min(self.splice_constraints.max_intron_length,
                                                           self.search_params.max_closure_intron_length)
        
        self.search_params.max_segment_intron_length = min(self.splice_constraints.max_intron_length,
                                                           self.search_params.max_segment_intron_length)

        self.search_params.max_coverage_intron_length = min(self.splice_constraints.max_intron_length,
                                                            self.search_params.max_coverage_intron_length)
        
    def cmd(self):
        cmd = ["--min-anchor", str(self.splice_constraints.min_anchor_length),
               "--splice-mismatches", str(self.splice_constraints.splice_mismatches),
               "--min-report-intron", str(self.splice_constraints.min_intron_length),
               "--max-report-intron", str(self.splice_constraints.max_intron_length),
               "--min-isoform-fraction", str(self.splice_constraints.min_isoform_fraction),
               "--output-dir", output_dir,
               "--max-multihits", str(self.max_hits),
               "--segment-length", str(self.segment_length),
               "--segment-mismatches", str(self.segment_mismatches),
               "--min-closure-exon", str(self.search_params.min_closure_exon_length),
               "--min-closure-intron", str(self.search_params.min_closure_intron_length),
               "--max-closure-intron", str(self.search_params.max_closure_intron_length),
               "--min-coverage-intron", str(self.search_params.min_coverage_intron_length),
               "--max-coverage-intron", str(self.search_params.max_coverage_intron_length),
               "--min-segment-intron", str(self.search_params.min_segment_intron_length),
               "--max-segment-intron", str(self.search_params.max_segment_intron_length),
               "--sam-header", sam_header]
        
               
        if self.read_params.mate_inner_dist != None:
            cmd.extend(["--inner-dist-mean", str(self.read_params.mate_inner_dist),
                        "--inner-dist-std-dev", str(self.read_params.mate_inner_dist_std_dev)])
        if self.gff_annotation != None:
            cmd.extend(["--gtf-annotations", str(self.gff_annotation)])
        if self.closure_search == False:
            cmd.append("--no-closure-search")
        if self.coverage_search == False:
            cmd.append("--no-coverage-search")
        if self.microexon_search == False:
            cmd.append("--no-microexon-search")            
        if self.butterfly_search == True:
            cmd.append("--butterfly-search")
        if self.read_params.solexa_quals == True:
            cmd.append("--solexa-quals")
        if self.read_params.quals == True:
            cmd.append("--quals")
        if self.read_params.integer_quals == True:
            cmd.append("--integer-quals")
        if self.read_params.color == True:
            cmd.append("--color")
            if self.read_params.color_out == True:
                cmd.append("--color-out")
        if self.read_params.library_type != "":
            cmd.extend(["--library-type", self.read_params.library_type])
        if self.read_params.phred64_quals == True:
            cmd.append("--phred64-quals")
        return cmd
    
    # This is the master options parsing routine, which calls parse_options for
    # the delegate classes (e.g. SpliceConstraints) that handle certain groups
    # of options.
    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "hvp:m:F:a:i:I:G:r:o:j:g:QC", 
                                        ["version",
                                         "help",  
                                         "output-dir=",
                                         "solexa-quals",
                                         "solexa1.3-quals",
                                         "phred64-quals",
                                         "quals",
                                         "integer-quals",
                                         "color",
                                         "color-out",
                                         "library-type=",
                                         "num-threads=",
                                         "splice-mismatches=",
                                         "max-multihits=",
                                         "min-isoform-fraction=",
                                         "min-anchor-length=",
                                         "min-intron-length=",
                                         "max-intron-length=",
                                         "GTF=",
                                         "raw-juncs=",
                                         "no-novel-juncs",
                                         "no-gtf-juncs",
                                         "skip-check-reads",
                                         "mate-inner-dist=",
                                         "mate-std-dev=",
                                         "no-closure-search",
                                         "no-coverage-search",
                                         "closure-search",
                                         "coverage-search",
                                         "microexon-search",
                                         "min-closure-exon=",
                                         "min-closure-intron=",
                                         "max-closure-intron=",
                                         "min-coverage-intron=",
                                         "max-coverage-intron=",
                                         "min-segment-intron=",
                                         "max-segment-intron=",
                                         "segment-length=",
                                         "segment-mismatches=",
                                         "butterfly-search",
                                         "no-butterfly-search",
                                         "keep-tmp",
                                         "rg-id=",
                                         "rg-sample=",
                                         "rg-library=",
                                         "rg-description=",
                                         "rg-platform-unit=",
                                         "rg-center=",
                                         "rg-date=",
                                         "rg-platform=",
                                         "tmp-dir="])
        except getopt.error, msg:
            raise Usage(msg)
        
        self.splice_constraints.parse_options(opts)
        self.system_params.parse_options(opts)
        self.read_params.parse_options(opts)
        self.search_params.parse_options(opts)
        
        global output_dir
        global logging_dir
        global tmp_dir
        global sam_header
        
        custom_tmp_dir = None
        custom_out_dir = None
        # option processing
        for option, value in opts:
            if option in ("-v", "--version"):
                print "TopHat v%s" % (get_version())
                sys.exit(0)
            if option in ("-h", "--help"):
                raise Usage(use_message)
            if option in ("-g", "--max-gene-family"):
                self.max_hits = int(value)
            if option in ("-G", "--GTF"):
                self.gff_annotation = value
            if option in ("-j", "--raw-juncs"):
                self.raw_junctions = value
            if option == "--no-novel-juncs":
                self.find_novel_juncs = False
            if option == "--no-gtf-juncs":
                self.find_GFF_juncs = False
            if option == "--skip-check-reads":
                self.skip_check_reads = True
            if option == "--no-coverage-search":
                self.coverage_search = False
            if option == "--no-closure-search":
                self.closure_search = False
            if option == "--coverage-search":
                self.coverage_search = True
            if option == "--closure-search":
                self.closure_search = True
            if option == "--microexon-search":
                self.microexon_search = True
            if option == "--butterfly-search":
                self.butterfly_search = True
            if option == "--no-butterfly-search":
                self.butterfly_search = False
            if option == "--segment-length":
                self.segment_length = int(value)
            if option == "--segment-mismatches":
                self.segment_mismatches = int(value)
            if option in ("-o", "--output-dir"):
                custom_out_dir = value + "/"
            if option == "--tmp-dir":
                custom_tmp_dir = value + "/"
                
        if custom_out_dir != None:
            output_dir = custom_out_dir + "/"
            logging_dir = output_dir + "logs/"
            tmp_dir = output_dir + "tmp/"
            sam_header = tmp_dir + "stub_header.sam" 
        if custom_tmp_dir != None:
            tmp_dir = custom_tmp_dir
            sam_header = tmp_dir + "stub_header.sam" 
    
        if len(args) < 2:
            raise Usage(use_message)
        return args

# Returns the current time in a nice format
def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

# Ensures that the output, logging, and temp directories are present. If not, 
# they are created
def prepare_output_dir():
    
    print >> sys.stderr, "[%s] Preparing output location %s" % (right_now(), output_dir)
    if os.path.exists(output_dir):
        pass
    else:        
        os.mkdir(output_dir)
        
    if os.path.exists(logging_dir):
        pass
    else:        
        os.mkdir(logging_dir)
        
    if os.path.exists(tmp_dir):
        pass
    else:        
        os.mkdir(tmp_dir)

# Check that the Bowtie index specified by the user is present and all files
# are there.
def check_bowtie_index(idx_prefix):
    print >> sys.stderr, "[%s] Checking for Bowtie index files" % right_now()
    
    idx_fwd_1 = idx_prefix + ".1.ebwt"
    idx_fwd_2 = idx_prefix + ".2.ebwt"
    idx_rev_1 = idx_prefix + ".rev.1.ebwt"
    idx_rev_2 = idx_prefix + ".rev.2.ebwt"
    
    if os.path.exists(idx_fwd_1) and \
       os.path.exists(idx_fwd_2) and \
       os.path.exists(idx_rev_1) and \
       os.path.exists(idx_rev_2):
        return 
    else:
        bowtie_idx_env_var = os.environ.get("BOWTIE_INDEXES")
        if bowtie_idx_env_var == None:
            print >> sys.stderr, "Error: Could not find Bowtie index files " + idx_prefix + ".*"
            sys.exit(1)
        idx_prefix = bowtie_idx_env_var + idx_prefix 
        idx_fwd_1 = idx_prefix + ".1.ebwt"
        idx_fwd_2 = idx_prefix + ".2.ebwt"
        idx_rev_1 = idx_prefix + ".rev.1.ebwt"
        idx_rev_2 = idx_prefix + ".rev.2.ebwt"
        
        if os.path.exists(idx_fwd_1) and \
           os.path.exists(idx_fwd_2) and \
           os.path.exists(idx_rev_1) and \
           os.path.exists(idx_rev_2):
            return 
        else:
            print >> sys.stderr, "Error: Could not find Bowtie index files " + idx_prefix + ".*"
            sys.exit(1)

# Reconstructs the multifasta file from which the Bowtie index was created, if 
# it's not already there.
def bowtie_idx_to_fa(idx_prefix):
    idx_name = idx_prefix.split('/')[-1]
    print >> sys.stderr, "[%s] Reconstituting reference FASTA file from Bowtie index" % (right_now())
    
    try:    
        tmp_fasta_file_name = tmp_dir + idx_name + ".fa"
        tmp_fasta_file = open(tmp_fasta_file_name, "w")

        inspect_log = open(logging_dir + "bowtie_inspect_recons.log", "w")

        inspect_cmd = ["bowtie-inspect",
                       idx_prefix]

        #print >> sys.stderr, "Executing: " + " ".join(inspect_cmd) + " > " + tmp_fasta_file_name   
        ret = subprocess.call(inspect_cmd, 
                              stdout=tmp_fasta_file,
                              stderr=inspect_log)

        # Bowtie reported an error
        if ret != 0:
           print >> sys.stderr, fail_str, "Error: bowtie-inspect returned an error"
           sys.exit(1)
           
    # Bowtie not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie-inspect not found on this system.  Did you forget to include it in your PATH?"
  
    return tmp_fasta_file_name

# Checks whether the multifasta file for the genome is present alongside the 
# Bowtie index files for it.
def check_fasta(idx_prefix):
    print >> sys.stderr, "[%s] Checking for reference FASTA file" % right_now()
    idx_fasta = idx_prefix + ".fa"
    if os.path.exists(idx_fasta):
        return idx_fasta
    else:
        idx_name = idx_prefix.split('/')[-1]
        bowtie_idx_env_var = os.environ.get("BOWTIE_INDEXES")
        if bowtie_idx_env_var != None:
            idx_fasta = bowtie_idx_env_var + idx_prefix + ".fa" 
            if os.path.exists(idx_fasta):
                return idx_fasta
        
        print >> sys.stderr, "\tWarning: Could not find FASTA file " + idx_fasta
        idx_fa = bowtie_idx_to_fa(idx_prefix)
        return idx_fa
        #print >> sys.stderr, "Error: Could not find Maq binary fasta file " + idx_bfa
        #sys.exit(1)
    
# Check that both the Bowtie index and the genome's fasta file are present
def check_index(idx_prefix):
    check_bowtie_index(idx_prefix)
    ref_fasta_file = check_fasta(idx_prefix)
    
    return (ref_fasta_file, None)

# Retrive a tuple containing the system's version of Bowtie.  Parsed from 
# `bowtie --version`
def get_bowtie_version():
    try:
        # Launch Bowtie to capture its version info
        proc = subprocess.Popen(['bowtie', '--version'],stdout=subprocess.PIPE)
        stdout_value = proc.communicate()[0]
        bowtie_version = None
        bowtie_out = repr(stdout_value)

        # Find the version identifier
        version_str = "bowtie version "
        ver_str_idx = bowtie_out.find(version_str)
        if ver_str_idx != -1:
            nl = bowtie_out.find("\\n", ver_str_idx)
            version_val = bowtie_out[ver_str_idx + len(version_str):nl]
            bowtie_version = [int(x) for x in version_val.split('.')]
        if len(bowtie_version) == 3:
            bowtie_version.append(0)
        
        return bowtie_version
    except OSError, o:
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: bowtie not found on this system"
       sys.exit(1)
       
# Retrive a tuple containing the system's version of Bowtie.  Parsed from 
# `bowtie --version`
def get_index_sam_header(read_params, idx_prefix):
    try:
        # Launch Bowtie to capture its version info
        bowtie_sam_header_filename = tmp_name()
        bowtie_sam_header_file = open(bowtie_sam_header_filename,"w")

        
        bowtie_header_cmd = ['bowtie', '--sam']
        if read_params.color == True:
            bowtie_header_cmd.append('-C')
        bowtie_header_cmd += [idx_prefix, '/dev/null']
        proc = subprocess.call(bowtie_header_cmd,stdout=bowtie_sam_header_file, stderr=open('/dev/null'))

        bowtie_sam_header_file.close()
        bowtie_sam_header_file = open(bowtie_sam_header_filename,"r")
        
        sam_header_file = open(sam_header, "w")
        
        preamble = []
        sq_dict_lines = []
        CL_header_line = []
        
        for line in bowtie_sam_header_file.readlines():
            line = line.strip()
            if line.find("@SQ") != -1:
                # Sequence dictionary record
                cols = line.split('\t')
                seq_name = None
                for col in cols: 
                    fields = col.split(':')
                    #print fields
                    if len(fields) > 0 and fields[0] == "SN":
                        seq_name = fields[1]
                if seq_name == None:
                    print >> sys.stderr, "Error: malformed sequence dictionary in sam header"
                    sys.exit(1)
                sq_dict_lines.append([seq_name,line])
            elif line.find("CL"):
                continue
            else:
                preamble.append(line)
        #for line in preamble:
        #    print >> sam_header_file, line

        
        print >> sam_header_file, "@HD\tVN:1.0\tSO:sorted"
    
        if read_params.read_group_id and read_params.sample_id:
            rg_str = "@RG\tID:%s\tSM:%s" % (read_params.read_group_id,
                                            read_params.sample_id)
            if read_params.library_id:
                rg_str += "\tLB:%s" % read_params.library_id
            if read_params.description:
                rg_str += "\tDS:%s" % read_params.description
            if read_params.seq_platform_unit:
                rg_str += "\tPU:%s" % read_params.seq_platform_unit
            if read_params.seq_center:
                rg_str += "\tCN:%s" % read_params.seq_center
            if read_params.mate_inner_dist:
                rg_str += "\tPI:%s" % read_params.mate_inner_dist
            if read_params.seq_run_date:
                rg_str += "\tDT:%s" % read_params.seq_run_date
            if read_params.seq_platform:
                rg_str += "\tPL:%s" % read_params.seq_platform
            
            print >> sam_header_file, rg_str
        
        sq_dict_lines.sort(lambda x,y: cmp(x[0],y[0]))
        for [name, line] in sq_dict_lines:
            print >> sam_header_file, line
        print >> sam_header_file, "@PG\tID:TopHat\tVN:%s\tCL:%s" % (get_version(), run_cmd)
        
        sam_header_file.close()
        
        return sam_header
        
    except OSError, o:
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: bowtie not found on this system"
       sys.exit(1)

# Make sure Bowtie is installed and is recent enough to be useful
def check_bowtie():
    print >> sys.stderr, "[%s] Checking for Bowtie" % right_now()
    bowtie_version = get_bowtie_version()
    if bowtie_version == None:
        print >> sys.stderr, "Error: Bowtie not found on this system"
        sys.exit(1)
    # daehwan - check
    elif bowtie_version[1] < 12 or bowtie_version[2] < 3:
        print >> sys.stderr, "Error: TopHat requires Bowtie 0.12.3 or later"
        sys.exit(1)
    print >> sys.stderr, "\tBowtie version:\t\t\t %s" % ".".join([str(x) for x in bowtie_version])
    

# Retrive a tuple containing the system's version of samtools.  Parsed from 
# `samtools`
def get_samtools_version():
    try:
        # Launch Bowtie to capture its version info
        proc = subprocess.Popen(['samtools'],stderr=subprocess.PIPE)
        stderr_value = proc.communicate()[1]
        samtools_version = None
        samtools_out = repr(stderr_value)
        #print >> sys.stderr, samtools_out
        # Find the version identifier
        version_str = "Version: "
        ver_str_idx = samtools_out.find(version_str)
        if ver_str_idx != -1:
            nl = samtools_out.find("\\n", ver_str_idx)
            ws = samtools_out.find(" ", ver_str_idx + len("Version: "))
            end = min(nl, ws)
            version_val = samtools_out[ver_str_idx + len(version_str):end]
            #print >> sys.stderr, ws, nl, end, samtools_out[ver_str_idx + len(version_str):ws]
            #samtools_version = [int(x) for x in version_val.split('.')]
            samtools_version = [int(x.split('-')[0]) for x in version_val.split('.')]
        if len(samtools_version) == 3:
            samtools_version.append(0)
        
        return samtools_version
    except OSError, o:
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: samtools not found on this system"
       sys.exit(1)

# Make sure the SAM tools are installed and are recent enough to be useful
def check_samtools():
    print >> sys.stderr, "[%s] Checking for Samtools" % right_now()
    samtools_version = get_samtools_version()
    if samtools_version == None:
        print >> sys.stderr, "Error: Samtools not found on this system"
        sys.exit(1)
    elif  samtools_version[1] < 1 or samtools_version[2] < 7:
        print >> sys.stderr, "Error: TopHat requires Samtools 0.1.7 or later"
        sys.exit(1)
    print >> sys.stderr, "\tSamtools version:\t\t %s" % ".".join([str(x) for x in samtools_version])
      


def fq_next(f, fname, color):
   '''
   basic fastq record iterator  
   as a function returning a tuple: (seqID, sequence_string, qv_string, seq_len)
   '''
   seqid,seqstr,qstr,seq_len='','','',0
   fline=f.readline #shortcut to save a bit of time
   line=fline()
   if not line : return (seqid, seqstr,qstr,seq_len)
   while len(line.rstrip())==0: # skip empty lines
      line=fline()
      if not line : return (seqid, seqstr,qstr,seq_len)
   try:
       if line[0] != "@":
          raise ValueError("Records in Fastq files should start with '@' character")
       seqid = line[1:].rstrip()
       seqstr = fline().rstrip()
       #There may now be more sequence lines, or the "+" quality marker line:
       while True:
          line = fline()
          if not line:
             raise ValueError("Premature end of file (missing quality values for "+seqid+")")
          if line[0] == "+":
             #sequence string ended  
             qtitle = line[1:].rstrip()
             if qtitle and qtitle != seqid:
                raise ValueError("Different read ID for sequence and quality (%s vs %s)" \
                                 % (seqid, qtitle))
             break
          seqstr += line.rstrip() #removes trailing newlines
          #loop until + found
       seq_len = len(seqstr)
       #at least one line of quality data should follow  
       qstrlen=0  
       #now read next lines as quality values until seq_len is reached
       while True:
          line=fline()
          if not line : break #end of file  
          qstr += line.rstrip()
          qstrlen=len(qstr)  
          if (not color and qstrlen >= seq_len) or (color and qstrlen + 1 >= seq_len):  
               break # qv string has reached the length of seq string
          #loop until qv has the same length as seq
       if (not color and seq_len != qstrlen) or (color and seq_len != qstrlen + 1):
           raise ValueError("Length mismatch between sequence and quality strings "+ \
                                "for %s (%i vs %i)." \
                                % (seqid, seq_len, qstrlen))
   except ValueError, err:
        print >> sys.stderr, "\nError encountered parsing file "+fname+":\n "+str(err)
        sys.exit(1)
   #return the record  
   return (seqid, seqstr, qstr, seq_len)

def fa_getline(fline):
    global fasta_linebuf, fasta_lastline
    if fasta_linebuf:
        fasta_lastline=fasta_linebuf
        fasta_linebuf=''
        return fasta_lastline
    fasta_lastline=fline()
    return fasta_lastline

def fa_ungetline():
    global fasta_linebuf, fasta_lastline
    fasta_linebuf=fasta_lastline;
    return fasta_lastline

def fa_init(f):
    global fasta_linebuf, fasta_lastline
    fasta_linebuf=''
    fasta_lastline=''
    
def fa_next(f, fname):
   '''
   basic fasta record iterator  
   implemented as a function returning a tuple: (seqID, sequence_string, seq_len)
   '''
   seqid,seqstr,seq_len='','',0
   fline=f.readline # shortcut to readline function of f
   line=fa_getline(fline) # this will use the buffer line if it's there
   if not line : return (seqid, seqstr, seq_len)
   while len(line.rstrip())==0: # skip empty lines
      line=fline()
      if not line : return (seqid, seqstr,qstr,seq_len)
   try:
       if line[0] != ">":
          raise ValueError("Records in Fasta files should start with '>' character")
       seqid = line[1:].split()[0]
       #more sequence lines, or the ">" quality marker line:
       while True:
          line = fa_getline(fline)
          if not line: break
          if line[0] == '>':
             #next sequence starts here  
             fa_ungetline()
             break
          seqstr += line.rstrip()
          #loop until '>' found
       seq_len = len(seqstr)
       if seq_len < 3:
          raise ValueError("Read %s too short (%i)." \
                           % (seqid, seq_len))
   except ValueError, err:
        print >> sys.stderr, "\nError encountered parsing fasta file "+fname+":\n "+str(err)
        sys.exit(1)
   #Return the record and then continue...  
   return (seqid, seqstr, seq_len)

# check_reads() has several jobs.  It examines the user's reads, one file at a 
# time, and determines the file format, read length, and other properties that 
# are used to set the junction search strategy later on.  
# TODO: When we add support for mixed read lengths, this routine 
# will need to set the seed length differently. 
def check_reads(params, reads_files):
    print >> sys.stderr, "[%s] Checking reads" % right_now()
    bowtie_version = get_bowtie_version()
    
    seed_len = params.seed_length
    format = params.reads_format

    observed_formats = set([])
    observed_scales = set([])
    min_seed_len = 99999
    max_seed_len = 0
    max_qual = -1
    files = reads_files.split(',')

    file_pos = 0
    for f_name in files:
        try:
            f = open(f_name)
        except IOError:
            print >> sys.stderr, "Error: could not open file", f_name
            sys.exit(1)

        # skip lines
        while True:
            file_pos = f.tell()
            first_line = f.readline()
            if first_line[0] in "@>":
                break
            
        if first_line[0] == "@":
            format = "fastq"
        elif first_line[0] == ">":
            format = "fasta"
        else:
            print >> sys.stderr, "Error: file %s does not appear to be a valid FASTA or FASTQ file" % f_name

        observed_formats.add(format)
        f.seek(file_pos)

        line_num = 0
        if format == "fastq":
            while True:
              seqid, seqstr, qstr, seq_len = fq_next(f, f_name, params.color)
              if not seqid: break
              if params.color:
                  seq_len -= 1
                  seqstr = seqstr[1:]
              if seq_len < 20:
                  print >> sys.stderr, "Warning: found a read < 20bp in", f_name
              else:
                  min_seed_len = min(seq_len, min_seed_len)
                  max_seed_len = max(seq_len, max_seed_len)
              max_line_qual = max([ord(x) for x in list(qstr)])
              max_qual = max(max_line_qual, max_qual)   
        elif format == "fasta":
            fa_init(f)
            while True:
                seqid, seqstr, seq_len = fa_next(f, f_name)
                if not seqid: break
                if params.color:
                  seq_len -= 1
                  seqstr = seqstr[1:]
                if seq_len < 20:
                     print >> sys.stderr, "Warning: found a read < 20bp in", f_name
                else:
                     min_seed_len = min(seq_len, min_seed_len)
                     max_seed_len = max(seq_len, max_seed_len)
            
    if len(observed_formats) > 1:
        print >> sys.stderr, "Error: TopHat requires all reads be either FASTQ or FASTA.  Mixing formats is not supported."
        sys.exit(1)

    if seed_len != None:
        seed_len = max(seed_len, min_seed_len)
    else:
        seed_len = max_seed_len
        
    print >> sys.stderr, "\tmin read length: %dbp, max read length: %dbp" % (min_seed_len, max_seed_len)
    print >> sys.stderr, "\tformat:\t\t %s" % format
    if format == "fastq":
        quality_scale = "phred33 (default)"
        if params.solexa_quals and not params.phred64_quals:
            quality_scale = "solexa33 (reads generated with GA pipeline version < 1.3)"
        elif params.phred64_quals:
            quality_scale = "phred64 (reads generated with GA pipeline version >= 1.3)"
        print >> sys.stderr, "\tquality scale:\t %s" % quality_scale
    elif format == "fasta":
        if params.color == True:
            params.integer_quals = True
    
    #print seed_len, format, solexa_scale
    return TopHatParams.ReadParams(params.solexa_quals,
                                   params.phred64_quals,
                                   params.quals,
                                   params.integer_quals,
                                   params.color,
                                   params.color_out,
                                   params.library_type,
                                   seed_len, 
                                   format, 
                                   params.mate_inner_dist, 
                                   params.mate_inner_dist_std_dev,
                                   params.read_group_id,
                                   params.sample_id,
                                   params.library_id,
                                   params.description,
                                   params.seq_platform_unit,
                                   params.seq_center,
                                   params.seq_run_date,
                                   params.seq_platform)

# Format a DateTime as a pretty string.  
# FIXME: Currently doesn't support days!
def formatTD(td):
  hours = td.seconds // 3600
  minutes = (td.seconds % 3600) // 60
  seconds = td.seconds % 60
  return '%02d:%02d:%02d' % (hours, minutes, seconds) 

# Calls the prep_reads executable, which prepares an internal read library.
# The read library features reads with monotonically increasing integer IDs.
# prep_reads also filters out very low complexy or garbage reads as well as 
# polyA reads.
def prep_reads(params, reads_list, quals_list, output_name):    
    reads_suffix = ".fq"
    kept_reads_filename = tmp_dir + output_name + reads_suffix
    
    if os.path.exists(kept_reads_filename):
        os.remove(kept_reads_filename)
    kept_reads = open(kept_reads_filename, "a")
    
    filter_log = open(logging_dir + "prep_reads.log", "w")
    
    filter_cmd = [prog_path("prep_reads")]
    filter_cmd.extend(params.cmd())
    if params.read_params.reads_format == "fastq":
        filter_cmd += ["--fastq"]
    elif params.read_params.reads_format == "fasta":
        filter_cmd += ["--fasta"]

    filter_cmd.append(reads_list)

    if params.read_params.quals == True:
        filter_cmd.append(quals_list)
       
    #print "\t executing: `%s'" % " ".join(filter_cmd)    
    # files = reads_list.split(',')
    # for reads_file in files:
    try:       
        print >> run_log, " ".join(filter_cmd)
        ret = subprocess.call(filter_cmd, 
                              stdout=kept_reads,
                              stderr=filter_log)
                              # Bowtie reported an error
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute prep_reads"
            sys.exit(1)
    # prep_reads not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: prep_reads not found on this system.  Did you forget to include it in your PATH?"
        sys.exit(1)

    return kept_reads_filename

# Call bowtie
def bowtie(params,
           bwt_idx_prefix,
           reads_list,
           reads_format,
           mapped_reads,
           unmapped_reads,
           reads_for_ordering = None,
           phred_thresh = 70,
           extra_output = ""):
    start_time = datetime.now()
    bwt_idx_name = bwt_idx_prefix.split('/')[-1]
    print >> sys.stderr, "[%s] Mapping reads against %s with Bowtie%s" % (start_time.strftime("%c"), bwt_idx_name, extra_output)
    
    # Setup Bowtie output redirects
    #bwt_map = output_dir + mapped_reads
    bwt_map = tmp_name()
    tmp_fname = bwt_map.split('/')[-1]
    bwt_log = open(logging_dir + tmp_fname + ".log", "w")
    
    # Launch Bowtie
    try:    
        bowtie_cmd = ["bowtie"]
        
        if reads_format == "fastq":
            bowtie_cmd += ["-q"]
        elif reads_format == "fasta":
            bowtie_cmd += ["-f"]

        if params.read_params.color:
            bowtie_cmd += ["-C", "--col-keepends"]

        if unmapped_reads != None:
            unmapped_reads_fasta_name = unmapped_reads
            bowtie_cmd += ["--un", unmapped_reads_fasta_name,
                           "--max", "/dev/null"]
        else:
            unmapped_reads_fasta_name = None

        # daehwan - check "-v" vs. "-n"
        bowtie_cmd += ["-n", str(params.segment_mismatches),
                         "-p", str(params.system_params.bowtie_threads),
                         "-k", str(params.max_hits),
                         "-m", str(params.max_hits),
                         bwt_idx_prefix, 
                         reads_list]
        
        bowtie_proc = subprocess.Popen(bowtie_cmd, stdout=subprocess.PIPE, stderr=bwt_log)
        
        # fix_map_cmd = [bin_dir + "fix_map_ordering"]
        fix_map_cmd = [prog_path("fix_map_ordering")]
        if reads_format == "fastq":
            fix_map_cmd += ["--fastq"]
        elif reads_format == "fasta":
            fix_map_cmd += ["--fasta"]
        
        if reads_for_ordering == None:
            reads_for_ordering = reads_list
        
        fix_map_cmd.extend([reads_for_ordering, "-"])
             
        print >> run_log, " ".join(bowtie_cmd),"|", " ".join(fix_map_cmd), ">", bwt_map
        
        bwt_map = mapped_reads 
        fix_order_proc = subprocess.Popen(fix_map_cmd, 
                                          stdin=bowtie_proc.stdout,
                                          stdout=open(mapped_reads, "w"))    
        
        # wait for the whole pipe to finish
        fix_order_proc.communicate()   
            
    # Bowtie not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: Bowtie not found on this system.  Did you forget to include it in your PATH?"
        sys.exit(1)
            
    # Success    
    finish_time = datetime.now()
    duration = finish_time - start_time
    #print >> sys.stderr, "\t\t\t[%s elapsed]" %  formatTD(duration)
    return (bwt_map, unmapped_reads_fasta_name)

# Generate a new temporary filename in the user's tmp directory
def tmp_name():
    tmp_root = tmp_dir
    if os.path.exists(tmp_root):
        pass
    else:        
        os.mkdir(tmp_root)
    return tmp_root + os.tmpnam().split('/')[-1] 

# Retrieve a .juncs file from a GFF file by calling the gtf_juncs executable
def get_gtf_juncs(gff_annotation):
    print >> sys.stderr, "[%s] Reading known junctions from GTF file" % (right_now())
    gtf_juncs_log = open(logging_dir + "gtf_juncs.log", "w")
    
    gff_prefix = gff_annotation.split('/')[-1].split('.')[0]
    
    gtf_juncs_out_name  = tmp_dir + gff_prefix + ".juncs"
    gtf_juncs_out = open(gtf_juncs_out_name, "w")
    
    #gtf_juncs_cmd = [bin_dir + "gtf_juncs", gff_annotation]
    gtf_juncs_cmd=[prog_path("gtf_juncs"), gff_annotation]                 
    try:    
        print >> run_log, " ".join(gtf_juncs_cmd)
        retcode = subprocess.call(gtf_juncs_cmd, 
                                  stderr=gtf_juncs_log,
                                  stdout=gtf_juncs_out)
        # cvg_islands returned an error
        if retcode == 1:
            print >> sys.stderr, "\tWarning: TopHat did not find any junctions in GTF file"
            return (False, gtf_juncs_out_name) 
        elif retcode < 0:
            print >> sys.stderr, fail_str, "Error: GTF junction extraction failed with err =", retcode
            sys.exit(1)
    # cvg_islands not found
    except OSError, o:
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: gtf_juncs not found on this system"
       sys.exit(1)
    return (True, gtf_juncs_out_name)

# Call bowtie-build on the FASTA file of sythetic splice junction sequences
def build_juncs_bwt_index(external_splice_prefix, color):
    print >> sys.stderr, "[%s] Indexing splices" % (right_now())
    bowtie_build_log = open(logging_dir + "bowtie_build.log", "w")
    
    #user_splices_out_prefix  = output_dir + "user_splices_idx"
    
    bowtie_build_cmd = ["bowtie-build"]
    if color == True:
        bowtie_build_cmd += ["-C"]
        
    bowtie_build_cmd += [external_splice_prefix + ".fa",
                         external_splice_prefix]            
    try:    
        print >> run_log, " ".join(bowtie_build_cmd)
        retcode = subprocess.call(bowtie_build_cmd, 
                                 stdout=bowtie_build_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Splice sequence indexing failed with err =", retcode
            sys.exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie-build not found on this system"
        sys.exit(1)
    return external_splice_prefix

# Build a splice index from a .juncs file, suitable for use with specified read
# (or read segment) lengths
def build_juncs_index(min_anchor_length, 
                      read_length,
                      juncs_prefix, 
                      external_juncs,  
                      reference_fasta,
                      color):
    print >> sys.stderr, "[%s] Retrieving sequences for splices" % (right_now())
    
    juncs_file_list = ",".join(external_juncs)
    juncs_db_log = open(logging_dir + "juncs_db.log", "w")
    
    external_splices_out_prefix  = tmp_dir + juncs_prefix
    external_splices_out_name = external_splices_out_prefix + ".fa"
    
    external_splices_out = open(external_splices_out_name, "w")
    # juncs_db_cmd = [bin_dir + "juncs_db",
    juncs_db_cmd = [prog_path("juncs_db"), 
                    str(min_anchor_length),
                    str(read_length),
                    juncs_file_list,
                    reference_fasta]            
    try:    
        print >> run_log, " ".join(juncs_db_cmd)
        retcode = subprocess.call(juncs_db_cmd, 
                                 stderr=juncs_db_log,
                                 stdout=external_splices_out)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Splice sequence retrieval failed with err =", retcode
            sys.exit(1)
    # juncs_db not found
    except OSError, o:
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: juncs_db not found on this system"
       sys.exit(1)
       
    external_splices_out_prefix = build_juncs_bwt_index(external_splices_out_prefix, color)
    return external_splices_out_prefix

# Print out the sam header, embedding the user's specified library properties.
# FIXME: also needs SQ dictionary lines
def write_sam_header(read_params, sam_file):
    print >> sam_file, "@HD\tVN:1.0\tSO:sorted"
    
    if read_params.read_group_id and read_params.sample_id:
        rg_str = "@RG\tID:%s\tSM:%s" % (read_params.read_group_id,
                                        read_params.sample_id)
        if read_params.library_id:
            rg_str += "\tLB:%s" % read_params.library_id
        if read_params.description:
            rg_str += "\tDS:%s" % read_params.description
        if read_params.seq_platform_unit:
            rg_str += "\tPU:%s" % read_params.seq_platform_unit
        if read_params.seq_center:
            rg_str += "\tCN:%s" % read_params.seq_center
        if read_params.mate_inner_dist:
            rg_str += "\tPI:%s" % read_params.mate_inner_dist
        if read_params.seq_run_date:
            rg_str += "\tDT:%s" % read_params.seq_run_date
        if read_params.seq_platform:
            rg_str += "\tPL:%s" % read_params.seq_platform
        
        print >> sam_file, rg_str
    print >> sam_file, "@PG\tID:TopHat\tVN:%s\tCL:%s" % (get_version(), run_cmd)
            
# Write final TopHat output, via tophat_reports and wiggles
def compile_reports(params, sam_header_filename, left_maps, left_reads, right_maps, right_reads, gff_annotation):
    print >> sys.stderr, "[%s] Reporting output tracks" % right_now()
    
    left_maps = [x for x in left_maps if (os.path.exists(x) and os.path.getsize(x) > 0)]
    left_maps = ','.join(left_maps)
    
    if len(right_maps) > 0:
        right_maps = [x for x in right_maps if (os.path.exists(x) and os.path.getsize(x) > 0)]
        right_maps = ','.join(right_maps)
    
    report_log = open(logging_dir + "reports.log", "w")
    junctions = output_dir + "junctions.bed"
    coverage =  "coverage.wig"
    accepted_hits_sam = tmp_dir + "accepted_hits.sam"
    accepted_hits_bam = output_dir + "accepted_hits.bam"
    report_cmdpath = which("tophat_reports")
    if report_cmdpath == None:
        print >> sys.stderr, err_find_prog+"tophat_reports"
        sys.exit(1)
    report_cmd = [report_cmdpath]
    report_cmd.extend(params.cmd())
        
    report_cmd.extend([junctions,
                       accepted_hits_sam,
                       left_maps,
                       left_reads])
    if len(right_maps) > 0 and right_reads != None:
        report_cmd.append(right_maps)
        report_cmd.append(right_reads)
                    
    try: 
    
        accepted_hits_sam_file = open(accepted_hits_sam, "w")
        #write_sam_header(params.read_params, sorted_map)
    
        header = open(sam_header_filename, "r")
        for line in header:
            print >> accepted_hits_sam_file, line,
    
        accepted_hits_sam_file.close()
        #accepted_hits_sam_file = open(accepted_hits_sam, "a")
    
        print >> run_log, " ".join(report_cmd)   
        retcode = subprocess.call(report_cmd, 
                                  stderr=report_log)
       
        # spanning_reads returned an error 
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Report generation failed with err =", retcode
            sys.exit(1)
        
        tmp_bam = tmp_name() 
        sam_to_bam_cmd = ["samtools", "view", "-S", "-b", accepted_hits_sam]
        print >> run_log, " ".join(sam_to_bam_cmd) + " > " + tmp_bam
        sam_to_bam_log = open(logging_dir + "accepted_hits_sam_to_bam.log", "w")
        tmp_bam_file = open(tmp_bam, "w")
        ret = subprocess.call(sam_to_bam_cmd, 
                              stdout=tmp_bam_file,
                              stderr=sam_to_bam_log)
        if ret != 0:
            print >> sys.stderr, "Error: could not convert to BAM with samtools"
            sys.exit(1)
        sort_cmd = ["samtools", "sort", tmp_bam, output_dir + "accepted_hits"]
        print >> run_log, " ".join(sort_cmd)
#        sorted_map_name = tmp_name()
#        sorted_map = open(sorted_map_name, "w")
#        write_sam_header(params.read_params, sorted_map)
#        sorted_map.close()
#        sorted_map = open(sorted_map_name, "a")
#        
#        sort_cmd =["sort",
#                    "-k",
#                    "3,3", 
#                    "-k", 
#                    "4,4n",
#                    "--temporary-directory="+tmp_dir,
#                    output_dir + accepted_hits]
#                    
                
        sort_bam_log = open(logging_dir + "accepted_hits_bam_sort.log", "w")
        ret = subprocess.call(sort_cmd, 
                        stdout=open('/dev/null'),
                        stderr=sort_bam_log)
        if ret != 0:
            print >> sys.stderr, "Error: could not sort BAM file with samtools"
            sys.exit(1)        
            
        os.remove(accepted_hits_sam)
#        print >> run_log, "mv %s %s" % (sorted_map, output_dir + accepted_hits)
#        os.rename(sorted_map_name, output_dir + accepted_hits) 

# FIXME: put wiggles back!
#        wig_cmd = [prog_path("wiggles"), output_dir + accepted_hits, output_dir + coverage]
#        print >> run_log, " ".join(wig_cmd)
#        subprocess.call(wig_cmd,
#                        stderr=open("/dev/null"))
                        
    # cvg_islands not found
    except OSError, o:
        print >>sys.stderr, o
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: tophat_reports not found on this system"
        sys.exit(1)
    return (coverage, junctions)


# Split up each read in a FASTQ file into multiple segments. Creates a FASTQ file
# for each segment  This function needs to be fixed to support mixed read length
# inputs
def split_reads(reads_filename, 
                prefix,
                fasta,
                color,
                segment_length):
    reads_file = open(reads_filename)
    output_files = []
    
    if fasta == True:
        extension = ".fa"
    else:
        extension = ".fq"

    def open_output_files(prefix, num_files_prev, num_files, output_files, extension):
        i = num_files_prev + 1
        while i <= num_files:
            output_files.append(open(prefix + ("_seg%d" % i) + extension, "w"))
            i += 1

    def convert_color_to_bp(color_seq):
        decode_dic = { 'A0':'A', 'A1':'C', 'A2':'G', 'A3':'T', 'A4':'N', 'A.':'N',
                           'C0':'C', 'C1':'A', 'C2':'T', 'C3':'G', 'C4':'N', 'C.':'N',
                           'G0':'G', 'G1':'T', 'G2':'A', 'G3':'C', 'G4':'N', 'G.':'N',
                           'T0':'T', 'T1':'G', 'T2':'C', 'T3':'A', 'T4':'N', 'T.':'N',
                           'N0':'N', 'N1':'N', 'N2':'N', 'N3':'N', 'N4':'N', 'N.':'N' }

        base = color_seq[0]
        bp_seq = base
        for ch in color_seq[1:]:
            base = decode_dic[base+ch]
            bp_seq += base
        return bp_seq

    def convert_bp_to_color(bp_seq):
        encode_dic = { 'AA':'0', 'CC':'0', 'GG':'0', 'TT':'0',
                       'AC':'1', 'CA':'1', 'GT':'1', 'TG':'1',
                       'AG':'2', 'CT':'2', 'GA':'2', 'TC':'2',
                       'AT':'3', 'CG':'3', 'GC':'3', 'TA':'3',
                       'A.':'4', 'C.':'4', 'G.':'4', 'T.':'4',
                       '.A':'4', '.C':'4', '.G':'4', '.T':'4',
                       '.N':'4', 'AN':'4', 'CN':'4', 'GN':'4',
                       'TN':'4', 'NA':'4', 'NC':'4', 'NG':'4',
                       'NT':'4', 'NN':'4', 'N.':'4', '..':'4' }

        base = bp_seq[0]
        color_seq = base
        for ch in bp_seq[1:]:
            color_seq += encode_dic[base + ch]
            base = ch

        return color_seq

    def split_record(read_name, read_seq, read_qual, output_files, offsets, color):
        if color == True:
            color_offset = 1
            read_seq_temp = convert_color_to_bp(read_seq)

            seg_num = 1
            while seg_num + 1 < len(offsets):
                if read_seq[offsets[seg_num]+1] not in ['0', '1', '2', '3']:
                    return
                seg_num += 1
        else:
            color_offset = 0

        seg_num = 0
        last_seq_offset = 0
        while seg_num + 1 < len(offsets):
            f = output_files[seg_num]
            seg_seq = read_seq[last_seq_offset+color_offset:offsets[seg_num + 1]+color_offset]
            print >> f, "%s|%d:%d:%d" % (read_name,last_seq_offset,seg_num, len(offsets) - 1)
            if color == True:
                print >> f, "%s%s" % (read_seq_temp[last_seq_offset], seg_seq)
            else:
                print >> f, seg_seq
            if fasta == False:
                seg_qual = read_qual[last_seq_offset:offsets[seg_num + 1]]
                print >> f, "+"
                print >> f, seg_qual
            seg_num += 1
            last_seq_offset = offsets[seg_num]

    line_state = 0
    read_name = ""
    read_seq = ""
    read_qual = ""
    num_segments = 0
    offsets = []
    for line in reads_file:
        if line.strip() == "":
            continue
        if line_state == 0:
            read_name = line.strip()
        elif line_state == 1:
            read_seq = line.strip()

            read_length = len(read_seq)
            tmp_num_segments = read_length / segment_length
            offsets = [segment_length * i for i in range(0, tmp_num_segments + 1)]

            # Bowtie's minimum read length here is 20bp, so if the last segment
            # is between 20 and segment_length bp long, go ahead and write it out
            if read_length % segment_length >= 20:
                offsets.append(read_length)
                tmp_num_segments += 1
            else:
                offsets[-1] = read_length

            if tmp_num_segments == 1:
                offsets = [0, read_length]

            if tmp_num_segments > num_segments:
                open_output_files(prefix, num_segments, tmp_num_segments, output_files, extension)
                num_segments = tmp_num_segments

            if fasta:
                split_record(read_name, read_seq, None, output_files, offsets, color)
        elif line_state == 2:
            line = line.strip()
        else:
            read_quals = line.strip()
            if not fasta:
                split_record(read_name, read_seq, read_quals, output_files, offsets, color)
                
        line_state += 1
        if fasta:
            line_state %= 2
        else:
            line_state %= 4
        
    for f in output_files:
        f.close()
    return [o.name for o in output_files]

# Find possible splice junctions using the "closure search" strategy, and report
# them in closures.juncs.  Calls the executable closure_juncs
def junctions_from_closures(params,
                            left_maps, 
                            right_maps, 
                            ref_fasta):
    #print >> sys.stderr, "[%s] " % right_now()
    print >> sys.stderr, "[%s] Searching for junctions via mate-pair closures" % right_now()
    #maps = [x for x in seg_maps if (os.path.exists(x) and os.path.getsize(x) > 0)]
    #if len(maps) == 0:
    #    return None
    #print >> sys.stderr, left_maps
    slash = left_maps[0].rfind('/')
    juncs_out = ""
    if slash != -1:
        juncs_out += left_maps[0][:slash+1]
    juncs_out += "closure.juncs"

    juncs_log = open(logging_dir + "closure.log", "w")
    juncs_cmdpath=prog_path("closure_juncs")
    juncs_cmd = [juncs_cmdpath]

    left_maps = ','.join(left_maps);
    right_maps = ','.join(right_maps);

    juncs_cmd.extend(params.cmd())
    juncs_cmd.extend([juncs_out,
                      ref_fasta,
                      left_maps,
                      right_maps])            
    try:
        print >> run_log, ' '.join(juncs_cmd)
        retcode = subprocess.call(juncs_cmd, 
                                 stderr=juncs_log)

        # spanning_reads returned an error 
        if retcode != 0:
           print >> sys.stderr, fail_str, "Error: closure-based junction search failed with err =", retcode
           sys.exit(1)
    # cvg_islands not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: closure_juncs not found on this system"
        sys.exit(1)
    return [juncs_out]

# Find possible junctions by examining coverage and split segments in the initial
# map and segment maps.  Report junctions in segment.juncs.  Calls the executable
# segment_juncs
def junctions_from_segments(params,
                            left_reads,
                            left_reads_map,
                            left_seg_maps, 
                            right_reads,
                            right_reads_map,
                            right_seg_maps,
                            unmapped_reads, 
                            reads_format, 
                            ref_fasta):
    print >> sys.stderr, "[%s] Searching for junctions via segment mapping" % right_now()
    slash = left_seg_maps[0].rfind('/')
    juncs_out = ""
    if slash != -1:
        juncs_out += left_seg_maps[0][:slash+1]
    juncs_out += "segment.juncs"

    left_maps = ','.join(left_seg_maps)
    align_log = open(logging_dir + "segment_juncs.log", "w")
    align_cmd = [prog_path("segment_juncs")]
    
    align_cmd.extend(params.cmd())
    
    align_cmd.extend(["--ium-reads", ",".join(unmapped_reads),
                      ref_fasta,
                      juncs_out,
                      left_reads,
                      left_reads_map,
                      left_maps])
    if right_seg_maps != None:
        right_maps = ','.join(right_seg_maps)
        align_cmd.extend([right_reads, right_reads_map, right_maps])            
    try:
        print >> run_log, " ".join(align_cmd)
        retcode = subprocess.call(align_cmd, 
                                 stderr=align_log)

        # spanning_reads returned an error 
        if retcode != 0:
           print >> sys.stderr, fail_str, "Error: segment-based junction search failed with err =",retcode
           sys.exit(1)
    # cvg_islands not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: segment_juncs not found on this system"
        sys.exit(1)

    return [juncs_out]

# Joins mapped segments into full-length read alignments via the executable
# long_spanning_reads
def join_mapped_segments(params,
                         sam_header_filename,
                         reads,
                         ref_fasta,
                         possible_juncs,
                         contig_seg_maps,
                         spliced_seg_maps,
                         alignments_out_name):
    print >> sys.stderr, "[%s] Joining segment hits" % right_now()
    contig_seg_maps = ','.join(contig_seg_maps)
    possible_juncs = ','.join(possible_juncs)
    
    align_log = open(logging_dir + "long_spanning_reads.log", "w")
    align_cmd = [prog_path("long_spanning_reads")]
    
    # if params.read_params.reads_format == "fastq":
    #     align_cmd += ["-q"]
    # elif params.read_params.reads_format == "fasta":
    #     align_cmd += ["-f"]
    
    alignments_out = open(alignments_out_name, "w")
    #write_sam_header(params.read_params, sorted_map)
    
    header = open(sam_header_filename, "r")
    for line in header:
        print >> alignments_out, line,
    
    alignments_out.close()
    alignments_out = open(alignments_out_name, "a")
    
    #alignments_out = open(alignments_out_name, "w")
    
    align_cmd.extend(params.cmd())

    align_cmd.append(ref_fasta)
    align_cmd.extend([ reads,
                       possible_juncs,
                       contig_seg_maps])
    if spliced_seg_maps != None:
        spliced_seg_maps = ','.join(spliced_seg_maps)
        align_cmd.append(spliced_seg_maps)
                
    try:    
        print >> run_log, " ".join(align_cmd),">", alignments_out_name
        retcode = subprocess.call(align_cmd, 
                                  stderr=align_log,
                                  stdout=alignments_out)

        # spanning_reads returned an error 
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Segment join failed with err =", retcode
            sys.exit(1)
     # cvg_islands not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: long_spanning_reads not found on this system"
        sys.exit(1)
      

# The main aligment routine of TopHat.  This function executes most of the 
# workflow producing a set of candidate alignments for each cDNA fragment in a 
# pair of SAM alignment files (for paired end reads).
def spliced_alignment(params,
                      bwt_idx_prefix,
                      sam_header_filename,
                      ref_fasta,
                      read_len,
                      segment_len,
                      left_reads,
                      right_reads,
                      user_supplied_junctions):
    
    possible_juncs = []
    possible_juncs.extend(user_supplied_junctions)
    
    left_maps = []
    right_maps = []
    maps = { left_reads : [], right_reads : [] }
    #single_segments = False

    # This class collects spliced and unspliced alignments for each of the 
    # left and right read files provided by the user.
    class Maps:
        def __init__(self, 
                     unspliced_bwt,
                     unspliced_sam,
                     seg_maps, 
                     unmapped_segs, 
                     segs):
            self.unspliced_bwt = unspliced_bwt
            self.unspliced_sam = unspliced_sam
            self.seg_maps = seg_maps
            self.unmapped_segs = unmapped_segs
            self.segs = segs
    
    # Perform the first part of the TopHat work flow on the left and right 
    # reads of paired ends separately - we'll use the pairing information later
    for reads in [left_reads, right_reads]:
        if reads == None:
            continue
        slash = reads.rfind("/")
        extension = reads.rfind(".")
        if extension != -1:
            prefix = reads[slash+1:extension]
        else:
            prefix = reads      
          
        extension = reads.rfind(".")
        assert extension != -1
        tmp = reads.rfind("/")

        unspliced_out = tmp_dir + reads[tmp+1:extension] + ".bwtout"  
        unmapped_unspliced = tmp_dir + reads[tmp+1:extension] + "_missing.fq"

        num_segs = read_len / segment_len
        if read_len % segment_len >= 20:
            num_segs += 1

        # Perform the initial Bowtie mapping of the full lenth reads
        (unspliced_map, unmapped) = bowtie(params,
                                           bwt_idx_prefix, 
                                           reads,
                                           "fastq",
                                           unspliced_out,
                                           unmapped_unspliced,
                                           reads)

        unspliced_sam = tmp_name()
        
        # Convert the initial Bowtie maps into SAM format.  
        # TODO: this step should be removable now that Bowtie supports SAM 
        # format
        join_mapped_segments(params,
                             sam_header_filename,
                             reads,
                             ref_fasta,
                             ["/dev/null"],
                             [unspliced_map],
                             [],
                             unspliced_sam)

        # Using the num_segs value returned by check_reads(), decide which 
        # junction discovery strategy to use
        if num_segs == 1:
            segment_len = read_len
        elif num_segs >= 3:
            # if we have at least three segments, just use split segment search,
            # which is the most sensitive and specific, fastest, and lightest-weight.
            if params.closure_search != True:
                params.closure_search = False
            if params.coverage_search != True:
                params.coverage_search = False
            if params.butterfly_search != True:
                params.butterfly_search = False
        if params.closure_search != True:
            params.closure_search = False
                                                       
        seg_maps = []
        unmapped_segs = []
        segs = []
        if num_segs > 1:
            # split up the IUM reads into segments
            read_segments = split_reads(unmapped_unspliced,
                                        tmp_dir + prefix, 
                                        False,
                                        params.read_params.color,
                                        segment_len)

            # Map each segment file independently with Bowtie
            for i in range(len(read_segments)):
                seg = read_segments[i]
                extension = seg.rfind(".")
                assert extension != -1
                tmp = seg.rfind("/")
                assert tmp != -1
                seg_out =  tmp_dir + seg[tmp+1:extension] + ".bwtout"
                unmapped_seg = tmp_dir + seg[tmp+1:extension] + "_missing.fq"
                extra_output = "(%d/%d)" % (i+1, len(read_segments))
                (seg_map, unmapped) = bowtie(params,
                                             bwt_idx_prefix, 
                                             seg,
                                             "fastq",
                                             seg_out,
                                             unmapped_seg,
                                             seg,
                                             70,
                                             extra_output)
                seg_maps.append(seg_map)
                unmapped_segs.append(unmapped)
                segs.append(seg)
            
            # Collect the segment maps for left and right reads together
            maps[reads] = Maps(unspliced_map, unspliced_sam, seg_maps, unmapped_segs, segs) 
        else:
            # if there's only one segment, just collect the initial map as the only
            # map to be used downstream for coverage-based junction discovery
            read_segments = [reads]
            maps[reads] = Maps(unspliced_map, unspliced_sam, [unspliced_map], [unmapped_unspliced], read_segments)
    
    # Unless the user asked not to discover new junctions, start that process
    # here
    if params.find_novel_juncs:
        left_reads_map = maps[left_reads].unspliced_bwt
        left_seg_maps = maps[left_reads].seg_maps
        unmapped_reads = maps[left_reads].unmapped_segs
        if right_reads != None:
            right_reads_map = maps[right_reads].unspliced_bwt
            right_seg_maps = maps[right_reads].seg_maps
            unmapped_reads.extend(maps[right_reads].unmapped_segs)
        else:
            right_reads_map = None
            right_seg_maps = None
        
        # Call segment_juncs to infer a list of possible splice junctions from
        # the regions of the genome covered in the initial and segment maps
        juncs = junctions_from_segments(params, 
                                        left_reads,
                                        left_reads_map,
                                        left_seg_maps, 
                                        right_reads,
                                        right_reads_map,
                                        right_seg_maps,
                                        unmapped_reads,
                                        "fastq", 
                                        ref_fasta)
        if os.path.getsize(juncs[0]) != 0:
            possible_juncs.extend(juncs)
    
        # Optionally, and for paired reads only, use a closure search to 
        # discover addtional junctions
        if params.closure_search and left_reads != None and right_reads != None:
            juncs = junctions_from_closures(params,
                                            [maps[left_reads].unspliced_bwt, maps[left_reads].seg_maps[-1]],
                                            [maps[right_reads].unspliced_bwt, maps[right_reads].seg_maps[-1]],
                                            ref_fasta)
            if os.path.getsize(juncs[0]) != 0:
                possible_juncs.extend(juncs)

    if len(possible_juncs) == 0:
        spliced_seg_maps = None
        junc_idx_prefix = None
        print >> sys.stderr, "Warning: junction database is empty!"
    else:  
        # index the junction sequences with bowtie-build
        junc_idx_prefix = "segment_juncs"
        build_juncs_index(3,
                          segment_len,
                          junc_idx_prefix, 
                          possible_juncs,
                          ref_fasta,
                          params.read_params.color)
    
    # Now map read segments (or whole IUM reads, if num_segs == 1) to the splice
    # index with Bowtie
    for reads in [left_reads, right_reads]:
        spliced_seg_maps = []
        if reads == None or reads not in maps:
            continue
        
        if junc_idx_prefix != None:
            i = 0    
            for seg in maps[reads].segs:
                extension = seg.rfind(".")
                assert extension != -1
                tmp = seg.rfind("/")
                assert tmp != -1
            
                ordering = maps[reads].segs[i]
                seg_out = tmp_dir + seg[tmp+1:extension] + "_to_spliced.bwtout"
                (seg_map, unmapped) = bowtie(params,
                                             tmp_dir + junc_idx_prefix, 
                                             seg,
                                             "fastq",
                                             seg_out,
                                             None,
                                             ordering)
                spliced_seg_maps.append(seg_map)
                i += 1
        
        # Join the contigous and spliced segment hits into full length read 
        # alignments
        mapped_reads = reads +".candidate_hits.sam"
        join_mapped_segments(params,
                             sam_header_filename,
                             reads,
                             ref_fasta,
                             possible_juncs,
                             maps[reads].seg_maps,
                             spliced_seg_maps,
                             mapped_reads)

        if num_segs > 1:
            # Merge the spliced and unspliced full length alignments into a 
            # single SAM file.
            # FIXME: we ought to be able to do this much more efficiently,
            # since the individual SAM files are all already sorted in 
            # increasing read ID order.  We should just be able to merge these
            # NOTE: We also should be able to address bug #134 here, by replacing
            # contiguous alignments that poke into an intron by a small amount by
            # the correct spliced alignment.
            
            # Using this requires converting to BAM upstream.  What a pain.
#            merged_map = tmp_name() + ".sam"
#
#            merge_sort_cmd =["samtools", "merge", "-n", "-h",
#                              sam_header_filename,
#                              maps[reads].unspliced_sam, 
#                              mapped_reads,
#                              merged_map]
                              
            merged_map = tmp_name()
            merge_sort_cmd =["sort",
                             "-k 1,1n",
                              "--temporary-directory="+tmp_dir,
                              maps[reads].unspliced_sam, 
                              mapped_reads]
            print >> run_log, " ".join(merge_sort_cmd), ">", merged_map

            print >> run_log, " ".join(merge_sort_cmd)
            subprocess.call(merge_sort_cmd,
                             stdout=open(merged_map,"w")) 
        else:
            merged_map = mapped_reads
        
        if params.system_params.keep_tmp == False:
            os.remove(maps[reads].unspliced_bwt)
            os.remove(maps[reads].unspliced_sam) 
        print >> run_log, "mv %s %s" % (merged_map, mapped_reads)           
        os.rename(merged_map, mapped_reads)
                          
        maps[reads] = [mapped_reads]
    return maps

# rough equivalent of the 'which' command to find external programs
# (current script path is tested first, then PATH envvar)
def which(program):
    def is_executable(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_executable(program):
            return program
    else:
        progpath = os.path.join(bin_dir, program); 
        if is_executable(progpath):
           return progpath
        for path in os.environ["PATH"].split(os.pathsep):
           progpath = os.path.join(path, program)
           if is_executable(progpath):
              return progpath
    return None

def prog_path(program):
    progpath=which(program)
    if progpath == None:
        print >> sys.stderr, "Error locating program: ", program
        sys.exit(1)
    return progpath

# FIXME: this should get set during the make dist autotools phase of the build
def get_version():
   return "1.1.4"

def main(argv=None):
    warnings.filterwarnings("ignore", "tmpnam is a potential security risk")
    
    # Initialize default parameter values
    params = TopHatParams()
    
    try:
        if argv is None:
            argv = sys.argv
            args = params.parse_options(argv)
            params.check()
        
        bwt_idx_prefix = args[0]
        left_reads_list = args[1]
        left_quals_list, right_quals_list = [], []
        if (params.read_params.quals != True and len(args) > 2) or (params.read_params.quals == True and len(args) > 3):
            if params.read_params.mate_inner_dist == None:
                print >> sys.stderr, "Error: you must set the mean inner distance between mates with -r"
                sys.exit(1)

            right_reads_list = args[2]
            if params.read_params.quals == True:
                left_quals_list = args[3]
                right_quals_list = args[4]
        else:
            right_reads_list = None
            if params.read_params.quals == True:
                left_quals_list = args[2]
            
        print >> sys.stderr
        print >> sys.stderr, "[%s] Beginning TopHat run (v%s)" % (right_now(), get_version())
        print >> sys.stderr, "-----------------------------------------------" 
        
        start_time = datetime.now()
        prepare_output_dir()
        
        global run_log
        run_log = open(logging_dir + "run.log", "w", 0)
        global run_cmd
        run_cmd = " ".join(argv)
        print >> run_log, run_cmd
        
        
        # Validate all the input files, check all prereqs before committing 
        # to the run
        (ref_fasta, ref_seq_dict) = check_index(bwt_idx_prefix)
        
        check_bowtie()
        check_samtools()
        
        sam_header_filename = get_index_sam_header(params.read_params, bwt_idx_prefix)
        
        if params.skip_check_reads == False:
            reads_list = left_reads_list
            if right_reads_list:
                reads_list = reads_list + "," + right_reads_list
            params.read_params = check_reads(params.read_params, reads_list)
            
        user_supplied_juncs = []
        if params.gff_annotation != None and params.find_GFF_juncs == True:
            (found_juncs, gtf_juncs) = get_gtf_juncs(params.gff_annotation)
            if found_juncs == True:
                user_supplied_juncs.append(gtf_juncs)
        if params.raw_junctions != None:
            user_supplied_juncs.append(params.raw_junctions)
                
        # Now start the time consuming stuff
        left_kept_reads = prep_reads(params,
                                     left_reads_list,
                                     left_quals_list,
                                     "left_kept_reads")

        if right_reads_list != None:
            right_kept_reads = prep_reads(params,
                                          right_reads_list,
                                          right_quals_list,
                                          "right_kept_reads")
        else:
            right_kept_reads = None

        # turn off integer-quals
        if params.read_params.integer_quals == True:
            params.read_params.integer_quals = False
            
        spliced_reads = []
        mapping = spliced_alignment(params, 
                                    bwt_idx_prefix,
                                    sam_header_filename,
                                    ref_fasta,
                                    params.read_params.seed_length,
                                    params.segment_length,
                                    left_kept_reads,
                                    right_kept_reads,
                                    user_supplied_juncs)
                                    
        left_maps = mapping[left_kept_reads]
        #left_unmapped_reads = mapping[1]
        right_maps = mapping[right_kept_reads]
        #right_unmapped_reads = mapping[3]

        compile_reports(params,
                        sam_header_filename,
                        left_maps,
                        left_kept_reads,
                        right_maps,
                        right_kept_reads,
                        params.gff_annotation)

        if params.system_params.keep_tmp == False:
            for m in left_maps:
                os.remove(m)
            if left_kept_reads != None:
                os.remove(left_kept_reads)
            for m in right_maps:
                os.remove(m)
            if right_kept_reads != None:
                os.remove(right_kept_reads)
            tmp_files = os.listdir(tmp_dir)
            for t in tmp_files:
                os.remove(tmp_dir+t)
            os.rmdir(tmp_dir)


        finish_time = datetime.now()
        duration = finish_time - start_time
        print >> sys.stderr,"-----------------------------------------------"
        print >> sys.stderr, "Run complete [%s elapsed]" %  formatTD(duration)
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "    for detailed help see http://tophat.cbcb.umd.edu/manual.html"
        return 2


if __name__ == "__main__":
    sys.exit(main())
