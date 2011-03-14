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

import getopt
import subprocess
import errno
import os
import warnings
import re
from datetime import datetime, date, time

use_message = '''
TopHat maps short sequences from spliced transcripts to whole genomes.

Usage:
    tophat [options] <bowtie_index> <reads1[,reads2,...]> [reads1[,reads2,...]] \\
                                    [quals1,[quals2,...]] [quals1[,quals2,...]]
    
Options:
    -v/--version
    -o/--output-dir                <string>    [ default: ./tophat_out     ]
    -a/--min-anchor                <int>       [ default: 8                ]
    -m/--splice-mismatches         <0-2>       [ default: 0                ]
    -i/--min-intron-length         <int>       [ default: 50               ]
    -I/--max-intron-length         <int>       [ default: 500000           ]
    -g/--max-multihits             <int>       [ default: 20               ]
    -F/--min-isoform-fraction      <float>     [ default: 0.15             ]
    --max-insertion-length         <int>       [ default: 3                ]
    --max-deletion-length          <int>       [ default: 3                ]
    --solexa-quals                          
    --solexa1.3-quals                          (same as phred64-quals)
    --phred64-quals                            (same as solexa1.3-quals)
    -Q/--quals
    --integer-quals
    -C/--color                                 (Solid - color space)
    --color-out
    --library-type                 <string>    (fr-unstranded, fr-firststrand, 
                                                fr-secondstrand)
    -p/--num-threads               <int>       [ default: 1                ]
    -G/--GTF                       <filename>
    -j/--raw-juncs                 <filename>
    --insertions                   <filename>
    --deletions                    <filename>
    -r/--mate-inner-dist           <int>       
    --mate-std-dev                 <int>       [ default: 20               ]
    --no-novel-juncs
    --allow-indels
    --no-novel-indels              
    --no-gtf-juncs                 
    --no-coverage-search
    --coverage-search              
    --no-closure-search
    --closure-search      
    --microexon-search
    --butterfly-search
    --no-butterfly-search
    --keep-tmp
    --tmp-dir                      <dirname>   [ default: <output_dir>/tmp ]
    -z/--zpacker                   <filepath>  [ default: gzip             ]
    
Advanced Options:
    --segment-mismatches           <int>       [ default: 2                ]
    --segment-length               <int>       [ default: 25               ]
    --bowtie-n                                 [ default: bowtie -v        ]
    --min-closure-exon             <int>       [ default: 100              ]
    --min-closure-intron           <int>       [ default: 50               ]
    --max-closure-intron           <int>       [ default: 5000             ]
    --min-coverage-intron          <int>       [ default: 50               ]
    --max-coverage-intron          <int>       [ default: 20000            ]
    --min-segment-intron           <int>       [ default: 50               ]
    --max-segment-intron           <int>       [ default: 500000           ]

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
tmp_dir = output_dir + "tmp/"
bin_dir = sys.path[0] + "/"
unmapped_reads_fifo = None # tricking bowtie into writing the unmapped reads into a compressed file
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
                die("Error: arg to --splice-mismatches must be 0, 1, or 2")
            if self.min_anchor_length < 4:
                die("Error: arg to --min-anchor-len must be greater than 4")
            if self.min_isoform_fraction < 0.0 or self.min_isoform_fraction > 1.0:
                die("Error: arg to --min-isoform-fraction must be between 0.0 and 1.0")
            if self.min_intron_length <= 0:
                die("Error: arg to --min-intron-length must be greater than 0")
            if self.max_intron_length <= 0:
                die("Error: arg to --max-intron-length must be greater than 0")
    
    # SystemParams is a group of runtime parameters that determine how to handle
    # temporary files produced during a run and how many threads to use for threaded
    # stages of the pipeline (e.g. Bowtie)
    
    class SystemParams:
        def __init__(self,
                     num_cpus,
                     keep_tmp):
            self.num_cpus = num_cpus
            self.keep_tmp = keep_tmp
            self.zipper = 'gzip'
            self.zipper_opts= []
            
        def parse_options(self, opts):
            for option, value in opts:
                if option in ("-p", "--num-threads"):
                    self.num_cpus = int(value)
                elif option == "--keep-tmp":
                    self.keep_tmp = True
                elif option in ("-z","--zpacker"):
                    self.zipper = value
                    if not self.zipper:
                       self.zipper='gzip'
        def cmd(self):
            cmdline=[]
            if self.zipper and self.zipper!='gzip':
                 cmdline.extend(['-z',self.zipper])
            if self.num_cpus>1:
                 cmdline.extend(['-p'+str(self.num_cpus)])
            return cmdline
        
        def check(self):
            if self.num_cpus<1 :
                 die("Error: arg to --num-threads must be greater than 0")
            if self.num_cpus>1 :
                if self.zipper.endswith("gzip"): # try pigz instead
                    pigz=which("pigz")
                    if (pigz is not None): 
                        self.zipper=pigz
                    #    self.zipper_opts.append('-p'+str(self.num_cpus))
                    #else:
                    #    print >> sys.stderr, "Consider installing 'pigz' for faster handling of compressed temporary files."
                elif self.zipper.endswith("bzip2") and not self.zipper.endswith("pbzip2"): # try pbzip2 instead
                    pbzip=which("pbzip2")
                    if (pbzip is not None): 
                        self.zipper=pbzip
                    #    self.zipper_opts.append('-p'+str(self.num_cpus))
                    #else:
                    #    print >> sys.stderr, "Consider installing 'pigz' or 'pbzip2' for faster handling of temporary files."
        
            if self.zipper:  #this should ALWAYS be the case
                xzip=which(self.zipper)
                if not xzip:
                    die("Error: cannot find compression program "+xzip)
                if self.num_cpus>1 and not self.zipper_opts:
                    if self.zipper.endswith('pbzip2') or self.zipper.endswith('pigz'):
                         self.zipper_opts.append('-p'+str(self.num_cpus))
            else: 
                 die("Error: no compression program!") #this should NEVER be the case
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
                die("Error: arg to --seed-length must be at least 20")
            if self.mate_inner_dist_std_dev != None and self.mate_inner_dist_std_dev < 0:
                die("Error: arg to --mate-std-dev must at least 0")
            if (not self.read_group_id and self.sample_id) or (self.read_group_id and not self.sample_id):
                die("Error: --rg-id and --rg-sample must be specified or omitted together")
    
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
                die("Error: arg to --min-closure-exon must be at least 20")
            if self.min_closure_intron_length < 0:
                die("Error: arg to --min-closure-intron must be at least 20")
            if self.max_closure_intron_length < 0:
                die("Error: arg to --max-closure-intron must be at least 20")
            if self.min_coverage_intron_length < 0:
                die("Error: arg to --min-coverage-intron must be at least 20")
            if self.max_coverage_intron_length < 0:
                die("Error: arg to --max-coverage-intron must be at least 20")
            if self.min_segment_intron_length < 0:
                die("Error: arg to --min-segment-intron must be at least 20")
            if self.max_segment_intron_length < 0:
                die("Error: arg to --max-segment-intron must be at least 20")
                    
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
        
        self.system_params = self.SystemParams(1,               # bowtie_threads (num_cpus)
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
        self.find_novel_indels = False
        self.find_GFF_juncs = True
        self.skip_check_reads = False
        self.max_hits = 20
        self.segment_length = 25
        self.segment_mismatches = 2
        self.bowtie_alignment_option = "-v"
        self.max_insertion_length = 3
        self.max_deletion_length = 3
        self.raw_insertions = None
        self.raw_deletions = None
        self.closure_search = None
        self.coverage_search = None
        self.microexon_search = False
        self.butterfly_search = None
        
    def check(self):
        self.splice_constraints.check()
        self.read_params.check()
        self.system_params.check()
       
        if self.segment_length <= 4:
            die("Error: arg to --segment-length must at least 4")
        if self.segment_mismatches < 0 or self.segment_mismatches > 3:
            die("Error: arg to --segment-mismatches must in [0, 3]")

        if self.read_params.color and self.butterfly_search:
            die("Error: butterfly-search in colorspace is not yet supported")

        library_types = ["fr-unstranded", "fr-firststrand", "fr-secondstrand"]

        if self.read_params.library_type != "" and self.read_params.library_type not in library_types:
            die("Error: library-type should be one of: "+', '.join(library_types))
        
        self.search_params.max_closure_intron_length = min(self.splice_constraints.max_intron_length,
                                                           self.search_params.max_closure_intron_length)
        
        self.search_params.max_segment_intron_length = min(self.splice_constraints.max_intron_length,
                                                           self.search_params.max_segment_intron_length)

        self.search_params.max_coverage_intron_length = min(self.splice_constraints.max_intron_length,
                                                            self.search_params.max_coverage_intron_length)
        
        if self.max_insertion_length >= self.segment_length:
            die("Error: the max insertion length ("+self.max_insertion_length+") can not be equal to or greater than the segment length ("+self.segment_length+")")

        if self.max_insertion_length < 0:
            die("Error: the max insertion length ("+self.max_insertion_length+") can not be less than 0")

        if self.max_deletion_length >= self.splice_constraints.min_intron_length:
            die("Error: the max deletion length ("+self.max_deletion_length+") can not be equal to or greater than the min intron length ("+self.splice_constraints.min_intron_length+")")

        if self.max_deletion_length < 0:
           die("Error: the max deletion length ("+self.max_deletion_length+") can not be less than 0")

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
               "--sam-header", sam_header,
               "--max-insertion-length", str(self.max_insertion_length),
               "--max-deletion-length", str(self.max_deletion_length)]
        
        cmd.extend(self.system_params.cmd())
               
        if self.read_params.mate_inner_dist != None:
            cmd.extend(["--inner-dist-mean", str(self.read_params.mate_inner_dist),
                        "--inner-dist-std-dev", str(self.read_params.mate_inner_dist_std_dev)])
        if self.gff_annotation != None:
            cmd.extend(["--gtf-annotations", str(self.gff_annotation)])
        if not self.closure_search:
            cmd.append("--no-closure-search")
        if not self.coverage_search:
            cmd.append("--no-coverage-search")
        if not self.microexon_search:
            cmd.append("--no-microexon-search")            
        if self.butterfly_search:
            cmd.append("--butterfly-search")
        if self.read_params.solexa_quals:
            cmd.append("--solexa-quals")
        if self.read_params.quals:
            cmd.append("--quals")
        if self.read_params.integer_quals:
            cmd.append("--integer-quals")
        if self.read_params.color:
            cmd.append("--color")
            if self.read_params.color_out:
                cmd.append("--color-out")
        if self.read_params.library_type != "":
            cmd.extend(["--library-type", self.read_params.library_type])
        if self.read_params.phred64_quals:
            cmd.append("--phred64-quals")
        return cmd
    
    # This is the master options parsing routine, which calls parse_options for
    # the delegate classes (e.g. SpliceConstraints) that handle certain groups
    # of options.
    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "hvp:m:F:a:i:I:G:r:o:j:z:g:QC", 
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
                                         "allow-indels",
                                         "no-novel-indels",
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
                                         "bowtie-n",
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
                                         "tmp-dir=",
                                         "zpacker=",
                                         "max-insertion-length=",
                                         "min-insertion-length=",
                                         "insertions=",
                                         "deletions="])
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
            if option in ("-g", "--max-multihits"):
                self.max_hits = int(value)
            if option in ("-G", "--GTF"):
                self.gff_annotation = value
            if option in ("-j", "--raw-juncs"):
                self.raw_junctions = value
            if option == "--no-novel-juncs":
                self.find_novel_juncs = False
            if option == "--allow-indels":
                self.find_novel_indels = True
            if option == "--no-novel-indels":
                self.find_novel_indels = False
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
            if option == "--bowtie-n":
                self.bowtie_alignment_option = "-n"
            if option == "--max-insertion-length":
                self.max_insertion_length = int(value)
            if option == "--max-deletion-length": 
                self.max_deletion_length = int(value)
            if option == "--insertions":
                self.raw_insertions = value
            if option == "--deletions":
                self.raw_deletions = value 
            if option in ("-o", "--output-dir"):
                custom_out_dir = value + "/"
            if option == "--tmp-dir":
                custom_tmp_dir = value + "/"
                
        if custom_out_dir != None:
            output_dir = custom_out_dir
            logging_dir = output_dir + "logs/"
            tmp_dir = output_dir + "tmp/"
            sam_header = tmp_dir + "stub_header.sam" 
        if custom_tmp_dir != None:
            tmp_dir = custom_tmp_dir
            sam_header = tmp_dir + "stub_header.sam" 
        if len(args) < 2:
            raise Usage(use_message)
        return args


def getFileDir(filepath):
   #if fullpath, returns path including the ending /
   fpath, fname=os.path.split(filepath)
   if fpath: fpath+='/'
   return fpath

def getFileBaseName(filepath):
   fpath, fname=os.path.split(filepath)
   fbase, fext =os.path.splitext(fname)
   fx=fext.lower()
   if (fx in ['.fq','.txt','.seq','.bwtout'] or fx.find('.fa')==0) and len(fbase)>0:
      return fbase
   elif fx == '.z' or fx.find('.gz')==0 or fx.find('.bz')==0:
      fb, fext = os.path.splitext(fbase)
      fx=fext.lower()
      if (fx in ['.fq','.txt','.seq','.bwtout'] or fx.find('.fa')==0) and len(fb)>0:
         return fb
      else:
         return fbase
   else:
     if len(fbase)>len(fext):
        return fbase
     else:
        return fname
   
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
        try:
          os.makedirs(tmp_dir)
        except OSError, o:
          die("\nError creating directory %s (%s)" % (tmp_dir, o))

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
            die("Error: Could not find Bowtie index files " + idx_prefix + ".*")
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
            die("Error: Could not find Bowtie index files " + idx_prefix + ".*")

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

        print >> sys.stderr, "Executing: " + " ".join(inspect_cmd) + " > " + tmp_fasta_file_name
        ret = subprocess.call(inspect_cmd, 
                              stdout=tmp_fasta_file,
                              stderr=inspect_log)

        # Bowtie reported an error
        if ret != 0:
           die(fail_str+"Error: bowtie-inspect returned an error")
           
    # Bowtie not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            die(fail_str+"Error: bowtie-inspect not found on this system.  Did you forget to include it in your PATH?")
  
    return tmp_fasta_file_name

# Checks whether the multifasta file for the genome is present alongside the 
# Bowtie index files for it.
def check_fasta(idx_prefix):
    print >> sys.stderr, "[%s] Checking for reference FASTA file" % right_now()
    idx_fasta = idx_prefix + ".fa"
    if os.path.exists(idx_fasta):
        return idx_fasta
    else:
        bowtie_idx_env_var = os.environ.get("BOWTIE_INDEXES")
        if bowtie_idx_env_var != None:
            idx_fasta = bowtie_idx_env_var + idx_prefix + ".fa" 
            if os.path.exists(idx_fasta):
                return idx_fasta
        
        print >> sys.stderr, "\tWarning: Could not find FASTA file " + idx_fasta
        idx_fa = bowtie_idx_to_fa(idx_prefix)
        return idx_fa
    
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
       errmsg=fail_str+str(o)+"\n"
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           errmsg+="Error: bowtie not found on this system"
       die(errmsg)
       
def get_index_sam_header(read_params, idx_prefix):
    try:
        bowtie_sam_header_filename = tmp_name()
        bowtie_sam_header_file = open(bowtie_sam_header_filename,"w")

        
        bowtie_header_cmd = ['bowtie', '--sam']
        if read_params.color:
            bowtie_header_cmd.append('-C')
        bowtie_header_cmd.extend([idx_prefix, '/dev/null'])
        subprocess.call(bowtie_header_cmd,stdout=bowtie_sam_header_file, stderr=open('/dev/null'))

        bowtie_sam_header_file.close()
        bowtie_sam_header_file = open(bowtie_sam_header_filename,"r")
        
        sam_header_file = open(sam_header, "w")
        
        preamble = []
        sq_dict_lines = []
        
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
                    die("Error: malformed sequence dictionary in sam header")
                sq_dict_lines.append([seq_name,line])
            elif line.find("CL"):
                continue
            else:
                preamble.append(line)
        #for line in preamble:
        #    print >> sam_header_file, line

#       print >> sam_header_file, "@HD\tVN:1.0\tSO:sorted"
        print >> sam_header_file, "@HD\tVN:1.0\tSO:coordinate"
    
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
       errmsg=fail_str+str(o)+"\n"
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           errmsg+="Error: bowtie not found on this system"
       die(errmsg)

# Make sure Bowtie is installed and is recent enough to be useful
def check_bowtie():
    print >> sys.stderr, "[%s] Checking for Bowtie" % right_now()
    bowtie_version = get_bowtie_version()
    if bowtie_version == None:
        die("Error: Bowtie not found on this system")
    # daehwan - check
    elif bowtie_version[1] < 12 or bowtie_version[2] < 3:
        die("Error: TopHat requires Bowtie 0.12.3 or later")
    print >> sys.stderr, "\tBowtie version:\t\t\t %s" % ".".join([str(x) for x in bowtie_version])
    

# Retrive a tuple containing the system's version of samtools.  Parsed from 
# `samtools`
def get_samtools_version():
    try:
        # Launch Bowtie to capture its version info
        proc = subprocess.Popen(['samtools'],stderr=subprocess.PIPE)
        samtools_out = proc.communicate()[1]

        # Find the version identifier
        version_match = re.search(r'Version:\s+(\d+)\.(\d+).(\d+)([a-zA-Z]?)', samtools_out)
        samtools_version_arr = [int(version_match.group(x)) for x in [1,2,3]]
        if version_match.group(4):
            samtools_version_arr.append(version_match.group(4))
        else:
            samtools_version_arr.append(0)
            
        return version_match.group(), samtools_version_arr
    except OSError, o:
       errmsg=fail_str+str(o)+"\n"
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           errmsg+="Error: samtools not found on this system"
       die(errmsg)

# Make sure the SAM tools are installed and are recent enough to be useful
def check_samtools():
    print >> sys.stderr, "[%s] Checking for Samtools" % right_now()
    samtools_version_str, samtools_version_arr = get_samtools_version()
    if samtools_version_str == None:
        die("Error: Samtools not found on this system")
    elif  samtools_version_arr[1] < 1 or samtools_version_arr[2] < 7:
        die("Error: TopHat requires Samtools 0.1.7 or later")
    print >> sys.stderr, "\tSamtools %s" % samtools_version_str
      


class FastxReader:
  def __init__(self, i_file, is_color=0, fname=''):
    self.bufline=None
    self.format=None
    self.ifile=i_file
    self.nextRecord=None
    self.eof=None
    self.fname=fname
    self.lastline=None
    self.numrecords=0
    self.isColor=0
    if is_color : self.isColor=1
    # determine file type
    #no records processed yet, skip custom header lines if any
    hlines=10 # allow maximum 10 header lines
    self.lastline=" "
    while hlines>0 and self.lastline[0] not in "@>" :
       self.lastline=self.ifile.readline()
       hlines-=1
    if self.lastline[0] == '@':
      self.format='fastq'
      self.nextRecord=self.nextFastq
    elif self.lastline[0] == '>':
      self.format='fasta'
      self.nextRecord=self.nextFasta
    else: 
      die("Error: cannot determine record type in input file %s" % fname)
    self.bufline=self.lastline
    self.lastline=None
    
  def nextFastq(self):
    # returning tuple: (seqID, sequence_string, seq_len, qv_string)
    seqid,seqstr,qstr,seq_len='','','',0
    if self.eof: return (seqid, seqstr, seq_len, qstr)
    fline=self.getLine #shortcut to save a bit of time
    line=fline()

    if not line : return (seqid, seqstr, seq_len, qstr)
    while len(line.rstrip())==0: # skip empty lines
      line=fline()
      if not line : return (seqid, seqstr,seq_len, qstr)
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
             # if qtitle and qtitle != seqid:
             #   raise ValueError("Different read ID for sequence and quality (%s vs %s)" \
             #                    % (seqid, qtitle))
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
          if qstrlen + self.isColor >= seq_len :
               break # qv string has reached the length of seq string
          #loop until qv has the same length as seq
      if seq_len != qstrlen+self.isColor :
           raise ValueError("Length mismatch between sequence and quality strings "+ \
                                "for %s (%i vs %i)." \
                                % (seqid, seq_len, qstrlen))
    except ValueError, err:
        die("\nError encountered parsing file "+self.fname+":\n "+str(err))
    #return the record 
    self.numrecords+=1
    if self.isColor :
        seq_len-=1
        seqstr = seqstr[1:]
    return (seqid, seqstr, seq_len, qstr)

  def nextFasta(self):
    # returning tuple: (seqID, sequence_string, seq_len)
    seqid,seqstr,seq_len='','',0
    fline=self.getLine # shortcut to readline function of f
    line=fline() # this will use the buffer line if it's there
    if not line : return (seqid, seqstr, seq_len, None)
    while len(line.rstrip())==0: # skip empty lines
      line=fline()
      if not line : return (seqid, seqstr, seq_len, None)
    try:
       if line[0] != ">":
          raise ValueError("Records in Fasta files must start with '>' character")
       seqid = line[1:].split()[0]
       #more sequence lines, or the ">" quality marker line:
       while True:
          line = fline()
          if not line: break
          if line[0] == '>':
             #next sequence starts here  
             self.ungetLine()
             break
          seqstr += line.rstrip()
          #loop until '>' found
       seq_len = len(seqstr)
       if seq_len < 3:
          raise ValueError("Read %s too short (%i)." \
                           % (seqid, seq_len))
    except ValueError, err:
        die("\nError encountered parsing fasta file "+self.fname+"\n "+str(err))
    #return the record and continue
    self.numrecords+=1
    if self.isColor :
        seq_len-=1
        seqstr = seqstr[1:]
    return (seqid, seqstr, seq_len, None)

  def getLine(self):
      if self.bufline: #return previously buffered line
         r=self.bufline
         self.bufline=None
         return r
      else: #read a new line from stream and return it
         if self.eof: return None
         self.lastline=self.ifile.readline()
         if not self.lastline:
            self.eof=1
            return None
         return self.lastline
  def ungetLine(self):
      if self.lastline is None:
         print >> sys.stderr, "Warning: FastxReader called ungetLine() with no prior line!"
      self.bufline=self.lastline
      self.lastline=None
#< class FastxReader

class ZReader:
    def __init__(self, filename, sysparams, guess=True):
        self.fname=filename
        self.file=None
        self.fsrc=None 
        self.popen=None
        pipecmd=[]
        if guess:
           s=filename.lower()
           rgz=s.rfind(".gz")
           if rgz>0 and rgz>len(s)-6:
                  pipecmd=['gzip']
           else:
                rgz=s.rfind(".z")
                if rgz>0 and rgz==len(s)-3:
                   pipecmd=['gzip']
                else:
                   rgz=s.rfind(".bz")
                   if rgz>0 and rgz>len(s)-7:
                       pipecmd=['bzip2']
           if len(pipecmd)>0 and which(pipecmd[0]) is None:
               die("Error: cannot find %s to decompress input file %s " % (pipecmd, filename))

           if len(pipecmd)>0:
              if pipecmd[0]=='gzip' and sysparams.zipper.endswith('pigz'):
                 pipecmd[0]=sysparams.zipper
                 pipecmd.extend(sysparams.zipper_opts)
              elif pipecmd[0]=='bzip2' and sysparams.zipper.endswith('pbzip2'):
                 pipecmd[0]=sysparams.zipper
                 pipecmd.extend(sysparams.zipper_opts)
        else:
           pipecmd=[sysparams.zipper]
           pipecmd.extend(sysparams.zipper_opts)

        if pipecmd:
           pipecmd+=['-cd']
           try:
              self.fsrc=open(self.fname, 'rb')
              self.popen=subprocess.Popen(pipecmd, 
                    stdin=self.fsrc, 
                    stdout=subprocess.PIPE, close_fds=True)
           except Exception:
              die("Error: could not open pipe "+' '.join(pipecmd)+' < '+ self.fname)
           self.file=self.popen.stdout
        else: 
           self.file=open(filename)
    def close(self):
       if self.fsrc: self.fsrc.close()
       self.file.close()
       if self.popen: self.popen.wait()
       self.popen=None

class ZWriter:
   def __init__(self, filename, sysparams):
      pipecmd=[sysparams.zipper,'-cf', '-']
      self.ftarget=open(filename, 'wb')
      try:
         self.popen=subprocess.Popen(pipecmd, 
               stdin=subprocess.PIPE, 
               stdout=self.ftarget, close_fds=True)
      except Exception:
          die("Error: could not open writer pipe "+' '.join(pipecmd)+' < '+ self.fname)
      self.fname=filename
      self.file=self.popen.stdin # client writes to this end of the pipe
      
   def close(self):
      self.file.close()
      self.ftarget.close()
      self.popen.wait() #! required to actually flush the pipes (grumble)
      self.popen=None

# check_reads() has several jobs.  It examines the user's reads, one file at a 
# time, and determnes the file format, read length, and other properties that 
# are used to set the junction search strategy later on.  
# TODO: When we add support for mixed read lengths, this routine 
# will need to set the seed length differently. 
def check_reads(params, reads_files):
    print >> sys.stderr, "[%s] Checking reads" % right_now()
    get_bowtie_version()
    
    seed_len = params.read_params.seed_length
    fileformat = params.read_params.reads_format

    observed_formats = set([])
    # observed_scales = set([])
    min_seed_len = 99999
    max_seed_len = 0
    max_qual = -1
    files = reads_files.split(',')

    for f_name in files:
        #try:
        zf = ZReader(f_name, params.system_params)
        #except IOError:
        #   die("Error: could not open file "+f_name)
        freader=FastxReader(zf.file, params.read_params.color, zf.fname)
        while True:
            seqid, seqstr, seq_len, qstr = freader.nextRecord()
            if not seqid: break
            if seq_len < 20:
                  print >> sys.stderr, "Warning: found a read < 20bp in", f_name
            else:
                min_seed_len = min(seq_len, min_seed_len)
                max_seed_len = max(seq_len, max_seed_len)
            if freader.format == "fastq":
                max_line_qual = max([ord(x) for x in list(qstr)])
                max_qual = max(max_line_qual, max_qual)
        zf.close()
        observed_formats.add(freader.format)
    if len(observed_formats) > 1:
        die("Error: TopHat requires all reads be either FASTQ or FASTA.  Mixing formats is not supported.")
    fileformat=list(observed_formats)[0]
    if seed_len != None:
        seed_len = max(seed_len, min_seed_len)
    else:
        seed_len = max_seed_len
        
    print >> sys.stderr, "\tmin read length: %dbp, max read length: %dbp" % (min_seed_len, max_seed_len)
    print >> sys.stderr, "\tformat:\t\t %s" % fileformat
    if fileformat == "fastq":
        quality_scale = "phred33 (default)"
        if params.read_params.solexa_quals and not params.read_params.phred64_quals:
            quality_scale = "solexa33 (reads generated with GA pipeline version < 1.3)"
        elif params.read_params.phred64_quals:
            quality_scale = "phred64 (reads generated with GA pipeline version >= 1.3)"
        print >> sys.stderr, "\tquality scale:\t %s" % quality_scale
    elif fileformat == "fasta":
        if params.read_params.color:
            params.integer_quals = True
    
    #print seed_len, format, solexa_scale
    return TopHatParams.ReadParams(params.read_params.solexa_quals,
                                   params.read_params.phred64_quals,
                                   params.read_params.quals,
                                   params.read_params.integer_quals,
                                   params.read_params.color,
                                   params.read_params.color_out,
                                   params.read_params.library_type,
                                   seed_len, 
                                   fileformat, 
                                   params.read_params.mate_inner_dist, 
                                   params.read_params.mate_inner_dist_std_dev,
                                   params.read_params.read_group_id,
                                   params.read_params.sample_id,
                                   params.read_params.library_id,
                                   params.read_params.description,
                                   params.read_params.seq_platform_unit,
                                   params.read_params.seq_center,
                                   params.read_params.seq_run_date,
                                   params.read_params.seq_platform)

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
    #reads_suffix = ".fq"
    reads_suffix = ".fq.z"
    kept_reads_filename = tmp_dir + output_name + reads_suffix
    
    if os.path.exists(kept_reads_filename):
        os.remove(kept_reads_filename)
    kept_reads = open(kept_reads_filename, "wb")
    
    filter_log = open(logging_dir + "prep_reads.log", "w")
    
    filter_cmd = [prog_path("prep_reads")]
    filter_cmd.extend(params.cmd())
    if params.read_params.reads_format == "fastq":
        filter_cmd += ["--fastq"]
    elif params.read_params.reads_format == "fasta":
        filter_cmd += ["--fasta"]

    filter_cmd.append(reads_list)

    if params.read_params.quals:
        filter_cmd.append(quals_list)
       
    #finally, add the compression pipe
    zip_cmd=[ params.system_params.zipper ]
    zip_cmd.extend(params.system_params.zipper_opts)
    zip_cmd.extend(['-c','-'])
    shell_cmd = ' '.join(filter_cmd)+' | '+' '.join(zip_cmd) +' >' +kept_reads_filename
    try:       
        print >> run_log, shell_cmd
        
        filter_proc = subprocess.Popen(filter_cmd, 
                              stdout=subprocess.PIPE,
                              stderr=filter_log)
        zip_proc=subprocess.Popen(zip_cmd,
                              stdin=filter_proc.stdout,
                              stdout=kept_reads)
        filter_proc.stdout.close() #as per http://bugs.python.org/issue7678
        zip_proc.communicate()
    # prep_reads not found
    except OSError, o:
        errmsg=fail_str+str(o)
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            errmsg+="\nprep_reads not found on this system.  Did you forget to include it in your PATH?"
        die("Error: "+errmsg)
    kept_reads.close()
    return kept_reads_filename

# Call bowtie
def bowtie(params,
           bwt_idx_prefix,
           reads_list,
           reads_format,
           mapped_reads,
           unmapped_reads,
           reads_for_ordering = None,
           extra_output = ""):
    start_time = datetime.now()
    bwt_idx_name = bwt_idx_prefix.split('/')[-1]
    print >> sys.stderr, "[%s] Mapping reads against %s with Bowtie %s" % (start_time.strftime("%c"), bwt_idx_name, extra_output)
    
    # Setup Bowtie output redirects
    #bwt_map = output_dir + mapped_reads
    #bwt_map = tmp_name()
    #tmp_fname = bwt_map.split('/')[-1]
    readfile_basename=getFileBaseName(reads_list)
    #tmp_fname = 'bowtie.'+reads_list+'.fixmap.log'
    bwt_log = open(logging_dir + 'bowtie.'+readfile_basename+'.fixmap.log', "w")
    bwt_mapped=mapped_reads
    # Launch Bowtie
    try:    
        bowtie_cmd = ["bowtie"]
        
        if reads_format == "fastq":
            bowtie_cmd += ["-q"]
        elif reads_format == "fasta":
            bowtie_cmd += ["-f"]

        if params.read_params.color:
            bowtie_cmd += ["-C", "--col-keepends"]


        unzip_cmd=[ params.system_params.zipper ]
        unzip_cmd.extend(params.system_params.zipper_opts)
        zip_cmd=unzip_cmd[:]
        unzip_cmd+=['-cd']
        zip_cmd+=['-c']

        global unmapped_reads_fifo
           
        if unmapped_reads != None:
             bowtie_cmd += ["--un", unmapped_reads_fifo,
                           "--max", "/dev/null"]
             if os.fork()==0:
                subprocess.call(zip_cmd, 
                                stdin=open(unmapped_reads_fifo, "r"), 
                                stdout=open(unmapped_reads, "w"))
                os._exit(os.EX_OK)

        fix_map_cmd = [prog_path('fix_map_ordering')]

        unzip_proc = subprocess.Popen(unzip_cmd, 
                                 stdin=open(reads_list, "rb"),
                                 stdout=subprocess.PIPE)
        bowtie_cmd += [params.bowtie_alignment_option, str(params.segment_mismatches),
                         "-p", str(params.system_params.num_cpus),
                         "-k", str(params.max_hits),
                         "-m", str(params.max_hits),
                         bwt_idx_prefix, 
                         '-']
        bowtie_proc = subprocess.Popen(bowtie_cmd, 
                                 stdin=unzip_proc.stdout, 
                                 stdout=subprocess.PIPE, 
                                 stderr=bwt_log)
        unzip_proc.stdout.close() # see http://bugs.python.org/issue7678
        
        if reads_format == "fastq":
            fix_map_cmd += ["--fastq"]
        elif reads_format == "fasta":
            fix_map_cmd += ["--fasta"]
        
        if reads_for_ordering == None:
            reads_for_ordering = reads_list
        
        fix_map_cmd.extend(['-',reads_for_ordering])
        shell_cmd=' '.join(unzip_cmd) + "< '" +reads_list +"'|" + ' '.join(bowtie_cmd) + '|' + \
            ' '.join(fix_map_cmd) + "|"+ ' '.join(zip_cmd)+" >'" + bwt_mapped+"'"
        fix_order_proc = subprocess.Popen(fix_map_cmd, 
                                          stdin=bowtie_proc.stdout,
                                          stdout=subprocess.PIPE)
        bowtie_proc.stdout.close()
        zip_proc = subprocess.Popen(zip_cmd,
                                 stdin=fix_order_proc.stdout,
                                 stdout=open(mapped_reads, "wb"))
        fix_order_proc.stdout.close()
        print >> run_log, shell_cmd
        #subprocess.call(shell_cmd, shell=True)
        #os.system(shell_cmd)
        zip_proc.communicate()
        
    except OSError, o:
        die(fail_str+"Error: "+str(o))
            
    # Success    
    #finish_time = datetime.now()
    #duration = finish_time - start_time
    #print >> sys.stderr, "\t\t\t[%s elapsed]" %  formatTD(duration)
    return (bwt_mapped, unmapped_reads)

def bowtie_segment(params,
           bwt_idx_prefix,
           reads_list,
           reads_format,
           mapped_reads,
           unmapped_reads,
           reads_for_ordering = None,
           extra_output = ""):

    backup_bowtie_alignment_option = params.bowtie_alignment_option
    params.bowtie_alignment_option = "-v"
    params.max_hits *= 2
    
    result = bowtie(params, bwt_idx_prefix, reads_list, reads_format,
                    mapped_reads, unmapped_reads, reads_for_ordering,
                    extra_output)

    params.bowtie_alignment_option = backup_bowtie_alignment_option
    params.max_hits /= 2
    return result
    
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
            die(fail_str+"Error: GTF junction extraction failed with err ="+str(retcode))
    # cvg_islands not found
    except OSError, o:
       errmsg=fail_str+str(o)+"\n"
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           errmsg+="Error: gtf_juncs not found on this system"
       die(errmsg)
    return (True, gtf_juncs_out_name)

# Call bowtie-build on the FASTA file of sythetic splice junction sequences
def build_juncs_bwt_index(external_splice_prefix, color):
    print >> sys.stderr, "[%s] Indexing splices" % (right_now())
    bowtie_build_log = open(logging_dir + "bowtie_build.log", "w")
    
    #user_splices_out_prefix  = output_dir + "user_splices_idx"
    
    bowtie_build_cmd = ["bowtie-build"]
    if color:
        bowtie_build_cmd += ["-C"]
        
    bowtie_build_cmd += [external_splice_prefix + ".fa",
                         external_splice_prefix]            
    try:    
        print >> run_log, " ".join(bowtie_build_cmd)
        retcode = subprocess.call(bowtie_build_cmd, 
                                 stdout=bowtie_build_log)
       
        if retcode != 0:
            die(fail_str+"Error: Splice sequence indexing failed with err ="+ str(retcode))
    except OSError, o:
        errmsg=fail_str+str(o)+"\n"
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            errmsg+="Error: bowtie-build not found on this system"
        die(errmsg)
    return external_splice_prefix

# Build a splice index from a .juncs file, suitable for use with specified read
# (or read segment) lengths
def build_juncs_index(min_anchor_length, 
                      read_length,
                      juncs_prefix, 
                      external_juncs,
                      external_insertions,
                      external_deletions,  
                      reference_fasta,
                      color):
    print >> sys.stderr, "[%s] Retrieving sequences for splices" % (right_now())
    
    juncs_file_list = ",".join(external_juncs)
    insertions_file_list = ",".join(external_insertions)
    deletions_file_list = ",".join(external_deletions)


    juncs_db_log = open(logging_dir + "juncs_db.log", "w")
    
    external_splices_out_prefix  = tmp_dir + juncs_prefix
    external_splices_out_name = external_splices_out_prefix + ".fa"
    
    external_splices_out = open(external_splices_out_name, "w")
    # juncs_db_cmd = [bin_dir + "juncs_db",
    juncs_db_cmd = [prog_path("juncs_db"), 
                    str(min_anchor_length),
                    str(read_length),
                    juncs_file_list,
                    insertions_file_list,
                    deletions_file_list,
                    reference_fasta]            
    try:    
        print >> run_log, " ".join(juncs_db_cmd)
        retcode = subprocess.call(juncs_db_cmd, 
                                 stderr=juncs_db_log,
                                 stdout=external_splices_out)
       
        if retcode != 0:
            die(fail_str+"Error: Splice sequence retrieval failed with err ="+str(retcode))
    # juncs_db not found
    except OSError, o:
       errmsg=fail_str+str(o)+"\n"
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           errmsg+="Error: juncs_db not found on this system"
       die(errmsg)
       
    external_splices_out_prefix = build_juncs_bwt_index(external_splices_out_prefix, color)
    return external_splices_out_prefix

# Print out the sam header, embedding the user's specified library properties.
# FIXME: also needs SQ dictionary lines
def write_sam_header(read_params, sam_file):
    #print >> sam_file, "@HD\tVN:1.0\tSO:sorted"
    print >> sam_file, "@HD\tVN:1.0\tSO:coordinate"
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
    
    left_maps = [x for x in left_maps if (os.path.exists(x) and os.path.getsize(x) > 25)]
    left_maps = ','.join(left_maps)
    
    if len(right_maps) > 0:
        right_maps = [x for x in right_maps if (os.path.exists(x) and os.path.getsize(x) > 25)]
        right_maps = ','.join(right_maps)
    
    report_log = open(logging_dir + "reports.log", "w")
    junctions = output_dir + "junctions.bed"
    insertions = output_dir + "insertions.bed"
    deletions = output_dir + "deletions.bed" 
    coverage =  "coverage.wig"
    accepted_hits_sam = tmp_dir + "accepted_hits.sam"
    report_cmdpath = prog_path("tophat_reports")
    report_cmd = [report_cmdpath]
    report_cmd.extend(params.cmd())
        
    report_cmd.extend([junctions,
                       insertions,
                       deletions,
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
       
        if retcode != 0:
            die(fail_str+"Error: Report generation failed with err ="+str(retcode))
        
        tmp_bam = tmp_name()+".bam"
        sam_to_bam_cmd = ["samtools", "view", "-S", "-b", accepted_hits_sam]
        print >> run_log, " ".join(sam_to_bam_cmd) + " > " + tmp_bam
        sam_to_bam_log = open(logging_dir + "accepted_hits_sam_to_bam.log", "w")
        tmp_bam_file = open(tmp_bam, "w")
        ret = subprocess.call(sam_to_bam_cmd, 
                              stdout=tmp_bam_file,
                              stderr=sam_to_bam_log)
        if ret != 0:
            die("Error: could not convert to BAM with samtools")
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
            die("Error: could not sort BAM file with samtools")
            
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
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: tophat_reports not found on this system"
        die(str(o))
    return (coverage, junctions)


# Split up each read in a FASTQ file into multiple segments. Creates a FASTQ file
# for each segment  This function needs to be fixed to support mixed read length
# inputs

def open_output_files(prefix, num_files_prev, num_files, out_segf, extension, params):
       i = num_files_prev + 1
       while i <= num_files:
          #out_fhandles.append(open(prefix + ("_seg%d" % i) + extension, "w"))
          segfname=prefix+("_seg%d" % i)+extension
          out_segf.append(ZWriter(segfname,params.system_params));
          i += 1

def split_reads(reads_filename, 
                prefix,
                fasta,
                params,
                segment_length):
    #reads_file = open(reads_filename)
    zreads = ZReader(reads_filename, params.system_params, False)
    out_segfiles = []
    
    if fasta:
        extension = ".fa.z"
    else:
        extension = ".fq.z"

    def convert_color_to_bp(color_seq):
        decode_dic = { 'A0':'A', 'A1':'C', 'A2':'G', 'A3':'T', 'A4':'N', 'A.':'N', 'AN':'N',
                       'C0':'C', 'C1':'A', 'C2':'T', 'C3':'G', 'C4':'N', 'C.':'N', 'CN':'N',
                       'G0':'G', 'G1':'T', 'G2':'A', 'G3':'C', 'G4':'N', 'G.':'N', 'GN':'N',
                       'T0':'T', 'T1':'G', 'T2':'C', 'T3':'A', 'T4':'N', 'T.':'N', 'TN':'N',
                       'N0':'N', 'N1':'N', 'N2':'N', 'N3':'N', 'N4':'N', 'N.':'N', 'NN':'N',
                       '.0':'N', '.1':'N', '.2':'N', '.3':'N', '.4':'N', '..':'N', '.N':'N' }

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

    def split_record(read_name, read_seq, read_qual, out_segf, offsets, color):
        if color:
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
            f = out_segf[seg_num].file
            seg_seq = read_seq[last_seq_offset+color_offset:offsets[seg_num + 1]+color_offset]
            print >> f, "%s|%d:%d:%d" % (read_name,last_seq_offset,seg_num, len(offsets) - 1)
            if color:
                print >> f, "%s%s" % (read_seq_temp[last_seq_offset], seg_seq)
            else:
                print >> f, seg_seq
            if not fasta:
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
    for line in zreads.file:
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
                open_output_files(prefix, num_segments, tmp_num_segments, out_segfiles, extension, params)
                num_segments = tmp_num_segments

            if fasta:
                split_record(read_name, read_seq, None, out_segfiles, offsets, params.read_params.color)
        elif line_state == 2:
            line = line.strip()
        else:
            read_quals = line.strip()
            if not fasta:
                split_record(read_name, read_seq, read_quals, out_segfiles, offsets, params.read_params.color)
                
        line_state += 1
        if fasta:
            line_state %= 2
        else:
            line_state %= 4
    zreads.close()
    out_fnames=[]
    for zf in out_segfiles:
        zf.close()
        out_fnames.append(zf.fname)
    #return [o.fname for o in out_segfiles]
    return out_fnames

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

    left_maps = ','.join(left_maps)
    right_maps = ','.join(right_maps)

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
           die(fail_str+"Error: closure-based junction search failed with err ="+str(retcode))
    # cvg_islands not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: closure_juncs not found on this system"
        die(str(o))
    return [juncs_out]

# Find possible junctions by examining coverage and split segments in the initial
# map and segment maps.  Report junctions, insertions, and deletions in segment.juncs,
# segment.insertions, and segment.deletions.  Calls the executable
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
    #slash = left_seg_maps[0].rfind('/')
    #juncs_out = ""
    #if slash != -1:
    #    juncs_out += left_seg_maps[0][:slash+1]
    out_path=getFileDir(left_seg_maps[0])
    juncs_out=out_path+"segment.juncs"
    #insertions_out = ""
    insertions_out=out_path+"segment.insertions"
    deletions_out =out_path+"segment.deletions"

    left_maps = ','.join(left_seg_maps)
    align_log = open(logging_dir + "segment_juncs.log", "w")
    align_cmd = [prog_path("segment_juncs")]
    
    align_cmd.extend(params.cmd())
    
    align_cmd.extend(["--ium-reads", ",".join(unmapped_reads),
                      ref_fasta,
                      juncs_out,
                      insertions_out,
                      deletions_out,
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
           die(fail_str+"Error: segment-based junction search failed with err ="+str(retcode))
    # cvg_islands not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: segment_juncs not found on this system"
        die(str(o))

    return [juncs_out, insertions_out, deletions_out]

# Joins mapped segments into full-length read alignments via the executable
# long_spanning_reads
def join_mapped_segments(params,
                         sam_header_filename,
                         reads,
                         ref_fasta,
                         possible_juncs,
                         possible_insertions,
                         possible_deletions,
                         contig_seg_maps,
                         spliced_seg_maps,
                         alignments_out_name):
    print >> sys.stderr, "[%s] Joining segment hits" % right_now()
    contig_seg_maps = ','.join(contig_seg_maps)

    possible_juncs = ','.join(possible_juncs)
    possible_insertions = ",".join(possible_insertions)
    possible_deletions = ",".join(possible_deletions)

    
    align_log = open(logging_dir + "long_spanning_reads.log", "w")
    align_cmd = [prog_path("long_spanning_reads")]
    
    ## if params.read_params.reads_format == "fastq":
    ##     align_cmd += ["-q"]
    ## elif params.read_params.reads_format == "fasta":
    ##     align_cmd += ["-f"]
    
    alignments_out = open(alignments_out_name, "w")
    ##write_sam_header(params.read_params, sorted_map)
    #zip_cmd=[ params.system_params.zipper, '-c', '-' ]
    
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
                       possible_insertions,
                       possible_deletions,
                       contig_seg_maps])
    if spliced_seg_maps != None:
        spliced_seg_maps = ','.join(spliced_seg_maps)
        align_cmd.append(spliced_seg_maps)
    #skip_Join=True  # DEBUG
    #if not skip_Join:
    try:    
        print >> run_log, " ".join(align_cmd),">", alignments_out_name
        retcode = subprocess.call(align_cmd, 
                                  stderr=align_log,
                                  stdout=alignments_out)

        # spanning_reads returned an error 
        if retcode != 0:
            die(fail_str+"Error: Segment join failed with err ="+str(retcode))
     # cvg_islands not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: long_spanning_reads not found on this system"
        die(str(o))
      

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
                      user_supplied_junctions,
                      user_supplied_insertions,
                      user_supplied_deletions):
    
    possible_juncs = []
    possible_juncs.extend(user_supplied_junctions)

    possible_insertions = []
    possible_insertions.extend(user_supplied_insertions)

    possible_deletions = []
    possible_deletions.extend(user_supplied_deletions)
    
    maps = { left_reads : [], right_reads : [] }
    #single_segments = False

    # Perform the first part of the TopHat work flow on the left and right 
    # reads of paired ends separately - we'll use the pairing information later
    for reads in [left_reads, right_reads]:
        if reads == None or os.path.getsize(reads)<25 :
            continue
        fbasename=getFileBaseName(reads)
#        slash = reads.rfind("/")
#        extension = reads.rfind(".")
#        if extension != -1:
#             prefix = reads[slash+1:extension]
#        else:
#             prefix = reads      
#           
#        extension = reads.rfind(".")
#        assert extension != -1
#        tmp = reads.rfind("/")
          
#        unspliced_out = tmp_dir + reads[tmp+1:extension] + ".bwtout.z"  
#        unmapped_unspliced = tmp_dir + reads[tmp+1:extension] + "_missing.fq.z"
        unspliced_out = tmp_dir + fbasename + ".bwtout.z"
        unmapped_unspliced = tmp_dir + fbasename + "_missing.fq.z"
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
        #unspliced_sam = tmp_name()
        unspliced_sam = tmp_name()+'.unspl.sam'
        # TODO: should write BAM directly

        # Convert the initial Bowtie maps into SAM format.  
        # TODO: this step should be removable now that Bowtie supports SAM 
        # format
        join_mapped_segments(params,
                             sam_header_filename,
                             reads,
                             ref_fasta,
                             ["/dev/null"],
			     ["/dev/null"],
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
            if not params.closure_search:
                params.closure_search = False
            if not params.coverage_search:
                params.coverage_search = False
            if not params.butterfly_search:
                params.butterfly_search = False
        if not params.closure_search:
            params.closure_search = False
        seg_maps = []
        unmapped_segs = []
        segs = []
        if num_segs > 1:
            # split up the IUM reads into segments
            read_segments = split_reads(unmapped_unspliced,
                                        tmp_dir + fbasename, 
                                        False,
                                        params,
                                        segment_len)

            # Map each segment file independently with Bowtie
            for i in range(len(read_segments)):
                seg = read_segments[i]
                #extension = seg.rfind(".")
                #assert extension != -1
                #tmp = seg.rfind("/")
                #assert tmp != -1
                #seg_out =  tmp_dir + seg[tmp+1:extension] + ".bwtout"
                #unmapped_seg = tmp_dir + seg[tmp+1:extension] + "_missing.fq"
                fbasename=getFileBaseName(seg)
                seg_out =  tmp_dir + fbasename + ".bwtout.z"
                unmapped_seg = tmp_dir + fbasename + "_missing.fq.z"
                extra_output = "(%d/%d)" % (i+1, len(read_segments))
                (seg_map, unmapped) = bowtie_segment(params,
                                                     bwt_idx_prefix, 
                                                     seg,
                                                     "fastq",
                                                     seg_out,
                                                     unmapped_seg,
                                                     seg,
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
    
    if params.find_novel_juncs or params.find_novel_indels:
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
        if params.find_novel_juncs:
                if os.path.getsize(juncs[0]) != 0:
                    possible_juncs.append(juncs[0])
        if params.find_novel_indels:
                if os.path.getsize(juncs[1]) != 0:
                    possible_insertions.append(juncs[1])
                if os.path.getsize(juncs[2]) != 0:
                    possible_deletions.append(juncs[2])
        # Optionally, and for paired reads only, use a closure search to 
        # discover addtional junctions
        if params.find_novel_juncs and params.closure_search and left_reads != None and right_reads != None:
            juncs = junctions_from_closures(params,
                                            [maps[left_reads].unspliced_bwt, maps[left_reads].seg_maps[-1]],
                                            [maps[right_reads].unspliced_bwt, maps[right_reads].seg_maps[-1]],
                                            ref_fasta)
            if os.path.getsize(juncs[0]) != 0:
                possible_juncs.extend(juncs)

    if len(possible_insertions) == 0 and len(possible_deletions) == 0 and len(possible_juncs) == 0:
        spliced_seg_maps = None
        junc_idx_prefix = None
    else:
        junc_idx_prefix = "segment_juncs"
    if len(possible_insertions) == 0:
        possible_insertions.append(os.devnull)
        # print >> sys.stderr, "Warning: insertions database is empty!"
    if len(possible_deletions) == 0:
        possible_deletions.append(os.devnull)
        # print >> sys.stderr, "Warning: deletions database is empty!"
    if len(possible_juncs) == 0:
        possible_juncs.append(os.devnull)
        print >> sys.stderr, "Warning: junction database is empty!"
    if junc_idx_prefix != None:  
        build_juncs_index(3,
                          segment_len,
                          junc_idx_prefix, 
                          possible_juncs,
                          possible_insertions,
                          possible_deletions,
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
                #fsegdir=getFileDir(seg)
                fsegname=getFileBaseName(seg)
                ordering = maps[reads].segs[i]
                seg_out = tmp_dir + fsegname + ".to_spliced.bwtout.z"
                (seg_map, unmapped) = bowtie_segment(params,
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
        rfname=getFileBaseName(reads)
        rfdir=getFileDir(reads)
        mapped_reads = rfdir+rfname+".candidate_hits.sam"
        join_mapped_segments(params,
                             sam_header_filename,
                             reads,
                             ref_fasta,
                             possible_juncs,
                             possible_insertions,
                             possible_deletions,
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

            subprocess.call(merge_sort_cmd,
                             stdout=open(merged_map,"w")) 
        else:
            merged_map = mapped_reads
        
        if not params.system_params.keep_tmp:
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
        progpath = os.path.join(bin_dir, program)
        if is_executable(progpath):
           return progpath
        for path in os.environ["PATH"].split(os.pathsep):
           progpath = os.path.join(path, program)
           if is_executable(progpath):
              return progpath
    return None

def die(msg=None):
 if msg is not None: 
    print >> sys.stderr, msg
    sys.exit(1)

def prog_path(program):
    progpath=which(program)
    if progpath == None:
        die("Error locating program: "+program)
    return progpath

# FIXME: this should get set during the make dist autotools phase of the build
def get_version():
   return "1.2.1"

def mlog(msg):
  print >> sys.stderr, "[DBGLOG]:"+msg

def test_input_file(filename):
    try:
        test_file = open(filename, "r")
    # Bowtie not found
    except IOError, o: 
        die("Error: Opening file %s" % filename)
    return
    
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
        if (not params.read_params.quals and len(args) > 2) or (params.read_params.quals and len(args) > 3):
            if params.read_params.mate_inner_dist == None:
                die("Error: you must set the mean inner distance between mates with -r")

            right_reads_list = args[2]
            if params.read_params.quals:
                left_quals_list = args[3]
                right_quals_list = args[4]
        else:
            right_reads_list = None
            if params.read_params.quals:
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
        
        if not params.skip_check_reads:
            reads_list = left_reads_list
            if right_reads_list:
                reads_list = reads_list + "," + right_reads_list
            params.read_params = check_reads(params, reads_list)
            
        user_supplied_juncs = []
        user_supplied_insertions = []
        user_supplied_deletions = []

        if params.gff_annotation and params.find_GFF_juncs:
            test_input_file(params.gff_annotation)
            (found_juncs, gtf_juncs) = get_gtf_juncs(params.gff_annotation)
            if found_juncs:
                user_supplied_juncs.append(gtf_juncs)
        if params.raw_junctions:
            test_input_file(params.raw_junctions)
            user_supplied_juncs.append(params.raw_junctions)

        if params.raw_insertions:
            test_input_file(params.raw_insertions)
            user_supplied_insertions.append(params.raw_insertions)

        if params.raw_deletions:
            test_input_file(params.raw_deletions)
            user_supplied_deletions.append(params.raw_deletions)

        # create a named pipe so we can have bowtie write unmapped reads into a compress pipe (when needed)
        global unmapped_reads_fifo
        unmapped_reads_fifo = tmp_dir + str(os.getpid())+".bwt_unmapped.z.fifo"
        if os.path.exists(unmapped_reads_fifo):
           os.unlink(unmapped_reads_fifo)
        try:
           os.mkfifo(unmapped_reads_fifo)
        except OSError, o:
           die(fail_str+"Error at mkfifo("+unmapped_reads_fifo+'). '+str(o))
                
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
        if params.read_params.integer_quals:
            params.read_params.integer_quals = False
            
        mapping = spliced_alignment(params, 
                                    bwt_idx_prefix,
                                    sam_header_filename,
                                    ref_fasta,
                                    params.read_params.seed_length,
                                    params.segment_length,
                                    left_kept_reads,
                                    right_kept_reads,
                                    user_supplied_juncs,
                                    user_supplied_insertions,
                                    user_supplied_deletions)
                                    
        os.unlink(unmapped_reads_fifo)
        
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

        if not params.system_params.keep_tmp:
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
