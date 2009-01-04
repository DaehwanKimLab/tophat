#!/usr/bin/env python
# encoding: utf-8
"""
tophat.py

Created by Cole Trapnell on 2008-12-25.
Copyright (c) 2008 Cole Trapnell. All rights reserved.
"""

import sys
import getopt
import subprocess
import errno
import os
from datetime import datetime, date, time

use_message = '''
TopHat maps short sequences from spliced transcripts to whole genomes.

Usage:
    tophat [options] <bowtie_index> <reads1.fq[,reads2.fq,...,readsN.fq]>
    
Options:
    -a/--min-anchor     <int>
    -p/--num-threads    <int>
    -i/--min-intron     <int>
    -I/--max-intron     <int>
    -X/--use-solexa     <int>
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

bowtie_threads = 1
output_dir = "./tophat_out/"

#ok_str = "\t\t\t\t[OK]\n"
fail_str = "\t[FAILED]\n"

def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def initial_mapping(bwt_idx_prefix, reads_list, output_dir):
    
    print >> sys.stderr, "[%s] Mapping reads with Bowtie" % right_now()
    
    # Setup Bowtie output redirects
    bwt_map = output_dir + "unspliced_map.bwtout"
    bwt_log = open("/dev/null", "w")
    
    # Launch Bowtie
    try:    
        bowtie_cmd = ["bowtie", 
                      "-p", 
                      str(bowtie_threads), 
                      bwt_idx_prefix, 
                      reads_list, 
                      bwt_map]              
        
        subprocess.check_call(bowtie_cmd, stderr=bwt_log)
    
    # Bowtie not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: Bowtie not found on this system.  Did you forget to include it in your PATH?"
    
    # Bowtie reported an error
    except subprocess.CalledProcessError:
        print >> sys.stderr, fail_str, "Error: could not execute Bowtie"
        exit(1)
    
    # Success    
    #print >> sys.stderr, ok_str
  
def prepare_output_dir():
    
    print >> sys.stderr, "[%s] Prepare output location %s" % (right_now(), output_dir)
    if os.path.exists(output_dir):
        pass
    else:        
        os.mkdir(output_dir)

def check_bowtie_index():
    print >> sys.stderr, "[%s] Checking for Bowtie index files" % right_now()
    #print >> sys.stderr, ok_str

def check_bfa():
    print >> sys.stderr, "[%s] Checking for binary fasta" % right_now()
    #print >> sys.stderr, ok_str
    
def check_index():
    check_bowtie_index()
    check_bfa()
    
def get_maq_version():
    pass

def check_maq():
    print >> sys.stderr, "[%s] Checking for Maq" % right_now()
    #print >> sys.stderr, ok_str
    
def check_bowtie():
    print >> sys.stderr, "[%s] Checking for Bowtie" % right_now()
    #print >> sys.stderr, ok_str

def collect_unmapped_reads():
    print >> sys.stderr, "[%s] Collecting unmapped reads" % right_now()
    #print >> sys.stderr, ok_str

def convert_to_maq():
    print >> sys.stderr, "[%s] Coverting alignments to Maq format" % right_now()
    #print >> sys.stderr, ok_str
    
def assemble_islands():
    print >> sys.stderr, "[%s] Assembling coverage islands" % right_now()
    #print >> sys.stderr, ok_str
    
def extract_islands():
    print >> sys.stderr, "[%s] Extracting coverage islands" % right_now()
    #print >> sys.stderr, ok_str
   
def align_spliced_reads():
    print >> sys.stderr, "[%s] Aligning spliced reads" % right_now()
    #print >> sys.stderr, ok_str
    
def compile_reports():
    print >> sys.stderr, "[%s] Reporting junctions" % right_now()
    #print >> sys.stderr, ok_str    
    
def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:v", ["help", "output="])
        except getopt.error, msg:
            raise Usage(msg)
        
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                output = value
                
        if len(args) < 2:
            raise Usage(use_message)
            
        bwt_idx_prefix = args[0]
        reads_list = args[1]
        
        print >> sys.stderr
        print >> sys.stderr, "[%s] Beginning TopHat run" % right_now()
        print >> sys.stderr, "-----------------------------------------------" 
        
        check_index()
        check_maq()
        check_bowtie()
        prepare_output_dir()
        initial_mapping(bwt_idx_prefix, reads_list, output_dir)
        collect_unmapped_reads()
        convert_to_maq()
        assemble_islands()
        extract_islands()
        align_spliced_reads()
        compile_reports()
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for detailed help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
