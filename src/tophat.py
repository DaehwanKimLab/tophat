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
import tempfile

from datetime import datetime, date, time

use_message = '''
TopHat maps short sequences from spliced transcripts to whole genomes.

Usage:
    tophat [options] <bowtie_index> <reads1.fq[,reads2.fq,...,readsN.fq]>
    
Options:
    -a/--min-anchor       <int>
    -p/--num-threads      <int>
    -i/--min-intron       <int>
    -I/--max-intron       <int>
    -X/--solexa-quals     <int>
    -s/--seed-length      <int>
    -g/--max-gene-family  <int>
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

output_dir = "./tophat_out/"
logging_dir = output_dir + "logs/"
bin_dir = sys.path[0] + "/"

#ok_str = "\t\t\t\t[OK]\n"
fail_str = "\t[FAILED]\n"

def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def prepare_output_dir():
    
    print >> sys.stderr, "[%s] Prepare output location %s" % (right_now(), output_dir)
    if os.path.exists(output_dir):
        pass
    else:        
        os.mkdir(output_dir)
        
    if os.path.exists(logging_dir):
        pass
    else:        
        os.mkdir(logging_dir)

def check_bowtie_index():
    print >> sys.stderr, "[%s] Checking for Bowtie index files" % right_now()
    # TODO
    #print >> sys.stderr, ok_str

def check_bfa():
    print >> sys.stderr, "[%s] Checking for binary fasta" % right_now()
    # TODO
    #print >> sys.stderr, ok_str

def check_index():
    check_bowtie_index()
    check_bfa()

def get_maq_version():
    
    # Launch Maq to capture its version info
    proc = subprocess.Popen('maq',stderr=subprocess.PIPE)
    stderr_value = proc.communicate()[1]
    maq_version = None
    maq_out = repr(stderr_value)
    
    # Find the version identifier
    version_str = "Version: "
    ver_str_idx = maq_out.find(version_str)
    if ver_str_idx != -1:
        nl = maq_out.find("\\n", ver_str_idx)
        version_val = maq_out[ver_str_idx + len(version_str):nl]
        maq_version = [int(x) for x in version_val.split('.')]
        
    return maq_version

def get_bowtie_version():

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

    return bowtie_version

def check_maq():
    print >> sys.stderr, "[%s] Checking for Maq" % right_now()
    maq_version = get_maq_version()
    if maq_version == None:
        print >> sys.stderr, "Error: Maq not found on this system"
        exit(1)
    elif maq_version[1] >= 7:
        return True
    else:
        return False

def check_bowtie():
    print >> sys.stderr, "[%s] Checking for Bowtie" % right_now()
    bowtie_version = get_bowtie_version()
    if bowtie_version == None:
        print >> sys.stderr, "Error: Bowtie not found on this system"
        exit(1)
    elif bowtie_version[1] < 9 or bowtie_version[2] < 8:
        print >> sys.stderr, "Error: TopHat requires Bowtie 0.9.9 or later"
        exit(1)

def formatTD(td):
  hours = td.seconds // 3600
  minutes = (td.seconds % 3600) // 60
  seconds = td.seconds % 60
  return '%02d:%02d:%02d' % (hours, minutes, seconds) 

def initial_mapping(bwt_idx_prefix, 
                    reads_list, 
                    output_dir, 
                    bowtie_threads, 
                    solexa_scale,
                    seed_length,
                    max_hits):
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Mapping reads with Bowtie" % start_time.strftime("%c"),
    
    # Setup Bowtie output redirects
    bwt_map = output_dir + "unspliced_map.bwtout"
    bwt_log = open(logging_dir + "unspliced_bwt.log", "w")
    
    unmapped_reads_fasta_name = output_dir + "unmapped.fa"
    unmapped_repeat_fasta_name = output_dir + "unmapped_repeat.fa"
    # Launch Bowtie
    try:    
        qual_format = ""
        if solexa_scale:
            qual_format = "--solexa-quals"
        
        bowtie_cmd = ["bowtie"]
        
        if solexa_scale:
            bowtie_cmd += [qual_format]
        
        bowtie_cmd += ["-p", str(bowtie_threads),
                      "--unfa", unmapped_reads_fasta_name,
                      "-l", str(seed_length),
                      "-k", str(max_hits),
                      "-m", str(max_hits + 1),
                      "--maxfa", unmapped_repeat_fasta_name,
                      bwt_idx_prefix, 
                      reads_list, 
                      bwt_map]   
        #print "\t executing: `%s'" % " ".join(bowtie_cmd)           
        subprocess.check_call(bowtie_cmd, 
                              stdout=open("/dev/null"), 
                              stderr=bwt_log)
    # Bowtie not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: Bowtie not found on this system.  Did you forget to include it in your PATH?"
    # Bowtie reported an error
    except subprocess.CalledProcessError:
        print >> sys.stderr, fail_str, "Error: could not execute Bowtie"
        exit(1)
    # Success    
    finish_time = datetime.now()
    duration = finish_time - start_time
    print >> sys.stderr, "\t\t\t[%s elapsed]" %  formatTD(duration)
    return (bwt_map, unmapped_reads_fasta_name)

def collect_unmapped_reads():
    print >> sys.stderr, "[%s] Collecting unmapped reads" % right_now()
    #print >> sys.stderr, ok_str

def convert_to_maq(use_long_maq_maps, bwt_map, bwt_idx_prefix):
    print >> sys.stderr, "[%s] Coverting alignments to Maq format" % right_now()
    
    convert_log = open(logging_dir + "convert_to_maq.log", "w")
    maq_map = output_dir + "unspliced_map.maqout"
    format_option = "-o"
    if use_long_maq_maps == True:
        format_option = ""
    convert_cmd = ["bowtie-maqconvert", 
                   format_option,
                   bwt_map, 
                   maq_map,
                   bwt_idx_prefix + ".bfa", 
                   bwt_map]            
        
    try:    
        retcode = subprocess.call(convert_cmd, stderr=convert_log)
        # bowtie-maqconvert reported an error
        if retcode > 0:
            print >> sys.stderr, fail_str, "Error: Conversion to Maq map format failed"
            exit(1)
    # converter not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie-maqconvert not found on this system"
            exit(1)
            
    # Bowtie reported an error
    except subprocess.CalledProcessError:
        print >> sys.stderr, fail_str, "Error: could not execute Bowtie"
        exit(1)
    
    # Success
    return maq_map

def assemble_islands(maq_map, idx_prefix):
    print >> sys.stderr, "[%s] Assembling coverage islands" % right_now()
    asm_log = open(logging_dir + "maq_asm.log", "w")
    maq_cns = output_dir + "unspliced.cns"
    asm_cmd = ["maq",
               "assemble",
               "-s",
               maq_cns,
               idx_prefix + ".bfa",
               maq_map]            
    try:    
       retcode = subprocess.call(asm_cmd, stderr=asm_log)
       
       # Maq assembler returned an error 
       if retcode > 0:
           print >> sys.stderr, fail_str, "Error: Maq assembly failed"
           exit(1)
    # Maq not found
    except OSError, o:
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: Maq not found on this system"
       exit(1)
    return maq_cns

def extract_islands(maq_cns):
    print >> sys.stderr, "[%s] Extracting coverage islands" % right_now()
    extract_log = open(logging_dir + "extract_islands.log", "w")
    island_fasta = output_dir + "islands.fa"
    island_gff = output_dir + "islands.gff"
    
    extract_cmd = [bin_dir + "cvg_islands",
                   "-d", "0.0", # Minimum average depth of coverage threshold
                   "-b", "6", # Max gap length
                   "-e", "45", # Extension length
                   "-R", # Always take reference sequence (no SNP calls)
                   maq_cns,
                   island_fasta,
                   island_gff]            
    try:    
       retcode = subprocess.call(extract_cmd, stderr=extract_log)
       
       # cvg_islands returned an error 
       if retcode > 0:
           print >> sys.stderr, fail_str, "Error: Islands extraction failed"
           exit(1)
    # cvg_islands not found
    except OSError, o:
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: cvg_islands not found on this system"
       exit(1)
    return (island_fasta, island_gff)

def align_spliced_reads(islands_fasta, 
                        islands_gff, 
                        unmapped_reads,
                        seed_length,
                        splice_mismatches):
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Aligning spliced reads" % start_time.strftime("%c"),
    splice_log = open(logging_dir + "spliced_align.log", "w")
    spliced_reads_name = output_dir + "spliced_map.sbwtout"
    spliced_reads = open(spliced_reads_name,"w")
    splice_cmd = [bin_dir + "spanning_reads",
                  "-v",
                  "-a", "5", # Anchor length
                  "-m", str(splice_mismatches), # Mismatches allowed in extension
                  "-I", "20000", # Maxmimum intron length
                  "-i", "70", # Minimum intron length
                  "-s", str(seed_length), # Seed size for reads
                  "-S", "300", # Min normalized DoC for self island junctions
                  "-M", "256", # Small memory footprint for now
                  islands_fasta,
                  islands_gff,
                  unmapped_reads]   
    #print "\n"," ".join(splice_cmd)         
    try:    
       retcode = subprocess.call(splice_cmd, 
                                 stderr=splice_log,
                                 stdout=spliced_reads)
       
       # spanning_reads returned an error 
       if retcode > 0:
           print >> sys.stderr, fail_str, "Error: Spliced read alignment failed"
           exit(1)
    # cvg_islands not found
    except OSError, o:
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: spanning_reads not found on this system"
       exit(1)
       
    finish_time = datetime.now()
    duration = finish_time - start_time
    print >> sys.stderr, "\t\t\t[%s elapsed]" %  formatTD(duration)
    return spliced_reads_name

def compile_reports(contiguous_map, spliced_map):
    print >> sys.stderr, "[%s] Reporting output tracks" % right_now()
    report_log = open(logging_dir + "reports.log", "w")
    junctions = output_dir + "junctions.bed"
    coverage = output_dir + "coverage.wig"
    report_cmd = [bin_dir + "tophat_reports",
                  coverage,
                  junctions,
                  contiguous_map,
                  spliced_map]            
    try:    
       retcode = subprocess.call(report_cmd, 
                                 stderr=report_log)
       
       # spanning_reads returned an error 
       if retcode > 0:
           print >> sys.stderr, fail_str, "Error: Report generation failed"
           exit(1)
    # cvg_islands not found
    except OSError, o:
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: tophat_reports not found on this system"
       exit(1)
    return (coverage, junctions)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:vXp:s:m:", 
                                        ["help", 
                                         "output=", 
                                         "solexa-quals",
                                         "num-threads=",
                                         "seed-length=",
                                         "splice-mismatches=",
                                         "max-gene-family="])
        except getopt.error, msg:
            raise Usage(msg)
        
        bowtie_threads = 1
        solexa_scale = False
        seed_length = 28
        splice_mismatches = 0
        max_hits = 10
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                output = value
            if option in ("-X", "--solexa-quals"):
                solexa_scale = True
            if option in ("-p", "--num-threads"):
                bowtie_threads = int(value)
            if option in ("-s", "--seed-length"):
                seed_length = int(value)
            if option in ("-m", "--splice-mismatches"):
                splice_mismatches = int(value)
            if option in ("-g", "--max-gene-family"):
                 max_hits = int(value)
                
        if len(args) < 2:
            raise Usage(use_message)
            
        bwt_idx_prefix = args[0]
        reads_list = args[1]
        
        print >> sys.stderr
        print >> sys.stderr, "[%s] Beginning TopHat run" % right_now()
        print >> sys.stderr, "-----------------------------------------------" 
        
        check_index()
        use_long_maq_maps = check_maq()
        check_bowtie()
        prepare_output_dir()
        (bwt_map, unmapped_reads) = initial_mapping(bwt_idx_prefix, 
                                                    reads_list, 
                                                    output_dir,
                                                    bowtie_threads,
                                                    solexa_scale,
                                                    seed_length,
                                                    max_hits)
        
        maq_map = convert_to_maq(use_long_maq_maps, bwt_map, bwt_idx_prefix)
        maq_cns = assemble_islands(maq_map, bwt_idx_prefix)
        (islands_fasta, islands_gff) = extract_islands(maq_cns)
        spliced_reads = align_spliced_reads(islands_fasta, 
                                            islands_gff, 
                                            unmapped_reads,
                                            seed_length,
                                            splice_mismatches)
        compile_reports(bwt_map, spliced_reads)
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for detailed help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
