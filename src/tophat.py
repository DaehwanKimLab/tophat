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
import warnings

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
    
    print >> sys.stderr, "[%s] Preparing output location %s" % (right_now(), output_dir)
    if os.path.exists(output_dir):
        pass
    else:        
        os.mkdir(output_dir)
        
    if os.path.exists(logging_dir):
        pass
    else:        
        os.mkdir(logging_dir)

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
            exit(1)
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
            exit(1)

def bowtie_idx_to_bfa(idx_prefix):
    
    idx_name = idx_prefix.split('/')[-1]
    idx_bfa = output_dir + idx_name + ".bfa"
    print >> sys.stderr, "[%s] Building Maq binary fasta file %s \n\tfrom Bowtie index" % (right_now(), idx_bfa)
    
    try:    
        
        tmp_fasta_file_name = output_dir + idx_name + ".fa"
        tmp_fasta_file = open(tmp_fasta_file_name, "w")

        inspect_log = open(logging_dir + "bowtie_inspect.log", "w")

        inspect_cmd = ["bowtie-inspect",
                      idx_prefix]
        #print >> sys.stderr, "Executing: " + " ".join(inspect_cmd) + " > " + tmp_fasta_file_name   
        subprocess.check_call(inspect_cmd, 
                              stdout=tmp_fasta_file,
                              stderr=inspect_log)
    # Bowtie not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie-inspect not found on this system.  Did you forget to include it in your PATH?"
    # Bowtie reported an error
    except subprocess.CalledProcessError:
        print >> sys.stderr, fail_str, "Error: bowtie-inspect returned an error"
        exit(1) 

    try:    
        
        fasta2bfa_log = open(logging_dir + "fasta2bfa.log", "w")

        fasta2bfa_cmd = ["maq",
                         "fasta2bfa",
                         tmp_fasta_file_name,
                         idx_bfa]
                          
        subprocess.check_call(fasta2bfa_cmd, 
                              stderr=fasta2bfa_log)
    # Maq not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: Maq not found on this system.  Did you forget to include it in your PATH?"
            exit(1)
    # Maq reported an error
    except subprocess.CalledProcessError:
        print >> sys.stderr, fail_str, "Error: maq fasta2bfa returned an error"
        exit(1)   
    os.remove(tmp_fasta_file_name)
        
    return idx_bfa

def check_bfa(idx_prefix):
    print >> sys.stderr, "[%s] Checking for binary fasta" % right_now()
    idx_bfa = idx_prefix + ".bfa"
    if os.path.exists(idx_bfa):
        return idx_bfa
    else:
        idx_name = idx_prefix.split('/')[-1]
        bowtie_idx_env_var = os.environ.get("BOWTIE_INDEXES")
        if bowtie_idx_env_var != None:
            idx_bfa = bowtie_idx_env_var + idx_prefix + ".bfa" 
            if os.path.exists(idx_bfa):
                return idx_bfa
        
        print >> sys.stderr, "Warning: Could not find Maq binary fasta file " + idx_bfa
        idx_bfa = bowtie_idx_to_bfa(idx_prefix)
        return idx_bfa
        #print >> sys.stderr, "Error: Could not find Maq binary fasta file " + idx_bfa
        #exit(1)
    
def check_index(idx_prefix):
    check_bowtie_index(idx_prefix)
    return check_bfa(idx_prefix)

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
    if len(bowtie_version) == 3:
        bowtie_version.append(0)
        
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
    elif bowtie_version[1] < 9 or bowtie_version[2] < 8 or (bowtie_version == [0,9,8,0]):
        print >> sys.stderr, "Error: TopHat requires Bowtie 0.9.8.1 or later"
        exit(1)
    print >> sys.stderr, "\tBowtie version:\t\t %s" % ".".join([str(x) for x in bowtie_version])
        
def check_reads(reads_files, default_seed_len, default_format, solexa_scale):
    print >> sys.stderr, "[%s] Checking reads" % right_now()
    bowtie_version = get_bowtie_version()
    
    seed_len = default_seed_len
    format = default_format
    
    observed_formats = set([])
    observed_scales = set([])
    min_seed_len = 99999
    max_qual = -1
    files = reads_files.split(',')
    for f_name in files:
        try:
            f = open(f_name)
        except IOError:
            print >> sys.stderr, "Error: could not open file", f_name
            exit(1)
            
        first_line = f.readline()
        if first_line[0] == "@":
            format = "fastq"
        elif first_line[0] == ">":
            format = "fasta"
        else:
            print >> sys.stderr, "Error: file %s does not appear to be a valid FASTA or FASTQ file" % f_name
        observed_formats.add(format)
        f.seek(0)
        line_num = 0
        if format == "fastq":
            for line in f:
                if line_num % 4 == 1:
                    seq_len = len(line) - 1
                    if seq_len < 20:
                        print >> sys.stderr, "Warning: found a read < 20bp in", f_name
                    else:
                        min_seed_len = min(seq_len, min_seed_len)                
                elif line_num % 4 == 3:
                    max_line_qual = max([ord(x) for x in list(line.strip())])
                    max_qual = max(max_line_qual, max_qual)
                line_num += 1
            if max_qual > 90:
                solexa_scale = True

        elif format == "fasta":
            for line in f:
                if line_num % 2 == 1:
                    seq_len = len(line) - 1
                    if seq_len < 20:
                        print >> sys.stderr, "Warning: found a read < 20bp in", f_name
                    else:
                        min_seed_len = min(seq_len, min_seed_len)
                line_num += 1
            
    if len(observed_formats) > 1:
        print >> sys.stderr, "Error: TopHat requires all reads be either FASTQ or FASTA.  Mixing formats is not supported."
        exit(1)
          
    seed_len = min(seed_len, min_seed_len)
    print >> sys.stderr, "\tseed length:\t %dbp" % seed_len
    print >> sys.stderr, "\tformat:\t\t %s" % format
    if format == "fastq":
        print >> sys.stderr, "\tquality scale:\t %s" % (solexa_scale and "solexa" or "phred")
    
    #print seed_len, format, solexa_scale
    return seed_len, format, solexa_scale
    
    

def formatTD(td):
  hours = td.seconds // 3600
  minutes = (td.seconds % 3600) // 60
  seconds = td.seconds % 60
  return '%02d:%02d:%02d' % (hours, minutes, seconds) 

def filter_garbage(reads_list, reads_format):
    
    try:    
        #filter_cmd = ["filter_garbage"]
        if reads_format == "-f":
            reads_suffix = ".fa"
        else:
            reads_suffix = ".fq"
            
        kept_reads_filename = output_dir + "kept_reads" + reads_suffix
        kept_reads = open(kept_reads_filename, "w")
        
        filter_log = open(logging_dir + "filter_garbage.log", "w")
        
        filter_cmd = ["filter_garbage",
                      reads_format]   
        #print "\t executing: `%s'" % " ".join(bowtie_cmd)    
        files = reads_list.split(',')
        for reads_file in files:       
            subprocess.check_call(filter_cmd, 
                                  stdin=open(reads_file),
                                  stdout=kept_reads,
                                  stderr=filter_log)
    # Bowtie not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: Bowtie not found on this system.  Did you forget to include it in your PATH?"
    # Bowtie reported an error
    except subprocess.CalledProcessError:
        print >> sys.stderr, fail_str, "Error: could not execute Bowtie"
        exit(1)    
    
    return kept_reads_filename

def initial_mapping(bwt_idx_prefix, 
                    reads_list,
                    reads_format, 
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
        
        bowtie_cmd += [reads_format,
                       "-p", str(bowtie_threads),
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

def convert_chunk_to_maq(use_long_maq_maps, bwt_map, maq_map, idx_bfa, convert_log):
    #convert_log = open(logging_dir + "convert_to_maq.log", "w")
    format_option = "-o"
    if use_long_maq_maps == True:
        format_option = ""
    convert_cmd = ["bowtie-maqconvert", 
                   format_option,
                   bwt_map, 
                   maq_map,
                   idx_bfa, 
                   bwt_map]            
    print >> convert_log, "Converting %s to %s" % (bwt_map, maq_map)
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
        print >> sys.stderr, fail_str, "Error: could not execute bowtie-maqconvert"
        exit(1)

    # Success
    return maq_map

def convert_to_maq(use_long_maq_maps, bwt_map, idx_bfa, alignments_per_chunk=10000000):
    print >> sys.stderr, "[%s] Converting alignments to Maq format" % right_now()
    
    maq_map = output_dir + "unspliced_map.maqout"
    big_bwt_map = open(bwt_map)
    i = 0
    
    convert_log = open(logging_dir + "convert_to_maq.log", "w")
    #convert_log = open("/dev/null", "w")
    tmp_bwt_name = os.tmpnam()
    tmp_bwt = open(tmp_bwt_name,"w")
    #print tmp_bwt_name
    tmp_maq = os.tmpnam()
    #print tmp_maq
    tmp_maps = [tmp_maq]
    tmp_bwts = [tmp_bwt_name]
    num_chunks = 1
    for line in big_bwt_map:
        i += 1
        if i >= alignments_per_chunk: 
            #print "converting chunk", num_chunks
            tmp_bwt.flush()
            convert_chunk_to_maq(use_long_maq_maps, tmp_bwt_name, tmp_maq, idx_bfa, convert_log)
            tmp_bwt_name = os.tmpnam()
            tmp_bwts.append(tmp_bwt_name)
            tmp_bwt = open(tmp_bwt_name,"w")
            tmp_maq = os.tmpnam()
            tmp_maps.append(tmp_maq)
            num_chunks += 1
            i = 0
        print >> tmp_bwt, line,
    if i > 0:
        #print "converting chunk", num_chunks
        tmp_bwt.flush()
        convert_chunk_to_maq(use_long_maq_maps, tmp_bwt_name, tmp_maq, idx_bfa, convert_log)         

    convert_cmd = ["maq",
                   "mapmerge",
                   maq_map]
                   
    convert_cmd.extend(tmp_maps)            
    print >> convert_log, " ".join(convert_cmd)
    try:    
        retcode = subprocess.call(convert_cmd, stderr=convert_log)
        # bowtie-maqconvert reported an error
        if retcode > 0:
            print >> sys.stderr, fail_str, "Error: Conversion to Maq map format failed"
            exit(1)
    # converter not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: maq not found on this system"
            exit(1)

    # Bowtie reported an error
    except subprocess.CalledProcessError:
        print >> sys.stderr, fail_str, "Error: could not execute maq mapmerge"
        exit(1)
    
    #print " ".join(tmp_maps)
    for m in tmp_maps:
        #print m
        os.remove(m)
    for m in tmp_bwts:
        #print m
        os.remove(m)
        
    # Success
    return maq_map

def assemble_islands(maq_map, idx_bfa):
    print >> sys.stderr, "[%s] Assembling coverage islands" % right_now()
    
    
    asm_log = open(logging_dir + "maq_asm.log", "w")
    maq_cns = output_dir + "unspliced.cns"
    asm_cmd = ["maq",
               "assemble",
               "-s",
               maq_cns,
               idx_bfa,
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

def extract_islands(maq_cns, island_gap, extend_islands):
    print >> sys.stderr, "[%s] Extracting coverage islands" % right_now()
    
    
    extract_log = open(logging_dir + "extract_islands.log", "w")
    island_fasta = output_dir + "islands.fa"
    island_gff = output_dir + "islands.gff"
    
    extract_cmd = [bin_dir + "cvg_islands",
                   "-d", "0.0", # Minimum average depth of coverage threshold
                   "-b", str(island_gap), # Max gap length
                   "-e", str(extend_islands), # Extension length
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
                        min_anchor_len,
                        splice_mismatches,
                        min_intron_length,
                        max_intron_length,
                        max_mem):
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Aligning spliced reads" % start_time.strftime("%c"),
    
    
    splice_log = open(logging_dir + "spliced_align.log", "w")
    spliced_reads_name = output_dir + "spliced_map.sbwtout"
    spliced_reads = open(spliced_reads_name,"w")
    splice_cmd = [bin_dir + "spanning_reads",
                  "-v",
                  "-a", str(min_anchor_len), # Anchor length
                  "-m", str(splice_mismatches), # Mismatches allowed in extension
                  "-I", str(max_intron_length), # Maxmimum intron length
                  "-i", str(min_intron_length), # Minimum intron length
                  "-s", str(seed_length), # Seed size for reads
                  "-S", "300", # Min normalized DoC for self island junctions
                  "-M", str(max_mem), # Small memory footprint for now
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

def compile_reports(contiguous_map, spliced_map, min_isoform_fraction):
    print >> sys.stderr, "[%s] Reporting output tracks" % right_now()
    
    
    report_log = open(logging_dir + "reports.log", "w")
    junctions = output_dir + "junctions.bed"
    coverage = output_dir + "coverage.wig"
    report_cmd = [bin_dir + "tophat_reports",
                  "-F", str(min_isoform_fraction),
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
            opts, args = getopt.getopt(argv[1:], "hvXp:s:m:M:F:a:i:I:e:b:", 
                                        ["help",  
                                         "solexa-quals",
                                         "num-threads=",
                                         "seed-length=",
                                         "splice-mismatches=",
                                         "max-gene-family=",
                                         "max-mem=",
                                         "min-isoform-fraction=",
                                         "min-anchor-length=",
                                         "min-intron-length=",
                                         "max-intron-length=",
                                         "island-gap="
                                         "extend-islands="])
        except getopt.error, msg:
            raise Usage(msg)
        
        bowtie_threads = 1
        min_anchor_len = 5
        solexa_scale = False
        seed_length = 28
        splice_mismatches = 0
        max_hits = 10
        max_mem = 1024
        reads_format = "fastq"
        min_isoform_fraction = 0.15
        min_intron_length = 70
        max_intron_length = 20000
        
        island_gap = 6
        extend_islands = 45
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-a", "--min-anchor"):
                min_anchor_len = int(value)
                if min_anchor_len < 3 or min_anchor_len > 6:
                    print >> sys.stderr, "Error: arg to --min-anchor-len must be in [3,6]"
                    exit(1)
            if option in ("-X", "--solexa-quals"):
                solexa_scale = True
            if option in ("-p", "--num-threads"):
                bowtie_threads = int(value)
            if option in ("-s", "--seed-length"):
                seed_length = int(value)
                if seed_length < 20:
                    print >> sys.stderr, "Error: arg to --seed-length must be at least 20"
                    exit(1)
            if option in ("-m", "--splice-mismatches"):
                splice_mismatches = int(value)
                if not splice_mismatches in [0,1,2]:
                    print >> sys.stderr, "Error: arg to --splice-mismatches must be 0, 1, or 2"
                    exit(1)
            if option in ("-g", "--max-gene-family"):
                 max_hits = int(value)
            if option in ("-M", "--max-mem"):
                max_mem = int(value)
            if option in ("-F", "--min-isoform-fraction"):
                min_isoform_fraction = float(value)
                if min_isoform_fraction < 0.0 or min_isoform_fraction > 1.0:
                    print >> sys.stderr, "Error: arg to --min-isoform-fraction must be between 0.0 and 1.0"
                    exit(1)
            if option in ("-i", "--min-intron-length"):
                min_intron_length = int(value)
                if min_intron_length <= 0:
                    print >> sys.stderr, "Error: arg to --min-intron-length must be greater than 0"
                    exit(1)                    
            if option in ("-I", "--max-intron-length"):
                max_intron_length = int(value)
                if max_intron_length <= 0:
                    print >> sys.stderr, "Error: arg to --max-intron-length must be greater than 0"
                    exit(1)
            if option in ("-e", "--extend-islands"):
                island_extension = int(value)
                if island_extension < 0:
                    print >> sys.stderr, "Error: arg to --extend-islands must be at least 0"
                    exit(1)
            if option in ("-b", "--island-gap"):
                island_gap = int(value)
                if island_gap < 0:
                    print >> sys.stderr, "Error: arg to --island-gap must be at least 0"
                    exit(1)
        if len(args) < 2:
            raise Usage(use_message)
            
        bwt_idx_prefix = args[0]
        reads_list = args[1]
        
        print >> sys.stderr
        print >> sys.stderr, "[%s] Beginning TopHat run" % right_now()
        print >> sys.stderr, "-----------------------------------------------" 
        
        start_time = datetime.now()
        
        idx_bfa = check_index(bwt_idx_prefix)
        use_long_maq_maps = check_maq()
        check_bowtie()
        
        (seed_length, reads_format, solexa_scale) = check_reads(reads_list, 
                                                                seed_length, 
                                                                reads_format, 
                                                                solexa_scale) 
        if reads_format == "fastq":
            format_flag = "-q"
        elif reads_format == "fasta":
            format_flag = "-f"
        
        prepare_output_dir()
        
        kept_reads = filter_garbage(reads_list, format_flag)
        kept_reads = reads_list
        (bwt_map, unmapped_reads) = initial_mapping(bwt_idx_prefix, 
                                                    kept_reads,
                                                    format_flag, 
                                                    output_dir,
                                                    bowtie_threads,
                                                    solexa_scale,
                                                    seed_length,
                                                    max_hits)
        warnings.filterwarnings("ignore", "tmpnam is a potential security risk")
        maq_map = convert_to_maq(use_long_maq_maps, bwt_map, idx_bfa)
        maq_cns = assemble_islands(maq_map, idx_bfa)
        (islands_fasta, islands_gff) = extract_islands(maq_cns,
                                                       island_gap,
                                                       extend_islands)
                                                       
        spliced_reads = align_spliced_reads(islands_fasta, 
                                            islands_gff, 
                                            unmapped_reads,
                                            seed_length,
                                            min_anchor_len,
                                            splice_mismatches,
                                            min_intron_length,
                                            max_intron_length,
                                            max_mem)
        compile_reports(bwt_map, spliced_reads, min_isoform_fraction)
        
        finish_time = datetime.now()
        duration = finish_time - start_time
        print >> sys.stderr,"-----------------------------------------------"
        print >> sys.stderr, "Run complete [%s elapsed]" %  formatTD(duration)
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for detailed help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
