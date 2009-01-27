#!/usr/bin/env python
# encoding: utf-8
"""
junction-sensitivity.py

Created by Cole Trapnell on 2009-01-08.
Copyright (c) 2009 Cole Trapnell. All rights reserved.
"""

import sys
import getopt


help_message = '''
Takes output from ASTD/transcript-junctions.py and a UCSC Bed track produced
by TopHat, reports sensitivity and specificity

Usage: <ground> <predicted> <# true negatives>
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


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
        ground_truth = open(args[0])
        predicted = open(args[1])
        
        #true_negatives = int(args[2])
        
        juncs = set([])
        
        for line in ground_truth:
            line = line.strip()
            cols = line.split()
            chrom = cols[0]
            chrom = chrom[3:]
            strand = cols[2]
            start = cols[3]
            stop = cols[4]
            
            juncs.add((chrom, strand, start, stop))
            
        num_predicted = 0
        num_correct = 0
        for line in predicted:
            if "track" in line:
                continue
            num_predicted += 1
            
            line = line.strip()
            cols = line.split()
            
            chrom = cols[0][3:]
            
            left_overhang = int(cols[10].split(',')[0]) 
            right_offset = int(cols[11].split(',')[1])
            start = str(int(cols[1]) + left_overhang)
            stop = str(int(cols[1]) + right_offset + 1)
            strand  = cols[5]
            #print chrom, strand, start, stop
            if (chrom, strand, start, stop) in juncs:
                num_correct += 1
        
        tp = num_correct
        fp = num_predicted - num_correct
        
        #fpr = (num_predicted - num_correct) / float(true_negatives)
        #tpr = (num_correct) / float(len(juncs))
        #print tpr, fpr
        print num_correct, num_predicted, len(juncs)
    
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
