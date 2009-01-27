#!/usr/bin/env python
# encoding: utf-8
"""
count-false-junctions.py

Created by Cole Trapnell on 2008-09-09.
Copyright (c) 2008 Cole Trapnell. All rights reserved.
"""

import sys
import getopt


help_message = '''
Given a junctions file and a list of transcript ids, this script counts
how many exons in those transcripts are close enough that they could be 
spliced together, and yet do not appear in any of the specified transcripts.

    count-false-junctions.py [options] <transcript.ids> <ASTD_junctions>
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:vG:", ["help", 
                                                            "output=",
                                                            "max-intron-length"])
        except getopt.error, msg:
            raise Usage(msg)
    
        max_intron_length = 20000
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                output = value
            if option in ("-G", "--max-intron-length"):
                max_intron_length = int(value)
        
        if len(args) < 2:
            raise Usage(help_message)
        
        ids_file = open(args[0])
        ids = set([line.strip() for line in ids_file])
        
        juncs_file = open(args[1])
        
        exons_by_chr = {}
        true_junctions_by_chr = {}
        
        num_true_juncs = 0
        
        for line in juncs_file:
            
            line = line.strip()
            cols = line.split()
            if not cols[1] in ids:
                continue
                
            chromosome = cols[0]
            orient = cols[2]
            exons = exons_by_chr.setdefault((chromosome,orient), set([]))
            #true_juncs = true_junctions_by_chr.setdefault((chromosome, orient), set([]))
            #true_juncs.add((int(cols[3]), int(cols[4]), cols[5], cols[6]))
            num_true_juncs += 1
            exons.add((cols[5], int(cols[3])))
            exons.add((cols[6], int(cols[4])))
        
        
        
        accepted_pairs = set([])
        
        for ((chromosome,orient), exons) in exons_by_chr.iteritems():
           
            for exon1 in exons:
                #print "ex1:", exon1
                for exon2 in exons:
                    #print "\tex2:", exon2
                    dist = (exon2[1] - exon1[1])
                    if exon1[0] != exon2[0] and dist < max_intron_length and dist > 0:
                        #print "\t\taccepted pair:", (exon1,exon2)
                        accepted_pairs.add((exon1[0], exon2[0]))
        
        num_false_negatives = len(accepted_pairs) - num_true_juncs
        
        print >> sys.stderr, len(accepted_pairs), "Possible junctions,", num_true_juncs, "True junctions"
        print num_false_negatives
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
