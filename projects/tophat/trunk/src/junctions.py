#!/usr/bin/env python2.5
# encoding: utf-8
"""
junctions.py

Created by Cole Trapnell on 2008-10-16.
Copyright (c) 2008 Cole Trapnell. All rights reserved.
"""

import sys
import getopt


help_message = '''
Takes a .splices file and compiles the list of junctions as a UCSC genome 
browser BED track.

Usage:
    junctions.py <in.splices>
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
            
        if len(args) < 1:
            raise Usage(help_message)
        
        splice_file = open(args[0])
        
        junctions = set([])
        
        
        print 'track name=junctions description="TopHat junctions"'
        for line in splice_file:
            line = line.strip()
            cols = line.split()
            
            if len(cols) < 9:
                continue
                
            read_id = cols[0]
            contig = cols[1]
            strand = cols[2]
            left_island = cols[3]
            
            # This is the position of the entire island within the contig
            left_island_pos = int(cols[4])
            
            # This is the position of the junction within the left island
            left_junction_pos = int(cols[5])
            
            right_island = cols[6]
            
            # This is the position of the entire island within the contig
            right_island_pos = int(cols[7])
            
            # This is the position of the junction within the left island
            right_junction_pos = int(cols[8])
            
            junctions.add((left_island_pos + left_junction_pos, 
                           right_island_pos + right_junction_pos,
                           contig, 
                           strand))
        
        junc_id = 1
        
        print >> sys.stderr, len(junctions), "Total junctions"
        
        for junc in junctions:
            contig = junc[2]
            strand = junc[3]
            contig_start = junc[0]
            contig_end = junc[1]
            
            print "%s\t%d\t%d\tJUNC%08d\t0\t%s\t%d\t%d\t255,0,0\t2\t1,1\t0,%d" % ( contig,
                contig_start, contig_end, junc_id, strand, contig_start, contig_end, contig_end - contig_start - 1)
            junc_id += 1
            
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
