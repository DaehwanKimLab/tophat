//
//  gtfToFasta.h
//  TopHat
//
//  Created by Harold Pimentel on 10/26/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef GTFToFasta_H
#define GTFToFasta_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include "common.h"
#include "gff.h"
#include "GFaSeqGet.h"
#include "FastaTools.h"

std::string get_exonic_sequence(GffObj *p_trans, FastaRecord *rec);

class GTFToFasta {
public:
    GTFToFasta(std::string gtf_fname, std::string genome_fname);
    ~GTFToFasta();
    void make_transcriptome(std::string out_fname);
  
    // for debugging
    void print_mapping();
private:
    GffReader gtfReader_;
    GFaSeqGet genome_fhandle_;
    
    std::string gtf_fname_;
    std::string genome_fname_;
    
    FILE* gtf_fhandle_;
    
    // "contig" => vector(index_of_gff_obj)
//    typedef std::map<const char* , std::vector< size_t >* > ContigTransMap;
    typedef std::map<std::string, std::vector< size_t >* > ContigTransMap;
    ContigTransMap contigTransMap_;
    
    void transcript_map();
    
    GTFToFasta(); // Don't want anyone calling the constructor w/o options
};

#endif
