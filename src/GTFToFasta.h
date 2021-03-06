//
//  gtfToFasta.h
//  TopHat
//
//  Created by Harold Pimentel on 10/26/11.
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

std::string get_exonic_sequence(GffObj& p_trans, FastaRecord& rec, std::string& coords);

class GTFToFasta {
public:
    GTFToFasta(std::string gtf_fname, std::string genome_fname);
    ~GTFToFasta();
    void make_transcriptome(std::string out_fname);

    // for debugging
    void print_mapping();
private:
    GffReader gtfReader_;
    // The genome_fhandle_ isn't used anywhere after being
    // initialized, and double opening the fasta file means that pipes
    // cannot work.
    // GFaSeqGet genome_fhandle_;

    std::string gtf_fname_;
    std::string genome_fname_;

    FILE* gtf_fhandle_;

    // "contig" => vector(index_of_gff_obj)
//    typedef std::map<const char* , std::vector< size_t >* > ContigTransMap;
    typedef std::map<std::string, std::vector< int >* > ContigTransMap;

    ContigTransMap contigTransMap_;

    void transcript_map();

    GTFToFasta(); // Don't want anyone calling the constructor w/o options
};

#endif
