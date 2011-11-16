/*
 * Author: Harold Pimentel
 * Contact: http://cs.berkeley.edu/~pimentel
 * Date: June 10, 2011
 */

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <bam/bam.h>
#include <bam/sam.h>

#include <seqan/sequence.h>

#include <getopt.h>
#include <unistd.h>

#include "bwt_map.h"
#include "common.h"
#include "gff.h"

#define MAX_READ_NAME_LEN 2048

/*
 * XXX: This class currently assumes someone used the script in TopHat to map
 *      the reads already. It also depends on that same format.
 */
class Map2GTF
{
public:
    Map2GTF(std::string gtf_fname, std::string sam_fname);
    ~Map2GTF();
    // Write out to a BAM file
    void convert_coords(std::string out_fname, std::string sam_header);

private:
    GffReader gtfReader_;

    std::string gtf_fname_;
    std::string reads_fname_;

    FILE* gtf_fhandle_;
    FILE* reads_fhandle_;
    FZPipe reads_pipe_;

    ReadTable readTable_;
    RefSequenceTable refSeqTable_;

    HitFactory* hitFactory_;
    HitStream* hitStream_;

    Map2GTF(); // Don't want anyone calling the constructor w/o options
};

class TranscriptomeHit {
  public:
    BowtieHit hit;
    GffObj* trans;
  TranscriptomeHit(GffObj* t=NULL):hit() {
    trans=t;
    }
  bool operator==(const TranscriptomeHit& th) const {
    return (th.hit == hit);
    }
  bool operator<(const TranscriptomeHit& th) const {
    return (th.hit < hit);
    }

};


void trans_to_genomic_coords(const char* read_name, HitFactory* hitFactory,
        const BowtieHit& in, TranscriptomeHit& out);

bool get_read_start(GList<GffExon>* exon_list, size_t gtf_start,
        size_t& genome_start, int& exon_idx);

void print_trans(GffObj* trans, const BowtieHit& in, size_t rem_len,
        size_t match_len, size_t cur_pos, size_t start_pos);
