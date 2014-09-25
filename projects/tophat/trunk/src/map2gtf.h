/*
 * Author: Harold Pimentel
 * Contact: http://cs.berkeley.edu/~pimentel
 * Date: June 10, 2011
 */

#ifndef _MAP2GTF_H_
#define _MAP2GTF_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <bam.h>
#include <sam.h>

#include <seqan/sequence.h>

#include <getopt.h>
#include <unistd.h>

#include "bwt_map.h"
#include "common.h"

//#include "gff.h"
#include "GVec.hh"

#define MAX_READ_NAME_LEN 2048

//simplified version of GffObj
//parsed from a simple text line
struct GffTranscript: public GSeg {
	GVec<GSeg> exons;
	int numID; //numeric ID in tlst
	std::string gffID;
	std::string refID;
	char strand;
	GffTranscript():exons(1), numID(-1), gffID(),
			refID(), strand(0) { }

	string& getRefName() {
		return refID;
	}
	GffTranscript(const std::string& tline);
};



/*
 * XXX: This class currently assumes someone used the script in TopHat to map
 *      the reads already. It also depends on that same format.
 */
class TranscriptomeHit;
class Map2GTF
{
public:
    Map2GTF(const std::string& gtf_fname, const std::string& sam_fname);
    ~Map2GTF();
    // Write out to a BAM file
    bool next_read_hits(std::vector<bam1_t*>& hits, size_t& num_hits, long& read_id);
    void convert_coords(const std::string& out_fname, const std::string& sam_header);
    bool trans_to_genomic_coords(TranscriptomeHit& hit);

private:
    //GffReader gtfReader_;
    GPVec<GffTranscript> transcripts;
    map<int, GffTranscript*> tidx_to_t;
    std::string gtf_fname_;
    std::string in_fname_;

    //FILE* gtf_fhandle_; //actually a tlst handle
    std::ifstream tlststream;
    samfile_t* in_fhandle_;
    bam_header_t* in_sam_header_;

    map<string, int> ref_to_id_;
    bam_header_t* out_sam_header_;

    ReadTable readTable_;
    RefSequenceTable refSeqTable_;

    Map2GTF(); // Don't want anyone calling the constructor w/o options
};

class TranscriptomeHit
{
 public:
  bam1_t* hit;
  GffTranscript* trans;
  TranscriptomeHit(bam1_t* h = NULL, GffTranscript* t=NULL): //GffObj* t = NULL)
  	hit(h), trans(t)
    { }
  bool operator==(const TranscriptomeHit& th) const
  {
    if (hit->core.tid != th.hit->core.tid)
      return false;
    
    if (hit->core.pos != th.hit->core.pos)
      return false;
    
    if (hit->core.n_cigar != th.hit->core.n_cigar)
      return false;
    
    for (int i = 0; i < hit->core.n_cigar; ++i)
      {
	if (bam1_cigar(hit)[i] != bam1_cigar(th.hit)[i])
	  return false;
      }
    
    return true;
  }
  bool operator<(const TranscriptomeHit& th) const
  {
    if (hit->core.tid != th.hit->core.tid)
      return hit->core.tid < th.hit->core.tid;
    
    if (hit->core.pos != th.hit->core.pos)
      return hit->core.pos < th.hit->core.pos;
    
    if (hit->core.n_cigar != th.hit->core.n_cigar)
      return hit->core.n_cigar < th.hit->core.n_cigar;
    
    for (int i = 0; i < hit->core.n_cigar; ++i)
      {
	if (bam1_cigar(hit)[i] != bam1_cigar(th.hit)[i])
	  return bam1_cigar(hit)[i] < bam1_cigar(th.hit)[i];
      }
    
    return false;
  }  
};

//bool get_read_start(GList<GffExon>* exon_list, size_t gtf_start,
bool get_read_start(GVec<GSeg>& exon_list, size_t gtf_start,
		 size_t& genome_start, int& exon_idx);

void print_trans(GffTranscript* trans, const bam1_t* in, size_t rem_len,
		 size_t match_len, size_t cur_pos, size_t start_pos);

#endif /* _MAP2GTF_H_ */
