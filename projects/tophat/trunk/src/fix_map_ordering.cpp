/*
 *  fix_map_ordering.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 2/28/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */


#include <queue>
#include <cstring>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <getopt.h>
#include "common.h"
#include "reads.h"
#include "bwt_map.h"

using namespace seqan;
using namespace std;

struct TabSplitLine {
  //split a text line into an array of strings t[]
  //holds a copy of the text line
  int tcount;
  char **t;
  char *str; //text line, with \0s instead of tabs
  int tcap;
  TabSplitLine(const char* line) {
    tcount=0;
    t=NULL;
    tcap=0;
    str=NULL;
    if (line==NULL) return;
    str=strdup(line);
    //Notes: * destructive operation for s (replaces every \t with \0)
    //       * user must free t when no longer needed
    tcap=14;
    t=(char**)malloc(tcap*sizeof(char*));
    char prevch=0;
    for (char* p = str; *p!=0 ;p++) {
      if (*p=='\t') *p=0; //break the string here
      if (prevch==0) { //field start
        if (tcount==tcap) {
           tcap+=4;
           t = (char**)realloc(t,tcap*sizeof(char*));
           }
        t[tcount]=p;
        tcount++;
        } //field start
      prevch=*p;
      if (*p=='\n' || *p=='\r') {
          *p=0;
          break;
          }
      }//for each character on the line
    }
  ~TabSplitLine() {
    if (str!=NULL) {
        free(str);
        free(t);
        }
    }
};

struct MapOrdering {
	bool operator()(pair<uint64_t, TabSplitLine*>& lhs, pair<uint64_t, TabSplitLine*>& rhs)
	{
		uint64_t lhs_id = lhs.first;
		uint64_t rhs_id = rhs.first;
		return lhs_id > rhs_id;
	}
};

// "AS:i" (alignment score) is considered.
struct BamMapOrdering {
  bool operator()(pair<uint64_t, bam1_t*>& lhs, pair<uint64_t, bam1_t*>& rhs) {
    uint64_t lhs_id = lhs.first;
    uint64_t rhs_id = rhs.first;

    if (lhs_id != rhs_id || !bowtie2)
      return lhs_id > rhs_id;

    int lhs_score, rhs_score;
    lhs_score = rhs_score = numeric_limits<int>::min();
    uint8_t* ptr = bam_aux_get(lhs.second, "AS");
    if (ptr)
      lhs_score = bam_aux2i(ptr);

    ptr = bam_aux_get(rhs.second, "AS");
    if (ptr)
      rhs_score = bam_aux2i(ptr);

    if (lhs_score != rhs_score)
      return lhs_score < rhs_score;

    return lhs.second->core.flag & BAM_FSECONDARY;
  }
};

void writeSamLine(TabSplitLine& l, FILE* f) {
  if (l.tcount<10) {
     //fprintf(stderr, "Warning: skipping malformed SAM line %s\n",samline);
     return;
     }
  int flag=atoi(l.t[1]); //FLAG
  if ((flag & BAM_FUNMAP) != 0)
		  return;
  fprintf(f, "%s", l.t[0]);
  for (int i=1;i<l.tcount;i++)
	  fprintf(f,"\t%s",l.t[i]);
  fprintf(f, "\n");
}

void writeSam2Bam(TabSplitLine& l, GBamWriter& wbam) {
  //let SAM lines through, except for unmapped reads
  if (l.tcount<10) {
     //fprintf(stderr, "Warning: skipping malformed SAM line %s\n",samline);
     return;
     }
  int flag=atoi(l.t[1]); //FLAG
  if ((flag & BAM_FUNMAP) != 0)
	  return;
  int gpos=isdigit(l.t[3][0]) ? atoi(l.t[3]) : 0;
  int insert_size=atoi(l.t[8]); //TLEN
  int mate_pos=atoi(l.t[7]);
  int mapq=atoi(l.t[4]);
  GBamRecord *brec=wbam.new_record(l.t[0], flag, l.t[2], gpos, mapq,
          l.t[5], l.t[6], mate_pos, insert_size, l.t[9], l.t[10]);
  //now parse and add aux data:
  if (l.tcount>11) {
       for (int i=11;i<l.tcount;i++) {
         brec->add_aux(l.t[i]);
         }//for each aux field
       }
  wbam.write(brec);
  delete brec;
}

void driver_bam(string& fname, GBamWriter& bam_writer, GBamWriter* umbam) {
 tamFile fh=sam_open(sam_header.c_str());
 bam_header_t* header=sam_header_read(fh);
 sam_close(fh);
 priority_queue< pair<uint64_t,bam1_t*>,
					vector<pair<uint64_t, bam1_t*> >,
					BamMapOrdering > map_pq;
 GBamWriter* wmulti=NULL; //for multi-mapped prefiltering
 if (!aux_outfile.empty() && max_multihits>0) {
    wmulti=new GBamWriter(aux_outfile.c_str(), sam_header.c_str());
 }
 tamFile fp=sam_open(fname.c_str());
 bam1_t *b = bam_init1();
 //uint64_t last_id = 0;
 do {
   bool new_record = (sam_read1(fp, header, b) >= 0);
   if (new_record) {
     char *qname  = bam1_qname(b);
     uint64_t qid=(uint64_t)atol(qname);
     bam1_t* bamrec=bam_dup1(b);
     map_pq.push(make_pair(qid, bamrec));
   }

   while (map_pq.size() > 1000000 || (!new_record && map_pq.size() > 0)) {
     uint64_t rid=map_pq.top().first;
     bam1_t* tb=map_pq.top().second;
     bool unmapped = (tb->core.flag & BAM_FUNMAP) != 0;
     if (unmapped) { //unmapped read
       if (umbam!=NULL)
    	 umbam->write(tb);
       bam_destroy1(tb);
       map_pq.pop();
       }
     else { //mapped read
       vector<pair<uint64_t, bam1_t*> > read_hits;
      //all mappings of a read are dealt with here
      //if (!unmapped) {
       // mapped read
       //collect all hits for this read
       read_hits.push_back(map_pq.top());
       unsigned int mcount=0; //number of good scoring multi-hits
       int tbscore=0;
       if (bowtie2) {
          uint8_t* ptr = bam_aux_get(tb, "AS");
          if (ptr) {
        	    tbscore=bam_aux2i(ptr);
        	    if (tbscore>=bowtie2_min_score) {
        	       ++mcount;
        	       }
                }
          }
       else mcount++; //bowtie 1
       map_pq.pop();
       while (map_pq.size()>0 && map_pq.top().first==rid) {
    	 read_hits.push_back(map_pq.top());
    	 if (bowtie2) {
           uint8_t* ptr = bam_aux_get(map_pq.top().second, "AS");
           if (ptr) {
         	    int score=bam_aux2i(ptr);
         	    if (score>=bowtie2_min_score) {
         	       ++mcount;
         	       }
                 }
    	    }
    	 else mcount++;
    	 map_pq.pop();
       }
       int32_t num_hits=read_hits.size();
       if (wmulti && mcount>max_multihits) {
    	 if (num_hits>1) {
    	   bam_aux_append(tb, "NH", 'i', 4, (uint8_t*)&num_hits);
    	   }
    	 wmulti->write(tb);
       }
       else { //keep these mappings
       // In case of Bowtie2, some of the mapped reads against either transcriptome or genome
       // may have low alignment scores due to gaps, in which case we will remap those.
       // Later, we may have better alignments that usually involve splice junctions.
         if (bowtie2 && tbscore<bowtie2_min_score) {
			 //low score for best mapping, we want to map this read later as well
			 //unmapped = true;
			 if (umbam!=NULL) {
			   umbam->write(tb);
			   }
         }
        //-- write all hits for this read:
         for (vector<pair<uint64_t, bam1_t*> >::size_type i=0;i<read_hits.size();++i)
         {
           const pair<uint64_t, bam1_t*>& v = read_hits[i];
           if (num_hits>1)
	    	 bam_aux_append(v.second, "NH", 'i', 4, (uint8_t*)&num_hits);
           bam_writer.write(v.second, v.first);
         }
       }
       //free the read hits
       for (vector<pair<uint64_t, bam1_t*> >::size_type i=0;i<read_hits.size();++i)
       {
         const pair<uint64_t, bam1_t*>& v = read_hits[i];
         bam_destroy1(v.second);
       }
      } // mapped reads
   }
 } while (map_pq.size() > 0); //while SAM records

 bam_destroy1(b);
 bam_header_destroy(header);
 if (wmulti) delete wmulti;
}

void driver_headerless(FILE* map_file, FILE* f_out)
{
	
	char bwt_buf[4096];
	
	priority_queue< pair<uint64_t,TabSplitLine*>,
					vector<pair<uint64_t, TabSplitLine*> >,
					MapOrdering > map_pq;
	
	while (fgets(bwt_buf, sizeof(bwt_buf), map_file))
	{
		// Chomp the newline
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;
		if (*bwt_buf == 0)
			continue;
 		TabSplitLine* l=new TabSplitLine(bwt_buf);
		if (l->tcount<10 || l->t[0][0]=='@') {
		     delete l;
		     continue;
		     }
		//char* hitline = strdup(bwt_buf);
		uint64_t qid = (uint64_t)atol(l->t[0]);
		map_pq.push(make_pair(qid, l));
		
		if (map_pq.size() > 1000000)
		{
			const pair<uint64_t, TabSplitLine*>& t = map_pq.top();
			writeSamLine(*t.second, f_out);
			delete t.second;
			map_pq.pop();
		}
	}
	
	while (map_pq.size())
	{
		const pair<uint64_t, TabSplitLine*>& t = map_pq.top();
		writeSamLine(*t.second, f_out);
		delete t.second;
		map_pq.pop();
	}
}

void print_usage()
{
//    fprintf(stderr, "Usage:   fix_map_ordering <map.bwtout> [<reads.fa/.fq>]\n");
  fprintf(stderr,
    "Usage: \nfix_map_ordering [--sam-header <sam_header_file>] <in_SAM_file> <out_BAM_file> [<out_unmapped_BAM>]\n");
}

int main(int argc, char** argv)
{
    int parse_ret = parse_options(argc, argv, print_usage);
    if (parse_ret)
        return parse_ret;

    //if --sam_header option was given write BAM and expects (headerless) SAM input;
    //            else simply lets SAM lines through
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    string map_file_name = argv[optind++];
    string out_file_name("-");
    string out_unmapped_fname;
    if (optind<argc) {
       out_file_name = argv[optind++];
       }
    if (optind<argc) {
       out_unmapped_fname = argv[optind++];
       }
	//fix_map_ordering gets bowtie's SAM output at stdin, uncompressed
    //if the SAM input is headerless and no unmapped files
	if (sam_header.empty()) {
	 	//write only mapped reads as SAM lines
		FILE* f_out=NULL;
		if (out_file_name=="-") f_out=stdout;
		    else {
		      if ((f_out=fopen(out_file_name.c_str(),"w"))==NULL)
		         err_die("Error: cannot create output file %s\n",
		        		  out_file_name.c_str());
		      }
		FILE* f_in=NULL;
		if (map_file_name=="-") f_in=stdin;
		    else {
		      if ((f_in=fopen(map_file_name.c_str(),"r"))==NULL)
		         err_die("Error: cannot open input file %s\n",
		        		  map_file_name.c_str());
		      }
        driver_headerless(f_in, f_out);
	    }
	else {
		//BAM output, we have SAM header
		GBamWriter bam_writer(out_file_name.c_str(),sam_header.c_str(), index_outfile);
		GBamWriter *unmapped_bam_writer=NULL;
		if (!out_unmapped_fname.empty())
		    unmapped_bam_writer=new GBamWriter(out_unmapped_fname.c_str(),
		                                       sam_header.c_str());
		driver_bam(map_file_name, bam_writer, unmapped_bam_writer);
		delete unmapped_bam_writer;
	    }

    return 0;
}
