/*
 * Author: Harold Pimentel
 * Contact: http://cs.berkeley.edu/~pimentel
 * Date: June 10, 2011
 */
#include "map2gtf.h"
#include "tokenize.h"

void m2g_print_usage()
{
    std::cerr << "Usage: map2gtf\tannotation.gtf "
            << " alignments.bwtout out_file.sam" << std::endl;
}

Map2GTF::Map2GTF(const std::string& gtf_fname, const std::string& in_fname) :
    gtf_fname_(gtf_fname), in_fname_(in_fname), refSeqTable_(true)
{
  gtf_fhandle_ = fopen(gtf_fname_.c_str(), "r");
  if (gtf_fhandle_ == NULL)
    {
      std::cerr << "FATAL: Couldn't open annotation: " << gtf_fname_
                << std::endl;
      exit(1);
    }
  std::cout << "Reading the annotation file: " << gtf_fname_ << std::endl;
  gtfReader_.init(gtf_fhandle_, true); //only recognizable transcripts will be loaded
  gtfReader_.readAll();
  
  in_fhandle_=samopen(in_fname_.c_str(), "rb", 0);
  if (in_fhandle_ == NULL)
    {
      std::cerr << "FATAL: Couldn't open input bam file: " << in_fname_
		<< std::endl;
      exit(1);
    }
  in_sam_header_ = in_fhandle_->header;
  
  std::cout << "Done with initializtion. " << std::endl;
}

Map2GTF::~Map2GTF()
{
  std::cout << "map2gtf has completed. Cleaning up." << std::endl;
  if (gtf_fhandle_ != NULL && fclose(gtf_fhandle_))
    {
      std::cerr << "Warning: Error closing annotation: " << gtf_fname_
                << std::endl;
    }
  
  if (in_fhandle_ != NULL)
    {
      samclose(in_fhandle_);
    }
  
  std::cout << "Done. Thanks!" << std::endl;
}

//
bool Map2GTF::next_read_hits(vector<bam1_t*>& hits, size_t& num_hits)
{
  long read_id = 0;
  if (hits.size() > num_hits)
    {
      bam1_t* temp = hits[num_hits];
      hits[num_hits] = hits.front();
      hits.front() = temp;

      num_hits = 1;
      char* name = bam1_qname(hits.front());
      read_id = atol(name);
    }
  else
    num_hits = 0;

  while (true)
    {
      bam1_t* hit = NULL;
      if (num_hits >= hits.size())
	hits.push_back(bam_init1());

      hit = hits[num_hits];

      if (samread(in_fhandle_, hit) <= 0)
	{
	  for (size_t i = num_hits; i < hits.size(); ++i)
	      bam_destroy1(hits[i]);

	  hits.erase(hits.begin() + num_hits, hits.end());
	  break;
	}

      char* name = bam1_qname(hit);
      long temp_read_id = atol(name);
      if (num_hits == 0)
	{
	  read_id = temp_read_id;
	}
      else if (read_id != temp_read_id)
	{
	  break;
	}
	
      ++num_hits;
    }
  
  return num_hits > 0;
}

void Map2GTF::convert_coords(const std::string& out_fname, const std::string& sam_header)
{
  samfile_t* out_sam_header_file = samopen(sam_header.c_str(), "r", 0);
  if (out_sam_header_file == NULL)
    std::cerr << "Error opening sam header: " << sam_header << std::endl;

  out_sam_header_ = out_sam_header_file->header;
  
  samfile_t* bam_writer = samopen(out_fname.c_str(), "wb", out_sam_header_);
  if (bam_writer == NULL)
    std::cerr << "Error opening bam file for writing: " << out_fname << std::endl;

  ref_to_id_.clear();
  for (int i = 0; i < out_sam_header_->n_targets; ++i)
    {
      ref_to_id_[out_sam_header_->target_name[i]] = i;
    }
  
  std::vector<TranscriptomeHit> read_list;
  GffObj* p_trans = NULL;
  
  HitsForRead hit_group;
  std::vector<TranscriptomeHit>::iterator bh_it;
  std::vector<TranscriptomeHit>::iterator bh_unique_it;
  BowtieHit bwt_hit;

  vector<bam1_t*> hits;
  size_t num_hits = 0;
  // a hit group is a set of reads with the same name
  while (next_read_hits(hits, num_hits))
    {
      for (size_t i = 0; i < num_hits; ++i)
        {
	  bam1_t* hit = hits[i];
	  const char* target_name = in_sam_header_->target_name[hit->core.tid];
	  size_t trans_idx = atoi(target_name);
	    
	  p_trans = gtfReader_.gflst.Get(trans_idx);
	  TranscriptomeHit converted_out(hit, p_trans);
	  trans_to_genomic_coords(converted_out);
	  
	  read_list.push_back(converted_out);
        }
      
      // XXX: Fine for now... should come up with a more efficient way though
      // FIXME: Take frag length into consideration when filtering
      std::sort(read_list.begin(), read_list.end());
      bh_unique_it = std::unique(read_list.begin(), read_list.end());
      for (bh_it = read_list.begin(); bh_it != bh_unique_it; ++bh_it)
        {
	  samwrite(bam_writer, bh_it->hit);
        }
      read_list.clear();
    }

  for (size_t i = 0; i < hits.size(); ++i)
    {
      bam_destroy1(hits[i]);
    }
  hits.clear();

  samclose(bam_writer);
}

void Map2GTF::trans_to_genomic_coords(TranscriptomeHit& hit)
//out.trans must already have the corresponding GffObj*
{
    // read start in genomic coords
    size_t read_start = 0;

    GList<GffExon>& exon_list = hit.trans->exons;
    GffExon* cur_exon;
    GffExon* next_exon;
    int cur_pos;
    int match_length;
    int miss_length;
    int cur_intron_len = 0;
    int i = 0;

    static const int MAX_CIGARS = 1024;
    int cigars[MAX_CIGARS];
    int num_cigars = 0;

    // TODO: Check this return value
    bool ret_val = get_read_start(&exon_list, hit.hit->core.pos, read_start, i);
    cur_pos = read_start;
    for (int c = 0; c < hit.hit->core.n_cigar; ++c)
      {
	int opcode = bam1_cigar(hit.hit)[c] & BAM_CIGAR_MASK;
	int length = bam1_cigar(hit.hit)[c] >> BAM_CIGAR_SHIFT;

	if (opcode == BAM_CINS)
	  {
	    cigars[num_cigars] = opcode | (length << BAM_CIGAR_SHIFT);
	    ++num_cigars;
	  }
	
	if (opcode != BAM_CMATCH && opcode != BAM_CDEL)
	  continue;
	
	int remaining_length = length;
	for (; i < exon_list.Count(); ++i)
	  {
	    cur_exon = exon_list.Get(i);
	    if (cur_pos >= (int)cur_exon->start &&
		cur_pos + remaining_length - 1 <= (int)cur_exon->end) // read ends in this exon
	      {
		cigars[num_cigars] = opcode | (remaining_length << BAM_CIGAR_SHIFT);
		++num_cigars;
		cur_pos += remaining_length;
		break;
	      }
	    
	    // shouldn't need the check... can switch to a regular "else"
	    else if (cur_pos >= (int)cur_exon->start &&
		     cur_pos + remaining_length - 1 > (int)cur_exon->end)// read is spliced and overlaps this exon
	      {
		// XXX: This should _never_ go out of range.
		// get the max length that fits in this exon, go to next exon
		// cur_pos should be the next exon start
		// set assertion to check this
		
		// TODO: check this
		match_length = cur_exon->end - cur_pos + 1;
		if (match_length > 0)
		  {
		    cigars[num_cigars] = opcode | (match_length << BAM_CIGAR_SHIFT);
		    ++num_cigars;
		  }
		
		// XXX: DEBUG
		if (i + 1 >= exon_list.Count())
		  {
		    std::cerr << "trying to access: " << i + 1 << " when size is: "
			      << exon_list.Count() << std::endl;
		    print_trans(hit.trans, hit.hit, remaining_length, match_length, cur_pos,
				read_start);
		    exit(1);
		  }
		
		else
		  next_exon = exon_list.Get(i + 1);
		
		// and this
		miss_length = next_exon->start - cur_exon->end - 1;
		cur_intron_len += miss_length;

		cigars[num_cigars] = BAM_CREF_SKIP | (miss_length << BAM_CIGAR_SHIFT);
		++num_cigars;
		
		cur_pos += match_length + miss_length;
		remaining_length -= match_length;
		assert(cur_pos == next_exon->start);
	      }
	  }
      }

    hit.hit->core.tid = ref_to_id_[hit.trans->getRefName()];
    hit.hit->core.pos = read_start - 1;
    hit.hit->core.flag &= ~BAM_FSECONDARY;
    
    int old_n_cigar = hit.hit->core.n_cigar;
    if (num_cigars != old_n_cigar)
      {
	int data_len = hit.hit->data_len + 4 * (num_cigars - old_n_cigar);
	int m_data = max(data_len, hit.hit->m_data);
	kroundup32(m_data);
	
	uint8_t* data = (uint8_t*)calloc(m_data, 1);

	int copy1_len = (uint8_t*)bam1_cigar(hit.hit) - hit.hit->data;
	memcpy(data, hit.hit->data, copy1_len);

	int copy2_len = num_cigars * 4;
	memcpy(data + copy1_len, cigars, copy2_len);

	int copy3_len = hit.hit->data_len - copy1_len - (old_n_cigar * 4);
	memcpy(data + copy1_len + copy2_len, bam1_seq(hit.hit), copy3_len);

	hit.hit->core.n_cigar = num_cigars;
	
	free(hit.hit->data);
	hit.hit->data = data;
	hit.hit->data_len = data_len;
	hit.hit->m_data = m_data;
      }

    char strand = hit.trans->strand;
    uint8_t* ptr = bam_aux_get(hit.hit, "XS");
    if (ptr)
      bam_aux_del(hit.hit, ptr);

    if (strand == '+' || strand == '-')
      bam_aux_append(hit.hit, "XS", 'A', 1, (uint8_t*)&strand);
}

void print_trans(GffObj* trans, const bam1_t* in, size_t rem_len,
		 size_t match_len, size_t cur_pos, size_t start_pos)
{
    GffExon* p_exon;
    std::cerr << "\tCur_pos: " << cur_pos << " remaining: " << rem_len
            << " match_len: " << match_len << std::endl;
    std::cerr << "\tTranscript:\t" << trans->start << "\t" << trans->end
            << std::endl;
    for (int i = 0; i < trans->exons.Count(); ++i)
    {
        p_exon = trans->exons.Get(i);
        std::cerr << "\t\t" << p_exon->start << "\t" << p_exon->end
                << std::endl;
    }
    std::cerr << std::endl;

    std::cerr << "Read_id: " << bam1_qname(in) << std::endl;
    std::cerr << "\tgff_start: " << in->core.pos << "\tgenome_start: "
            << start_pos << std::endl;
}

// Returns false if not in this exon list
bool get_read_start(GList<GffExon>* exon_list, size_t gtf_start,
        size_t& genome_start, int& exon_idx)
{
    GffExon* cur_exon;
    size_t cur_intron_dist = 0;
    size_t trans_start = exon_list->First()->start;
    int trans_offset = 0;
    for (int i = 0; i < exon_list->Count(); ++i)
    {
        cur_exon = exon_list->Get(i);
        trans_offset = trans_start + cur_intron_dist;

        if (gtf_start >= cur_exon->start - trans_offset && gtf_start
                <= cur_exon->end - trans_offset)
        {
            genome_start = gtf_start + trans_start + cur_intron_dist;
            exon_idx = i;
            return true;
        }
        else
        {
            if (i + 1 < exon_list->Count())
                cur_intron_dist += exon_list->Get(i + 1)->start - cur_exon->end
                        - 1;
            else
                return false;
        }
    }

    return false;
}

int main(int argc, char *argv[])
{
  int parse_ret = parse_options(argc, argv, m2g_print_usage);
  if (parse_ret)
    return parse_ret;
  
  if (optind >= argc)
    {
      m2g_print_usage();
      return 1;
    }
  
  std::string gtf_file(argv[optind++]);
  std::string in_fname(argv[optind++]);
  std::string out_fname(argv[optind++]);
  
  if (gtf_file == "" || in_fname == "" || out_fname == "")
    {
      m2g_print_usage();
      exit(1);
    }
  
  Map2GTF gtfMapper(gtf_file, in_fname);
  gtfMapper.convert_coords(out_fname, sam_header);
  
  return 0;
}
