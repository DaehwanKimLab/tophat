/*
 * Author: Harold Pimentel
 * Contact: http://cs.berkeley.edu/~pimentel
 * Date: June 10, 2011
 */

#include "map2gtf.h"

void m2g_print_usage()
{
    std::cerr << "Usage: map2gtf\tannotation.gtf "
            << " alignments.bwtout out_file.sam" << std::endl;
}

Map2GTF::Map2GTF(std::string gtf_fname, std::string reads_fname) :
    gtf_fname_(gtf_fname), reads_fname_(reads_fname), refSeqTable_(true, true)
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

    std::cout << "Initializing the SAMHitFactory." << std::endl;
    hitFactory_ = new BowtieHitFactory(readTable_, refSeqTable_);
    std::cout << "Done with init samfactory" << std::endl;

    // TODO: Use correct piping depending on file format (i.e. (don't) use zpacker)
    if (zpacker != "")
    {
        std::string pipe_cmd = getUnpackCmd(reads_fname, true);
        if (!pipe_cmd.empty())
        {
            pipe_cmd.append(" -cd");
        }
        else
        {
            std::cerr << "Couldn't get the correct unzipper" << std::endl;
            exit(1);
        }
        reads_pipe_.openRead(reads_fname_, pipe_cmd);
        if (reads_pipe_.file == NULL)
        {
            std::cerr << "FATAL: Couldn't open read file: " << reads_fname_
                    << std::endl;
            exit(1);
        }
        std::cout << "Initializing the HitStream." << std::endl;
        hitStream_ = new HitStream(reads_pipe_, hitFactory_, false, true,
                false, true, true, true);
    }
    else // no zpacker, open regular file
    {
        std::cerr << "Not using compression!" << std::endl;
        reads_fhandle_ = fopen(reads_fname_.c_str(), "r");
        if (reads_fhandle_ == NULL)
        {
            std::cerr << "FATAL: Couldn't open annotation: " << gtf_fname_
                    << std::endl;
            exit(1);
        }

        std::cout << "Initializing the HitStream." << std::endl;
        hitStream_ = new HitStream(reads_fhandle_, hitFactory_, false, true,
                false, true, true, true);
    }
    
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

    if (zpacker != "")
    {
        std::cout << "Closing the zpacker. " << std::endl;
        reads_pipe_.close();
    }
    else
    {
        if (reads_fhandle_ != NULL && fclose(reads_fhandle_))
        {
            std::cerr << "Warning: Error closing reads: " << reads_fname_
            << std::endl;
        }
    }

    delete hitFactory_;
    delete hitStream_;
    std::cout << "Done. Thanks!" << std::endl;
}

void Map2GTF::convert_coords(std::string out_fname, std::string sam_header)
{
    // FIXME: come up with a more efficient layout for doing this
    GBamWriter bam_writer(out_fname.c_str(), sam_header.c_str());

    std::vector<TranscriptomeHit> read_list;

    GffObj* p_trans = NULL;

    HitsForRead hit_group;
    char read_name[MAX_READ_NAME_LEN]; // FIXME: use the macro size

    std::vector<TranscriptomeHit>::iterator bh_it;
    std::vector<TranscriptomeHit>::iterator bh_unique_it;
    BowtieHit bwt_hit;
    // a hit group is a set of reads with the same name
    while (hitStream_->next_read_hits(hit_group))
    {
        for (size_t i = 0; i < hit_group.hits.size(); ++i)
        {
            // TODO: check if read is unmapped.
            // if so, figure out what to do with it
            bwt_hit = hit_group.hits[i];

            sprintf(read_name, "%u", bwt_hit.insert_id());

            size_t trans_idx = atoi(refSeqTable_.get_name(bwt_hit.ref_id()));
            p_trans = gtfReader_.gflst.Get(trans_idx);
            TranscriptomeHit converted_out(p_trans);
            trans_to_genomic_coords(read_name, hitFactory_, bwt_hit,
                    converted_out);
            read_list.push_back(converted_out);
        }

        // XXX: Fine for now... should come up with a more efficient way though
        // FIXME: Take frag length into consideration when filtering
        std::sort(read_list.begin(), read_list.end());
        bh_unique_it = std::unique(read_list.begin(), read_list.end());

        for (bh_it = read_list.begin(); bh_it != bh_unique_it; ++bh_it)
        {
            print_bamhit(bam_writer, read_name, bh_it->hit, bh_it->trans->getRefName(),
                    bh_it->hit.seq().c_str(), bh_it->hit.qual().c_str(), true);
        }
        read_list.clear();
    }

    //    fclose(out_sam);
}

void trans_to_genomic_coords(const char* read_name, HitFactory* hitFactory,
        const BowtieHit& in, TranscriptomeHit& out)
   //out.trans must already have the corresponding GffObj*
{
    std::string read_id(read_name);
    std::string ref_name(out.trans->getRefName());
    bool spliced = false;
    std::vector<CigarOp> cig_list;

    // read start in genomic coords
    size_t read_start = 0;

    GList<GffExon>& exon_list = out.trans->exons;
    GffExon* cur_exon;
    GffExon* next_exon;
    size_t cur_pos;
    size_t match_length;
    size_t miss_length;
    size_t remaining_length = static_cast<size_t> (in.read_len());
    size_t cur_intron_len = 0;
    int i = 0;

    // TODO: Check this return value
    bool ret_val = get_read_start(&exon_list, in.left(), read_start, i);
    cur_pos = read_start;

    // FIXME: Start is off by one
    for (; i < exon_list.Count(); ++i)
    {
        cur_exon = exon_list.Get(i);

        assert(cur_pos > 0 && cur_pos + in.read_len() > 0);

        if (cur_pos >= cur_exon->start && cur_pos + remaining_length - 1
                <= cur_exon->end) // read ends in this exon
        {
            CigarOp match_cig(MATCH, remaining_length);
            cig_list.push_back(match_cig);
            break;
        }
        // shouldn't need the check... can switch to a regular "else"
        else if (cur_pos >= cur_exon->start && cur_pos + remaining_length - 1
                > cur_exon->end)// read is spliced and overlaps this exon
        {
            spliced = true;
            // XXX: This should _never_ go out of range.
            // get the max length that fits in this exon, go to next exon
            // cur_pos should be the next exon start
            // set assertion to check this

            // TODO: check this
            match_length = cur_exon->end - cur_pos + 1;
            CigarOp match_cig(MATCH, match_length);
            cig_list.push_back(match_cig);

            // XXX: DEBUG
            if (i + 1 >= exon_list.Count())
            {
                std::cerr << "trying to access: " << i + 1 << " when size is: "
                        << exon_list.Count() << std::endl;
                print_trans(out.trans, in, remaining_length, match_length, cur_pos,
                        read_start);
                exit(1);
            }

            else
                next_exon = exon_list.Get(i + 1);

            // and this
            miss_length = next_exon->start - cur_exon->end - 1;
            cur_intron_len += miss_length;
            CigarOp miss_cig(REF_SKIP, miss_length);
            cig_list.push_back(miss_cig);

            cur_pos += match_length + miss_length;
            remaining_length -= match_length;
            assert(cur_pos == next_exon->start);
        }
    }
    bool antisense_splice = (spliced && out.trans->strand=='-'); //transcript strand <=> splice strand (if spliced)
    read_start -= 1; // handle the off-by-one problem
    out.hit = hitFactory->create_hit(read_id, ref_name,
            static_cast<int> (read_start), cig_list, in.antisense_align(),
            antisense_splice, in.edit_dist(), in.splice_mms(), in.end());
    out.hit.seq(in.seq());
    out.hit.qual(in.qual());
}

void print_trans(GffObj* trans, const BowtieHit& in, size_t rem_len,
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

    std::cerr << "Read_id: " << in.insert_id() << std::endl;
    std::cerr << "\tgff_start: " << in.left() << "\tgenome_start: "
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
    std::string reads_file(argv[optind++]);
    std::string out_fname(argv[optind++]);

    if (gtf_file == "" || reads_file == "" || out_fname == "")
    {
        m2g_print_usage();
        exit(1);
    }

    Map2GTF gtfMapper(gtf_file, reads_file);
    //    gtfMapper.convert_coords(out_fname);
    gtfMapper.convert_coords(out_fname, sam_header);

    return 0;
}
