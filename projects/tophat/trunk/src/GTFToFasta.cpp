//
//  gtfToFasta.cpp
//  TopHat
//
//  Created by Harold Pimentel on 10/26/11.
//
#include "GTFToFasta.h"

std::string get_exonic_sequence(GffObj &p_trans,
                                FastaRecord &rec, std::string& coords)
{

    GList<GffExon>& exon_list = p_trans.exons;

    std::string exon_seq("");
    size_t length;
    coords.clear();
    std::stringstream ss;
    for (int i = 0; i < exon_list.Count(); ++i) {
      GffExon& cur_exon = *(exon_list.Get(i));
      length = cur_exon.end - cur_exon.start + 1;
      exon_seq += rec.seq_.substr(cur_exon.start - 1, length);
      ss << ',' << cur_exon.start << '-' << cur_exon.end;
    }
    coords = ss.str().substr(1);
    return exon_seq;
}


GTFToFasta::GTFToFasta(std::string gtf_fname, std::string genome_fname)
: genome_fhandle_(genome_fname.c_str(), false)
{
    gtf_fname_ = gtf_fname;
    gtf_fhandle_ = fopen(gtf_fname_.c_str(), "r");
    if (gtf_fhandle_ == NULL)
    {
        std::cerr << "FATAL: Couldn't open annotation: " << gtf_fname_
        << std::endl;
        exit(1);
    }
    std::cout << "Reading the annotation file: " << gtf_fname_ << std::endl;
    gtfReader_.init(gtf_fhandle_, true); //load recognizable transcript features only
    gtfReader_.readAll();

    genome_fname_ = genome_fname;

    // Make a map from the GffObj
    transcript_map();
}

GTFToFasta::~GTFToFasta()
{
    ContigTransMap::iterator it;

    for (it = contigTransMap_.begin(); it != contigTransMap_.end(); ++it) {
        delete it->second;
    }

}

void GTFToFasta::make_transcriptome(std::string out_fname)
{
    std::vector<int> *p_contig_vec;

    FastaReader fastaReader(genome_fname_);
    FastaWriter fastaWriter(out_fname);
    std::string tlst_fname(out_fname);
    tlst_fname.append(".tlst");
    std::ofstream tlst(tlst_fname.c_str());
    FastaRecord cur_contig;
    while (fastaReader.good()) {

        fastaReader.next(cur_contig);

        // If this contig isn't in the map, then there are no transcripts
        // associated with it. Skip it.
        if (contigTransMap_.find(cur_contig.id_) ==
            contigTransMap_.end())
        {
            continue;
        }

        p_contig_vec = contigTransMap_[cur_contig.id_];

        FastaRecord out_rec;
        for (size_t i = 0; i < p_contig_vec->size(); ++i) {
            int trans_idx = (*p_contig_vec)[i];
            GffObj *p_trans = gtfReader_.gflst.Get(trans_idx);
            //if (p_trans->isDiscarded() || p_trans->exons.Count()==0) continue;
            std::string coordstr;
            out_rec.seq_ = get_exonic_sequence(*p_trans, cur_contig, coordstr);
            if (out_rec.seq_.empty()) continue;
            std::stringstream ss;
            ss << trans_idx;
            out_rec.id_ = ss.str();
            //ss.str(std::string()); //clear ss
            out_rec.desc_=p_trans->getID();
            out_rec.desc_.push_back(' ');
            //ss << p_trans->getID() << ' ' << p_trans->getGSeqName() << p_trans->strand << '\t' << coordstr ;
            out_rec.desc_.append(cur_contig.id_);
            out_rec.desc_.push_back(p_trans->strand);
            out_rec.desc_.push_back(' ');
            out_rec.desc_.append(coordstr); //list of exon coordinates
            tlst << out_rec.id_ << ' ' << out_rec.desc_ << std::endl;
            //out_rec.desc_ = "";
            //out_rec.desc_ = ss.str();
            //out_rec.seq_ = exon_seq;
            fastaWriter.write(out_rec);
        }
    }
    tlst.close();
}

void GTFToFasta::transcript_map()
{
    GffObj *p_gffObj;
    const char *p_contig_name;
    std::vector<int> *p_contig_vec;

    for (int i = 0; i < gtfReader_.gflst.Count(); ++i) 
    {
        p_gffObj = gtfReader_.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count()==0)
           continue;

        p_contig_name = p_gffObj->getRefName();
        std::string contig_name(p_contig_name);

        // Check if the current contig exists in the map
        // If it doesn't, add it
        if (contigTransMap_.find(contig_name) == contigTransMap_.end())
        {
            p_contig_vec = new std::vector<int>;
            contigTransMap_[contig_name] = p_contig_vec;
        }
        else
        {
            p_contig_vec = contigTransMap_[contig_name];
        }

        p_contig_vec->push_back(i);
    }
}

void GTFToFasta::print_mapping()
{
    std::ofstream out_file("out.names");
    GffObj *p_gffObj;

    for (int i = 0; i < gtfReader_.gflst.Count(); ++i) 
    {
        p_gffObj = gtfReader_.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count()==0) continue;
        out_file << i << "\t" << p_gffObj->getID() << std::endl;
    }

    out_file.close();
}

void gtf2fasta_print_usage() 
{
    std::cerr << "Usage: gtf_to_fasta transcripts.gtf genome.fa out_file" << std::endl;
}

int main(int argc, char *argv[])
{
    int parse_ret = parse_options(argc, argv, gtf2fasta_print_usage);
    if (parse_ret)
        return parse_ret;
    
    if (optind >= argc)
    {
        gtf2fasta_print_usage();
        return 1;
    }
    
    std::string gtf_fname(argv[optind++]);
    std::string genome_fname(argv[optind++]);
    std::string out_fname(argv[optind++]);

    GTFToFasta gtfToFasta(gtf_fname, genome_fname);
    gtfToFasta.make_transcriptome(out_fname);
    //gtfToFasta.print_mapping();
    
    return 0;
}
