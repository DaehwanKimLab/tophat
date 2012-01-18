//
//  gtfToFasta.cpp
//  TopHat
//
//  Created by Harold Pimentel on 10/26/11.
//
#include "GTFToFasta.h"


std::string get_exonic_sequence(GffObj *p_trans,
                                FastaRecord *rec)
{
    // TODO: Write me.
    
    // TODO: Should I check if the transcript name matches the record name?
    //       For now, leave this up to the caller. 
    GList<GffExon>& exon_list = p_trans->exons;
    GffExon* cur_exon;
    std::string exon_seq("");
    size_t length;
    
    for (int i = 0; i < exon_list.Count(); ++i) {
        cur_exon = exon_list.Get(i);
        length = cur_exon->end - cur_exon->start + 1;
        exon_seq += rec->seq_.substr(cur_exon->start - 1, length);
    }
    
    
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
    GffObj *p_trans;
    
    std::vector<size_t> *p_contig_vec;
    
    FastaReader fastaReader(genome_fname_);
    FastaWriter fastaWriter(out_fname);
    
    FastaRecord *cur_contig;
    
    while (fastaReader.good()) {
        cur_contig = new FastaRecord();
        fastaReader.next(cur_contig);

        // If this contig isn't in the map, then there are no transcripts
        // associated with it. Skip it.
        if (contigTransMap_.find(cur_contig->id_) ==
            contigTransMap_.end())
        {
            delete cur_contig;                
            continue;
        }
        
        p_contig_vec = contigTransMap_[cur_contig->id_.c_str()];
        std::string exon_seq;
        
        FastaRecord out_rec;
        for (size_t i = 0; i < p_contig_vec->size(); ++i) {
            size_t trans_idx = (*p_contig_vec)[i];            
            p_trans = gtfReader_.gflst.Get(trans_idx);
            exon_seq = get_exonic_sequence(p_trans, cur_contig);
            if (exon_seq.empty()) continue;
            std::stringstream ss;
            ss << trans_idx;
            out_rec.id_ = ss.str();
            ss.str(std::string()); //clear ss
            ss << p_trans->getGSeqName() << ':' << p_trans->start << '-' << p_trans->end << ' ' << p_trans->getID();
            //out_rec.desc_ = "";
            out_rec.desc_ = ss.str();
            out_rec.seq_ = exon_seq;
            fastaWriter.write(&out_rec);
        }
        delete cur_contig;        
    }
}

void GTFToFasta::transcript_map()
{
    GffObj *p_gffObj;
    const char *p_contig_name;
    std::vector<size_t> *p_contig_vec;
    
    for (int i = 0; i < gtfReader_.gflst.Count(); ++i) 
    {        
        p_gffObj = gtfReader_.gflst.Get(i);
        p_contig_name = p_gffObj->getRefName();
        std::string contig_name(p_contig_name);
        
        // Check if the current contig exists in the map
        // If it doesn't, add it        
        if (contigTransMap_.find(contig_name) == contigTransMap_.end())
        {
            p_contig_vec = new std::vector<size_t>;
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
