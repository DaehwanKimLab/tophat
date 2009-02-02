/*
 *  genes.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 2/1/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "genes.h"
#include "FSA/gff.h"
#include "FSA/sequence.h"

using namespace fsa;
using namespace std;

bool exon_lt_by_pos(const Exon* lhs, const Exon* rhs)
{
    return lhs->end() <= rhs->start();
}

Gene* Exon::gene()
{
    assert(_transcript);
    return _transcript->_gene;
}

uint32_t Gene::exonic_length() const
{
    uint32_t gene_exonic_length = 0;
    
    pair<uint32_t, uint32_t> gene_coords = coords();
    uint32_t start = gene_coords.first;
    uint32_t end = gene_coords.second;
    
    vector<bool> exonic_coords(end - start + 1);
    for (map<std::string, Transcript*>::const_iterator t_itr = transcripts.begin();
         t_itr != transcripts.end();
         ++t_itr)
    {
        Transcript& mRNA = *(t_itr->second);
        for (vector<Exon*>::const_iterator exon_itr = mRNA.exons.begin();
             exon_itr != mRNA.exons.end();
             ++exon_itr)
        {
            const Exon* exon = *exon_itr;
            assert (exon->start() >= start);
            for (uint32_t i = exon->start();
                 i < exon->end();
                 ++i)
            {
                // Did we already count this one?
                if (!exonic_coords[i - start] /*&&
                 (!ref_str || !ref_str->is_pos_hardmasked(i - 1))*/)
                {
                    gene_exonic_length++;
                }
                exonic_coords[i - start] = true;
            }
        }
    }
    return gene_exonic_length;
}

pair<uint32_t, uint32_t> Gene::coords() const
{
    uint32_t start = 0xFFFFFFFF;
    uint32_t end = 0;
    for (map<std::string, Transcript*>::const_iterator t_itr = transcripts.begin();
         t_itr != transcripts.end();
         ++t_itr)
    {
        if (t_itr->second->exons.front()->start() < start)
        {
            start = t_itr->second->exons.front()->start();
        }
        if (t_itr->second->exons.back()->end() < end)
        {
            end = t_itr->second->exons.back()->end();
        }
    }
    return make_pair<uint32_t, uint32_t>(start, end);
}

uint32_t Gene::exonic_depth(const vector<short>& DoC) const
{
    uint32_t gene_DoC = 0;
    
    pair<uint32_t, uint32_t> gene_coords = coords();
    uint32_t start = gene_coords.first;
    uint32_t end = gene_coords.second;
    
    vector<bool> exonic_coords(end - start + 1);
    for (map<std::string, Transcript*>::const_iterator t_itr = transcripts.begin();
         t_itr != transcripts.end();
         ++t_itr)
    {
        Transcript& mRNA = *(t_itr->second);
        for (vector<Exon*>::const_iterator exon_itr = mRNA.exons.begin();
             exon_itr != mRNA.exons.end();
             ++exon_itr)
        {
            const Exon* exon = *exon_itr;
            assert (exon->start() >= start);
            for (uint32_t i = exon->start();
                 i < exon->end();
                 ++i)
            {
                if (i - 1 >= DoC.size())
                    break;
                // Did we already count this one?
                if (!exonic_coords[i - start] /*&&
                    (!ref_str || !ref_str->is_pos_hardmasked(i - 1))*/)
                {
                    gene_DoC += DoC[i - 1];
                }
                
                exonic_coords[i - start] = true;
            }
        }
    }
    return gene_DoC;
}

Expression* Gene::expression(const vector<short>& DoC, 
                             uint32_t total_map_depth) const 
{
    uint32_t length = exonic_length();
    uint32_t depth = exonic_depth(DoC);
    double gene_avg_DoC = depth / (double) length;
    
    double gene_rpkm = 1000000000 * (gene_avg_DoC / (double)total_map_depth);
    return new Expression(gene_rpkm, 0.0);
}

void GeneFactory::get_genes(const GFF_database& gffdb, GeneTable& genes)
{
    map<string, Transcript*> mRNAs;
    
    for(GFF_database::const_iterator gff_itr = gffdb.begin();
        gff_itr != gffdb.end();
        ++gff_itr)
    {
        const GFF& gff_rec = *gff_itr;
        if (gff_rec.type == "gene")
        {
            GFF::AttributeTable::const_iterator att_itr;
            att_itr = gff_rec.attributes.find("ID");
            if (att_itr == gff_rec.attributes.end() ||
                att_itr->second.size() != 1)
            {
                cerr << "Malformed gene record " << gff_rec << endl;
                continue;
            }
            const string& id = att_itr->second.front();
            
            Gene* gene = new Gene();
            genes[id] = gene;
        }
    }
    
    for(GFF_database::const_iterator gff_itr = gffdb.begin();
        gff_itr != gffdb.end();
        ++gff_itr)
    {
        const GFF& gff_rec = *gff_itr;
        if (gff_rec.type == "mRNA")
        {
            GFF::AttributeTable::const_iterator att_itr;
            att_itr = gff_rec.attributes.find("ID");
            if (att_itr == gff_rec.attributes.end() ||
                att_itr->second.size() != 1)
            {
                cerr << "Malformed transcript record " << gff_rec << endl;
                continue;
            }
            const string& id = att_itr->second.front();
            
            att_itr = gff_rec.attributes.find("Parent");
            if (att_itr == gff_rec.attributes.end() ||
                att_itr->second.size() != 1)
            {
                cerr << "Malformed transcript record " << gff_rec << endl;
                continue;
            }
            const string& gene_id = att_itr->second.front();
            
            GeneTable::iterator gene_itr = genes.find(gene_id);
            if (gene_itr != genes.end())
            {
                Transcript* mRNA = new Transcript();
                mRNA->_gene = gene_itr->second;
                gene_itr->second->transcripts[id] = mRNA;
                mRNAs[id] = mRNA;
            }
            else
            {
                cerr << "No gene record "  << gene_id 
                     << " for transcript " << id << endl;
                continue;
            }
        }
    }
    
    for(GFF_database::const_iterator gff_itr = gffdb.begin();
        gff_itr != gffdb.end();
        ++gff_itr)
    {
        const GFF& gff_rec = *gff_itr;
        if (gff_rec.type == "exon")
        {
            GFF::AttributeTable::const_iterator att_itr;
            att_itr = gff_rec.attributes.find("Parent");
            if (att_itr == gff_rec.attributes.end())
            {
                cerr << "Malformed exon record " << gff_rec << endl;
                continue;
            }
            const string& mRNA_id = att_itr->second.front();
            
            att_itr = gff_rec.attributes.find("ID");
            if (att_itr == gff_rec.attributes.end() ||
                att_itr->second.size() != 1)
            {
                cerr << "Malformed transcript record " << gff_rec << endl;
                continue;
            }
            const string& id = att_itr->second.front();
            
            map<string, Transcript*>::iterator mRNA_itr = mRNAs.find(mRNA_id);
            if (mRNA_itr != mRNAs.end())
            {
                Exon* exon = new Exon();
                exon->_start = gff_rec.start;
                exon->_end = gff_rec.end;
                exon->_transcript = mRNA_itr->second;
                mRNA_itr->second->exons.push_back(exon);
            }
            else
            {
                cerr << "No transcript record "  << mRNA_id 
                << " for exon " << id << endl;
                continue;
            }
            
        }
    }
    
    for (map<string,Transcript*>::iterator itr = mRNAs.begin(); 
         itr != mRNAs.end();
         ++itr)
    {
        Transcript& mRNA = *(itr->second);
        sort(mRNA.exons.begin(), mRNA.exons.end(),exon_lt_by_pos); 
    }
}

uint32_t total_exonic_depth(const GeneTable& genes,
                            const vector<short>& DoC)
{
    uint32_t total_depth = 0;
    for (GeneTable::const_iterator gene_itr = genes.begin();
         gene_itr != genes.end();
         ++gene_itr)
    {
        total_depth += gene_itr->second->exonic_depth(DoC);
    }
    
    return total_depth;
}

void print_gene_expression(FILE* expr_out, 
                           const GeneTable& genes,
                           const vector<short>& DoC)
{
    
    uint32_t total_depth = total_exonic_depth(genes, DoC);
    for (GeneTable::const_iterator itr = genes.begin();
         itr != genes.end();
         ++itr)
    {
        Expression* expr = itr->second->expression(DoC, total_depth);
        fprintf(expr_out, "%s\t%lf\t%lf\n", itr->first.c_str(), expr->rpkm, expr->mend);
    }
}
