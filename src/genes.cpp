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

uint32_t Gene::exonic_length(const fsa::Sequence* ref_str) const
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
                 i <= exon->end();
                 ++i)
            {
                // Did we already count this one?
                if (!exonic_coords[i - start] &&
                 (!ref_str || !ref_str->is_pos_hardmasked(i - 1)))
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
        assert(t_itr->second->exons.size() > 0);
            
        if (t_itr->second->exons.front()->start() < start)
        {
            start = t_itr->second->exons.front()->start();
        }
        if (t_itr->second->exons.back()->end() > end)
        {
            end = t_itr->second->exons.back()->end();
        }
    }
    return make_pair<uint32_t, uint32_t>(start, end);
}

uint32_t Gene::exonic_depth(const vector<short>& DoC,
                            const fsa::Sequence* ref_str) const
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
                 i <= exon->end();
                 ++i)
            {
                if (i - 1 >= DoC.size())
                    break;
                // Did we already count this one?
                if (!exonic_coords[i - start] &&
                    (!ref_str || !ref_str->is_pos_hardmasked(i - 1)))
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
                             uint32_t total_map_depth,
                             const fsa::Sequence* ref_str) const 
{
    uint32_t length = exonic_length(ref_str);
    if (length == 0)
        return NULL;
    
    uint32_t depth = exonic_depth(DoC,ref_str);
    double gene_avg_DoC = depth / (double) length;
    
    double gene_rpkm = 1000000000 * (gene_avg_DoC / (double)total_map_depth);
    return new Expression(gene_rpkm, 0.0);
}

void GeneFactory::get_genes(const string& ref, 
                            const fsa::GFF_database& gffdb, 
                            GeneTable& genes,
                            const std::set<string>* gene_name_filter)
{
    GFF_database ref_gff_db = gffdb.chromosome_features(ref);
    get_genes(ref_gff_db, genes, gene_name_filter);
}

void GeneFactory::get_genes(const GFF_database& gffdb, 
                            GeneTable& genes,
                            const std::set<string>* gene_name_filter)
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
            
            string short_name;
            att_itr = gff_rec.attributes.find("Name");
            
            // No name record?  Use the ID as the short name, too.
            if (att_itr == gff_rec.attributes.end())
            {
                short_name = id;
            }
            else
            {
                short_name = att_itr->second.front();
            }
            
            Gene* gene = new Gene();
            gene->ID = id;
            gene->short_name = short_name;
            if (!gene_name_filter || 
                gene_name_filter->find(short_name) != gene_name_filter->end())
                genes[id] = gene;
        }
        
        //    GFF::AttributeTable::const_iterator att_itr;

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
                mRNA->ID = id;
                mRNAs[id] = mRNA;
            }
            else
            {
                //cerr << "No gene record "  << gene_id 
                //     << " for transcript " << id << endl;
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
            //const string& id = att_itr->second.front();
            
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
                //cerr << "No transcript record "  << mRNA_id 
                //<< " for exon " << id << endl;
                continue;
            }
            
        }
    }
    
    for (map<string,Transcript*>::iterator itr = mRNAs.begin(); 
         itr != mRNAs.end();
         ++itr)
    {
        Transcript& mRNA = *(itr->second);
        if (mRNA.exons.size() == 0)
        {
            Gene* parent = mRNA.gene();
            parent->transcripts.erase(mRNA.ID);
            // FIXME: leaking transcripts here.
        }
        else
        {
            sort(mRNA.exons.begin(), mRNA.exons.end(),exon_lt_by_pos); 
        }
    }
}

uint32_t total_exonic_depth(const GeneTable& genes,
                            const vector<short>& DoC,
                            const fsa::Sequence* ref_str)
{
    uint32_t total_depth = 0;
    for (GeneTable::const_iterator gene_itr = genes.begin();
         gene_itr != genes.end();
         ++gene_itr)
    {
        total_depth += gene_itr->second->exonic_depth(DoC, ref_str);
    }
    
    return total_depth;
}

void calculate_gene_expression(const GeneTable& genes,
                               const vector<short>& DoC,
                               const fsa::Sequence* ref_str,
                               uint32_t total_map_depth,
                               map<string, Expression*>& gene_expression)
{
    for (GeneTable::const_iterator itr = genes.begin();
         itr != genes.end();
         ++itr)
    {
        Expression* expr = itr->second->expression(DoC, total_map_depth, ref_str);
        if (expr != NULL)
            gene_expression[itr->second->short_name] = expr;
        
        //fprintf(expr_out, "%s\t%lf\t%lf\n", itr->first.c_str(), expr->rpkm, expr->mend);
    }
}

void print_gene_expression(FILE* expr_out, 
                           const map<string, Expression*>& gene_expression)
{
    for (map<string, Expression*>::const_iterator itr = gene_expression.begin();
         itr != gene_expression.end();
         ++itr)
    {
        Expression* expr = itr->second;
        fprintf(expr_out, "%s\t%lf\t%lf\n", itr->first.c_str(), expr->rpkm, expr->mend);
    }
}
