/*
 *  genes.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 2/1/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include <stdint.h>
#include <map>
#include <string>
#include <set>
#include <vector>
#include "FSA/gff.h"
#include "FSA/sequence.h"

struct Gene;
struct Transcript;

struct Exon
{
    Exon() : _start(0), _end(0), _transcript(NULL) {}
    
    uint32_t start() const      { return _start;      }
    uint32_t end()   const      { return _end;       }
    Gene* gene();                

    Transcript* transcript()    { return _transcript; }
    
    uint32_t            _start;
    uint32_t            _end;
    Transcript*         _transcript;
};

struct Transcript
{
    Transcript() : _gene(NULL) {}
    std::vector<Exon*> exons;
    
    uint32_t start() const
    { 
        assert(exons.size()); 
        return exons.front()->start(); 
    }
    
    uint32_t end() const
    { 
        assert(exons.size()); 
        return exons.back()->end(); 
    }
    
    Gene* gene() const
    { 
        return _gene;    
    }
    
    Gene* _gene;
    std::string ID;
};

struct Expression
{
    Expression(double R = 0.0, double M = 0.0) : rpkm(R), mend(M) {}
    double rpkm;
    double mend;
};

struct Gene
{
    std::string ID;
    std::string short_name;
    //std::string ref_name;
    std::map<std::string, Transcript*> transcripts;
    
    Expression* expression(const std::vector<short>& DoC, 
                           uint32_t total_map_depth,
                           const fsa::Sequence* ref_str) const;
    
    uint32_t exonic_length(const fsa::Sequence* ref_str) const;
    uint32_t exonic_depth(const std::vector<short>& DoC, 
                          const fsa::Sequence* ref_str) const;
    
    std::pair<uint32_t, uint32_t> coords() const;
};

// This simple container should be keyed by ID and short names
typedef std::map<std::string, Gene*> GeneTable;

struct GeneFactory
{
    void get_genes(const fsa::GFF_database& gffdb, 
                   GeneTable& genes,
                   const std::set<std::string>* gene_name_filter = NULL);
    
    void get_genes(const std::string& ref, 
                   const fsa::GFF_database& gffdb, 
                   GeneTable& genes,
                   const std::set<std::string>* gene_name_filter = NULL);
};

uint32_t total_exonic_depth(const GeneTable& genes,
                            const std::vector<short>& DoC,
                            const fsa::Sequence* ref_str);

void print_gene_expression(FILE* expr_out, 
                           const std::map<std::string, Expression*>& gene_expression);

void calculate_gene_expression(const GeneTable& genes,
                               const std::vector<short>& DoC,
                               const fsa::Sequence* ref_seq,
                               uint32_t total_map_depth,
                               std::map<std::string, Expression*>& gene_expression);
