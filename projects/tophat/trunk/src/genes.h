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
#include <vector>
#include "FSA/gff.h"

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
    std::string ref_name;
    std::map<std::string, Transcript*> transcripts;
    
    Expression* expression(const std::vector<short>& DoC, 
                           uint32_t total_map_depth) const;
    
    uint32_t exonic_length() const;
    uint32_t exonic_depth(const std::vector<short>& DoC) const;
    
    std::pair<uint32_t, uint32_t> coords() const;
};

// This simple container should be keyed by ID and short names
typedef std::map<std::string, Gene*> GeneTable;

class GeneFactory
{
    void get_genes(const fsa::GFF_database& gffdb, GeneTable& genes);
};
