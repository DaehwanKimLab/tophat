/*
 *  smith-waterman.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 2/5/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */
#include <seqan/sequence.h>
#include <string>
//typedef seqan::Dna5 dna_char;
//typedef seqan::Dna5String dna_string;

typedef char dna_char;
typedef std::string dna_string;

void align(const std::string & S, 
           const std::string & T,
           int match_score,
           int mismatch_score,
           int gap_open,
           int gap_extension);