//
//  FastaTools.h
//  TopHat
//
//  Created by Harold Pimentel on 10/27/11.
//

#ifndef TopHat_FastaTools_h
#define TopHat_FastaTools_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
// #include <ifstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "common.h"

#define LINE_BUF_SIZE 1024
#define ID_BUF_SIZE 1024
#define DESC_BUF_SIZE 1024

struct FastaRecord {
    // The identifier after ">"
    std::string id_;
    
    // The description after the identifier
    std::string desc_;

    // The sequence
    std::string seq_;

};

class FastaReader {
public:
    FastaReader();    
    FastaReader(std::string fname);
    ~FastaReader();
    void init(std::string fname);
    bool good() const;
    void next(FastaRecord* rec);    
private:
    std::string fname_;
    std::string prev_line_;
    
    std::ifstream ifstream_;
    
    // We assume the description is no longer than 2048 characters long
    // TODO: Think of a better way of doing this in the future.
    char line_buf_[LINE_BUF_SIZE];    
    
    char id_buf_[ID_BUF_SIZE];
    
    char desc_buf_[DESC_BUF_SIZE];
    
    // variable to check if stream is primed (has already been initialized)
    bool isPrimed_;
};

class FastaWriter {
public:
    FastaWriter();
    FastaWriter(std::string fname);
    ~FastaWriter();
    void init(std::string fname);
    void write(FastaRecord* rec, size_t column_size = 60);    
private:
    std::string fname_;
    
    std::ofstream ofstream_;
    
    bool isPrimed_;
};
#endif
