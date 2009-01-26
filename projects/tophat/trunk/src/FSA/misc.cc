
/**
 * \file misc.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include "misc.h"

using namespace fsa;


Regexp Util::re_chr_stripper = Regexp ("^(chr)?(.+)");


bool Util::exists_file (const std::string& filename) {

  struct stat fileinfo;
  bool gotattr;

  // Attempt to get the file attributes
  int statresults = stat (filename.c_str(), &fileinfo);

  if (statresults == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    gotattr = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    gotattr = false;

  }
  
  return gotattr; 

}

void Util::strip_leading_chr (std::string& chromosome) {

  if (re_chr_stripper.Match (chromosome.c_str())) // strip off leading 'chr' if present
    chromosome = std::string (re_chr_stripper[2]);

}
