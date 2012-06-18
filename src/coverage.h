#ifndef COVERAGE_H
#define COVERAGE_H
/*
 *  coverage.h
 *  TopHat
 *
 *  Created by Daehwan Kim on 2/11/2012
 */

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>

#include "common.h"
#include "bwt_map.h"
using namespace std;

typedef int8_t cov_t;
typedef std::map<int, vector<cov_t> > PosCoverage;
typedef std::map<RefID, PosCoverage> GenomeCoverage;

static const cov_t cov_max_value = std::numeric_limits<cov_t>::max();
static const cov_t cov_min_value = std::numeric_limits<cov_t>::min();

class Coverage
{
 public:
  Coverage();
  ~Coverage();

 public:
  // this should be called before calculate_coverage() call.
  void add_coverage(RefID refid, int pos, int length);
  void merge_with(const Coverage& other);

  // calculate real coverage instead of difference
  void calculate_coverage();

  // this can be called after calculate_coverage() call.
  int get_coverage(RefID refid, int pos) const;

  void clear();

  // for debug purposes;
  void print_info() const;
  void print_info(const PosCoverage& posCoverage, int begin = 0, int end = std::numeric_limits<int>::max()) const;
  void print_info(int pos, const vector<cov_t>& cov) const;
  
 private:
  void merge_contig(int pos, vector<cov_t>& cov, int pos2, const vector<cov_t>& cov2);
  PosCoverage::iterator get_contig(PosCoverage& posCoverage, int pos);
  
 private:
  GenomeCoverage genomeCoverage;

 public:
  bool debug;
};

#endif
