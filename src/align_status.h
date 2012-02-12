#ifndef ALIGN_STATUS_H
#define ALIGN_STATUS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <vector>
#include <cassert>
#include <cstring>
#include <seqan/sequence.h>
#include <bam/sam.h>
#include "common.h"

#include "bwt_map.h"
#include "junctions.h"
#include "insertions.h"
#include "deletions.h"
#include "fusions.h"

using namespace std;

class Coverage;

/**
 */
struct AlignStatus
{
public:
  /**
   * Is there an alignment?
   */
  int _alignment_score;

public:
  AlignStatus();
  AlignStatus(const BowtieHit& bh,
	      const JunctionSet& gtf_junctions,
	      const JunctionSet& junctions,
	      const InsertionSet& insertions,
	      const DeletionSet& deletions,
	      const FusionSet& fusions,
	      const Coverage& coverage);

  bool operator<(const AlignStatus& rhs) const;
  bool operator==(const AlignStatus& rhs) const;
  bool operator!=(const AlignStatus& rhs) const;
};

#endif
