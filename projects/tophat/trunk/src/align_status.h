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
#include "junctions.h"


using namespace std;


/**
 * The main purpose of this struct is to provide a
 * (fairly primitive) method for ranking competing alignments
 * for a read. In general,
 * continguous > spliced > in/dels
 * Ties are broken with the number of mismatches
 */
struct AlignStatus
{


private:
	/**
	 * Is this alignment free of indels. Note, if there is no alignment
	 * this should still be false.
	 */
	bool _indelFreeAlignment;

	/**
	 * Is this alignment free of unannotated splices? Note, if there is no alignment
	 * this should still be false
	 */
	bool _unannotatedSpliceFreeAlignment;

	/**
	 * Is there an alignment?
	 */
	bool _aligned;

public:
  AlignStatus();
  AlignStatus(const BowtieHit& bh, const JunctionSet& gtf_junctions);

  bool operator<(const AlignStatus& rhs) const;
  bool operator==(const AlignStatus& rhs) const;
  bool operator!=(const AlignStatus& rhs) const;
};

#endif
