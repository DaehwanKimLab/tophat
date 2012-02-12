/*
 *  coverage.cpp
 *  TopHat
 *
 *  Created by Daehwan Kim on 2/11/2012 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
#define SVN_REVISION "XXX"
#endif

#include "coverage.h"

Coverage::Coverage()
{
  debug = false;
}

Coverage::~Coverage()
{
  clear();
}

void Coverage::add_coverage(RefID refid, int pos, int length)
{
  GenomeCoverage::iterator itr = genomeCoverage.find(refid);
  if (itr == genomeCoverage.end())
    {
      genomeCoverage[refid] = PosCoverage();
      itr = genomeCoverage.find(refid);
    }

  PosCoverage::iterator itr2 = get_contig(itr->second, pos);
  if (itr2 == itr->second.end())
    {
      itr->second[pos] = vector<int>(length + 1,  0);
      itr2 = itr->second.find(pos);

      vector<int>& contig_coverage = itr2->second;
      contig_coverage[0] = 1;
      contig_coverage[length] = -1;
    }
  else
    {
      // daehwan - remove this
      if (debug)
	{
	  fprintf(stderr, "found2\n");
	}
      
      const size_t first = pos - itr2->first;
      const size_t last = first + length;
      const size_t resize = last + 1;

      if (resize > itr2->second.size())
	itr2->second.resize(resize, 0);

      itr2->second[first] += 1;
      itr2->second[last] -= 1;
    }

  PosCoverage::iterator itr_lower = itr2; ++itr_lower;
  PosCoverage::iterator itr_upper = itr->second.upper_bound(itr2->first + itr2->second.size() - 1);

  PosCoverage::iterator itr_temp = itr_lower;
  while (itr_temp != itr_upper)
    {
      merge_contig(itr2->first, itr2->second, itr_temp->first, itr_temp->second);
      ++itr_temp;
    }

  if  (itr_lower != itr_upper)
    itr->second.erase(itr_lower, itr_upper);
}

void Coverage::merge_with(const Coverage& other)
{
  GenomeCoverage::const_iterator other_itr = other.genomeCoverage.begin();
  for (; other_itr != other.genomeCoverage.end(); ++other_itr)
    {
      GenomeCoverage::iterator itr = genomeCoverage.find(other_itr->first);
      if (itr == genomeCoverage.end())
	{
	  genomeCoverage[other_itr->first] = other_itr->second;
	  continue;
	}
      PosCoverage::const_iterator other_pos_itr = other_itr->second.begin();
      for (; other_pos_itr != other_itr->second.end(); ++other_pos_itr)
	{
	  PosCoverage::iterator insert_itr = get_contig(itr->second, other_pos_itr->first);
	  if (insert_itr == itr->second.end())
	    {
	      itr->second[other_pos_itr->first] = other_pos_itr->second;
	      insert_itr = itr->second.find(other_pos_itr->first);
	    }
	  else
	    {
	      merge_contig(insert_itr->first, insert_itr->second,
			   other_pos_itr->first, other_pos_itr->second);
	    }

	  int pos = insert_itr->first;
	  size_t length = insert_itr->second.size();
	  PosCoverage::iterator itr_lower = insert_itr; ++itr_lower;
	  PosCoverage::iterator itr_upper = itr->second.upper_bound(pos + length - 1);
	  
	  PosCoverage::iterator itr_temp = itr_lower;
	  while (itr_temp != itr_upper)
	    {
	      merge_contig(insert_itr->first, insert_itr->second,
			   itr_temp->first, itr_temp->second);
	      ++itr_temp;
	    }
	  
	  if  (itr_lower != itr_upper)
	    itr->second.erase(itr_lower, itr_upper);
	}
    }
}

void Coverage::merge_contig(int pos, vector<int>& cov, int pos2, const vector<int>& cov2)
{
  assert (pos <= pos2);
  size_t resize = cov2.size() + pos2 - pos;
  if (resize > cov.size())
      cov.resize(resize, 0);

  for (size_t i = 0; i < cov2.size(); ++i)
    cov[i + pos2 - pos] += cov2[i];
}

PosCoverage::iterator Coverage::get_contig(PosCoverage& posCoverage, int pos)
{
  PosCoverage::iterator itr_contig = posCoverage.lower_bound(pos + 1);
  if (itr_contig != posCoverage.begin())
    {
      --itr_contig;
      if (pos >= itr_contig->first && pos < itr_contig->first + (int)itr_contig->second.size())
	return itr_contig;
    }
  
  return posCoverage.end();
}

void Coverage::calculate_coverage()
{
  GenomeCoverage::iterator itr = genomeCoverage.begin();
  for (; itr != genomeCoverage.end(); ++itr)
    {
      PosCoverage::iterator itr2 = itr->second.begin();
      for (; itr2 != itr->second.end(); ++itr2)
	{
	  vector<int>& contig_coverage = itr2->second;
	  for (size_t i = 1; i < contig_coverage.size(); ++i)
	    contig_coverage[i] += contig_coverage[i-1];
	}
    }
}

int Coverage::get_coverage(RefID refid, int pos) const
{
  assert (pos >= 0);
  
  int coverage = 0;
  GenomeCoverage::const_iterator itr = genomeCoverage.find(refid);
  if (itr != genomeCoverage.end())
    {
      PosCoverage::const_iterator itr_contig = itr->second.lower_bound(pos + 1);
      if (itr_contig != itr->second.begin())
	{
	  --itr_contig;
	  const vector<int>& contig_coverage = itr_contig->second;
	  if (pos >= itr_contig->first)
	    {
	      int index = pos - itr_contig->first;
	      if (index < contig_coverage.size())
		coverage = contig_coverage[index];
	    }
	}
    }
  
  return coverage;
}

void Coverage::clear()
{
  genomeCoverage.clear();
}

void Coverage::print_info() const
{
  size_t total_bases = 0;
  GenomeCoverage::const_iterator itr = genomeCoverage.begin();
  for (; itr != genomeCoverage.end(); ++itr)
    {
      fprintf(stderr, "Reference %d\n", itr->first);

      size_t bases = 0;
      PosCoverage::const_iterator itr2 = itr->second.begin();
      for (; itr2 != itr->second.end(); ++itr2)
	bases += (itr2->second.size() - 1);
	
      fprintf(stderr, "# of islands: %d, # of bases covered: %d\n", itr->second.size(), bases);
      total_bases += bases;
    }

  fprintf(stderr, "# of total bases: %d\n", total_bases);
  
  itr = genomeCoverage.begin();
  for (; itr != genomeCoverage.end(); ++itr)
    {
      fprintf(stderr, "Reference %d\n", itr->first);
      print_info(itr->second);
    }
}

void Coverage::print_info(const PosCoverage& posCoverage, int begin, int end) const
{
  PosCoverage::const_iterator itr = posCoverage.begin();
  for (; itr != posCoverage.end(); ++itr)
    if (itr->first >= begin && itr->first < end)
      print_info(itr->first, itr->second);
}

void Coverage::print_info(int pos, const vector<int>& cov) const
{
  fprintf(stderr, "\tPos: %d, size: %u\n", pos, cov.size());
  for (size_t i = 0; i < cov.size(); ++i)
    fprintf(stderr, "\t\t%d (%d)\n", cov[i], pos + (int)i);
}
