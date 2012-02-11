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

  PosCoverage::iterator itr2 = itr->second.find(pos);
  if (itr2 == itr->second.end())
    {
      itr->second[pos] = vector<int>();
      itr2 = itr->second.find(pos);
    }

  vector<int>& contig_coverage = itr2->second;
  if (contig_coverage.size() < length)
    contig_coverage.resize(length, 0);

  contig_coverage[0] += 1;
  if (contig_coverage.size() > length + 1)
    contig_coverage[length] -= 1;

  PosCoverage::iterator itr_lower = itr->second.lower_bound(pos + 1);
  PosCoverage::iterator itr_upper = itr->second.upper_bound(pos + length - 1);

  PosCoverage::iterator itr_temp = itr_lower;
  while (itr_temp != itr_upper)
    {
      merge_contig(itr2, itr_temp);
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

    }
}

void Coverage::merge_contig(PosCoverage::iterator l, PosCoverage::const_iterator r)
{
  assert (l->first <= r->first);
  size_t resize = r->second.size() + r->first - l->first;
  if (resize > l->second.size())
      l->second.resize(resize, 0);

  for (size_t i = 0; i < r->second.size(); ++i)
    l->second[i + r->first - l->first] += r->second[i];
}

PosCoverage::iterator Coverage::get_contig(PosCoverage& posCoverage, int pos)
{
  PosCoverage::iterator itr_contig = posCoverage.upper_bound(pos);
  if (itr_contig != posCoverage.begin())
    {
      --itr_contig;
      if (pos >= itr_contig->first && pos < itr_contig->first + itr_contig->second.size())
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

int Coverage::get_coverage(RefID refid, int pos)
{
  assert (pos >= 0);
  
  int coverage = 0;
  GenomeCoverage::iterator itr = genomeCoverage.find(refid);
  if (itr != genomeCoverage.end())
    {
      PosCoverage::iterator itr_contig = get_contig(itr->second, pos);
      if (itr_contig != itr->second.end())
	{
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
  GenomeCoverage::const_iterator itr = genomeCoverage.begin();
  for (; itr != genomeCoverage.end(); ++itr)
    {
      fprintf(stderr, "Reference %d\n", itr->first);
      print_info(itr->second);
    }
}

void Coverage::print_info(const PosCoverage& posCoverage) const
{
  PosCoverage::const_iterator itr = posCoverage.begin();
  for (; itr != posCoverage.end(); ++itr)
    {
      fprintf(stderr, "\tPos: %d\n", itr->first);
      const vector<int>& contig_coverage = itr->second;
      for (size_t i = 0; i < contig_coverage.size(); ++i)
	fprintf(stderr, "\t\t%d (%d)\n", contig_coverage[i], itr->first + i);
    }
}
