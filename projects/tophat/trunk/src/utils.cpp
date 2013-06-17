/*
 *  utils.cpp
 *  TopHat
 *
 *  Created by Daehwan Kim on 12/28/11.
 *  Copyright 2011 Daehwan Kim. All rights reserved.
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "utils.h"

bool calculate_offsets(const vector<string>& fnames,
		       vector<uint64_t>& ids,
		       vector<vector<int64_t> >& offsets)
{
  vector<vector<pair<uint64_t, int64_t> > > index_list(fnames.size());
  for (size_t i = 0; i < fnames.size(); ++i)
    {
      const string& fname = fnames[i];
      
      vector<string> temp_fnames;
      if (fname.substr(fname.length() - 4) == ".bam")
	{
	  temp_fnames.push_back(fname);
	}
      else
	{
	  size_t j = 0;
	  while (true)
	    {
	      char suffix[128];
	      sprintf(suffix, "%lu.bam", j);
	      string temp_fname = fname + suffix;
	      temp_fnames.push_back(temp_fname);

	      string temp_index_fname = temp_fname + ".index";
	      ifstream index_file(temp_index_fname.c_str());
	      if (!index_file.is_open())
		break;

	      ++j;
	    }
	}
      
      for (size_t j = 0; j < temp_fnames.size(); ++j)
	{
	  ifstream reads_index_file((temp_fnames[j] + ".index").c_str());
	  if (!reads_index_file.is_open())
	    continue;

	  string line;
	  while (getline(reads_index_file, line))
	    {
	      istringstream istream(line);
	      uint64_t read_id;
	      int64_t offset;
	      
	      istream >> read_id >> offset;
	      index_list[i].push_back(make_pair(read_id, offset));
	    }
	}
    }

  // too small for blocking
  for (size_t i = 0; i < index_list.size(); ++i)
    {
      if (index_list[i].size() < (size_t)num_threads)
	return false;
    }

  offsets.resize(num_threads - 1);
  for (int i = 1; i < num_threads; ++i)
    {
      size_t index = index_list.back().size() / num_threads * i;
      const pair<uint64_t, int64_t>& data = index_list.back()[index];
      uint64_t id = data.first;
      int64_t offset = data.second;

      ids.push_back(id);
      offsets[i-1].push_back(offset);

      // fprintf(stderr, "%lu)read %lu - offset %ld\n", index_list.size() - 1, id, offset);
  
      for (int j = index_list.size() - 2; j >= 0; --j)
	{
	  size_t other_index = index_list[j].size() / num_threads * i;
	  uint64_t other_id = index_list[j][other_index].first;

	  while (other_id > id && other_index > 0)
	    {
	      other_index -= 1;
	      other_id = index_list[j][other_index].first;
	    }

	  while (other_index + 1 < index_list[j].size() && index_list[j][other_index+1].first < id)
	    {
	      other_index += 1;
	      other_id = index_list[j][other_index].first;
	    }

	  int64_t other_offset = index_list[j][other_index].second;
	  
	  if (other_id > id)
	    {
	      other_id = 0;
	      other_offset = 0;
	    }

	  id = other_id;
	  offsets[i-1].push_back(other_offset);

	  // fprintf(stderr, "\t%d)read %lu - offset %lu\n", j, other_id, other_offset);
	}

      reverse(offsets[i-1].begin(), offsets[i-1].end());
    }
    
  return true;
}

void calculate_offsets_from_ids(const string& fname,
				const vector<uint64_t>& ids,
				vector<int64_t>& offsets)
{
  vector<pair<uint64_t, int64_t> > index_list;
  ifstream reads_index_file((fname + ".index").c_str());

  string line;
  uint64_t last_id = 0;
  int64_t last_offset = 0;
  for (size_t i = 0; i < ids.size(); ++i)
    {
      uint64_t ref_id = ids[i];
      while (getline(reads_index_file, line))
	{
	  istringstream istream(line);
	  uint64_t read_id;
	  int64_t offset;
	  
	  istream >> read_id >> offset;
	  if (read_id > ref_id)
	    {
	      assert(last_id <= ref_id);
	      offsets.push_back(last_offset);

	      // fprintf(stderr, "ref read: %lu\n", ref_id);
	      // fprintf(stderr, "\tread %lu - offset %lu\n", last_id, last_offset);
	    }

	  last_id = read_id;
	  last_offset = offset;

	  if (last_id > ref_id)
	    break;
	}
    }

  if (ids.size() != offsets.size())
    offsets.clear();
}
