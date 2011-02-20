 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: graph_utility_parsing.h 1911 2008-05-02 09:28:04Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_UTILITY_PARSING_H
#define SEQAN_HEADER_GRAPH_UTILITY_PARSING_H


#include <fstream>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File reading
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TPath, typename TStringSet, typename TNames>
inline unsigned int
_loadSequences(TPath const& in_path,
			   TStringSet& origStrSet,
			   TNames& names)
{
	typedef typename Size<TStringSet>::Type TSize;
	
	// Count sequences and read names
	TSize seqCount = 0;
	std::ifstream file;
	std::stringstream input;
	input << in_path;
	file.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
	if (!file.is_open()) return 0;
	while (!_streamEOF(file)) {
		String<char> id;
		readID(file, id, Fasta());
		appendValue(names, id);
		goNext(file, Fasta());
		++seqCount;
	}

	// Load sequences
	file.clear();
	file.seekg(0, std::ios_base::beg);
	resize(origStrSet, seqCount);
	TSize count = 0;
	for(TSize i = 0; (i < seqCount) && !_streamEOF(file); ++i) 	{
		read(file, origStrSet[i], Fasta());
		count += length(origStrSet[i]);
	}
    file.close();

	return count;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
