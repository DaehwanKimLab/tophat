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
  $Id:
 ==========================================================================*/

#ifndef SEQAN_HEADER_PSEUDOCOUNT_BASE_H
#define SEQAN_HEADER_PSEUDOCOUNT_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Class.Pseudocount:
..summary:Holds the pseudocounts for each residue of a given sequence alphabet.
..cat:Motif Search
..signature:Pseudocount<TValue, TSpec>
..param.TValue:The type of sequence which is considered.
...metafunction:Metafunction.Value
...type:Spec.Dna
...type:Spec.AminoAcid
..param.TSpec:Specialization tag for determining the pseudocount method.
...type:Spec.CMode 
...type:Spec.PMode
*/

template<typename TValue, typename TSpec>
class Pseudocount;

//////////////////////////////////////////////////////////////////////////////
//Metafunctions
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Pseudocount

template<typename TValue, typename TSpec>
struct Value< Pseudocount<TValue, TSpec> >
{
	typedef TValue Type;
};
template<typename TValue, typename TSpec>
struct Value< Pseudocount<TValue, TSpec> const>
{
	typedef TValue const Type;
};

/////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
