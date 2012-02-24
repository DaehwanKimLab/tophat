// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_ESA_DRAWING_H
#define SEQAN_HEADER_INDEX_ESA_DRAWING_H

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TFile, typename TText, typename TESASpec>
void write(TFile & file, 
	   Index<TText, IndexEsa<TESASpec> > & stree,
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	typedef Index<TText, IndexEsa<TESASpec> > TIndex;
	
	_streamWrite(file, "digraph G {\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Graph Attributes */\n");
	_streamWrite(file, "graph [rankdir = LR];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Node Attributes */\n");
	_streamWrite(file, "node [shape = ellipse, fillcolor = lightgrey, style = filled, fontname = \"Times-Italic\"];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Edge Attributes */\n");
	_streamWrite(file, "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n");
	_streamPut(file, '\n');

	_streamWrite(file, "/* Edges */\n");
	typedef typename Iterator<TIndex, TopDown<ParentLinks<Preorder> > >::Type TIterator;
	typedef typename Iterator<TIndex, TopDown<> >::Type TIteratorSimple;
	TIterator it(stree);

	for(;!atEnd(it);++it) 
	{
		// dump node
       		_streamWrite(file, "\"[");
 		_streamPutInt(file, value(it).range.i1);
		_streamPut(file, ':');
		_streamPutInt(file, value(it).range.i2);
       		_streamWrite(file, ")\"");
       		if (!isRightTerminal(it))
			_streamWrite(file, " [style = dashed]");
       		_streamWrite(file, ";\n");

		// dump edge from parent (if not root)
		if (!isRoot(it)) {
			TIteratorSimple src(container(it), nodeUp(it));

			_streamWrite(file, "\"[");
			_streamPutInt(file, value(src).range.i1);
			_streamPut(file, ':');
			_streamPutInt(file, value(src).range.i2);
			_streamWrite(file, ")\"");

			_streamWrite(file, " -> ");

			_streamWrite(file, "\"[");
			_streamPutInt(file, value(it).range.i1);
			_streamPut(file, ':');
			_streamPutInt(file, value(it).range.i2);
			_streamWrite(file, ")\"");

			_streamWrite(file, " [label = \"");
			_streamWrite(file, parentEdgeLabel(it));
			_streamWrite(file, "\"];\n");
		}
	}
	_streamPut(file, '\n');

	_streamWrite(file, "}\n");
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
