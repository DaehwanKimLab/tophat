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
  $Id: graph_align_tcoffee_guidetree.h 1813 2008-03-31 16:06:14Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_GUIDETREE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_GUIDETREE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Guide Tree
//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
// Neighbor Joining
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


/**
.Function.slowNjTree:
..summary:Computes a guide tree from a distance matrix.
..cat:Graph
..signature:
njTree(mat, graph)
..param.mat:A string of pairwise distance values, representing a square matrix.
...type:Class.String
...remarks: String must use double values because the algorithm recycles this string to store intermediate possibly fractional results.
..param.graph:Out-parameter:The guide tree.
...type:Spec.Tree
..returns:void
*/
template<typename TValue, typename TStringSpec, typename TCargo, typename TSpec>
inline void
njTree(String<TValue, TStringSpec>& mat, 
	   Graph<Tree<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	typedef String<TValue, TStringSpec> TMatrix;
	typedef typename Size<TMatrix>::Type TSize;
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	TSize nseq = (TSize) std::sqrt((double)length(mat));

	//for(TSize i=0;i<nseq;++i) {
	//	for(TSize j=0;j<nseq;++j) {
	//		std::cout << getValue(mat, i*nseq+j) << ",";
	//	}
	//	std::cout << std::endl;
	//}

	// First initialization
	clearVertices(g);
	if (nseq == 1) {
		g.data_root = addVertex(g);
		return;
	} else if (nseq == 2) {
		TVertexDescriptor v1 = addVertex(g);
		TVertexDescriptor v2 = addVertex(g);
		TVertexDescriptor internalVertex = addVertex(g);
		addEdge(g, internalVertex, v1, (TCargo) mat[1] / 2);
		addEdge(g, internalVertex, v2, (TCargo) mat[1] / 2);
		g.data_root = internalVertex;
		return;
	}
	String<TCargo> av;    // Average branch length to a combined node
	fill(av,nseq,0);

	String<TVertexDescriptor> connector;   // Nodes that need to be connected
	resize(connector, nseq);

	for(TSize i=0;i<nseq;++i) {
		addVertex(g);  // Add all the nodes that correspond to sequences
		connector[i] = i;
		mat[i*nseq+i] = 0;
	}

	// Main cycle
	TCargo fnseqs=(TCargo) nseq;
	for(TSize nc=0; nc<(nseq-3); ++nc) {
		TCargo sumOfBranches = 0;

		// Determine the sum of all branches and
		// copy upper triangle matrix to lower triangle
		for(TSize col=1; col<nseq; ++col) 
			for(TSize row=0; row<col; ++row) 
				sumOfBranches += mat[col*nseq+row] = mat[row*nseq+col];

		// Compute the sum of branch lengths for all possible pairs
		bool notFound = true;
		TCargo tmin = 0;	
		TSize mini = 0;  // Next pair of seq i and j to join
		TSize minj = 0;
		TCargo diToAllOthers = 0;
		TCargo djToAllOthers = 0;
		TCargo total = 0;
		for(TSize col=1; col<nseq; ++col)  {
			if (connector[col] != nilVertex) {
				for(TSize row=0; row<col; ++row) {
					if (connector[row] != nilVertex) {
						diToAllOthers = 0;
						djToAllOthers = 0;
						
						for(TSize i=0; i<nseq; ++i) {
							diToAllOthers += mat[i*nseq+row];
							djToAllOthers += mat[i*nseq+col];
						}

						total = diToAllOthers + djToAllOthers + (fnseqs - 2) * mat[row*nseq+col] + 2 * (sumOfBranches - diToAllOthers - djToAllOthers);
						total /= (2*(fnseqs - 2));

						if ((notFound) || (total < tmin)) {
							notFound = false;
							tmin = total;
							mini = row;
							minj = col;
						}
					}
				}
			}
		}

		// Print nodes that are about to be joined
		//std::cout << mini << std::endl;
		//std::cout << minj << std::endl;
		//std::cout << tmin << std::endl;
		//std::cout << std::endl;
		
		// Compute branch lengths
		TCargo dMinIToOthers = 0;
		TCargo dMinJToOthers = 0;
		for(TSize i=0; i<nseq; ++i) {
			dMinIToOthers += mat[i*nseq + mini];
			dMinJToOthers += mat[i*nseq + minj];
		}
		TCargo dmin = mat[mini*nseq + minj];
		dMinIToOthers = dMinIToOthers / (fnseqs - 2);
		dMinJToOthers = dMinJToOthers / (fnseqs - 2);
		TCargo iBranch = (dmin + dMinIToOthers - dMinJToOthers) / 2;
		TCargo jBranch = dmin - iBranch;
		iBranch -= av[mini];
		jBranch -= av[minj];
		
		// Set negative branch length to zero
		if( iBranch < 0) iBranch = 0;
		if( jBranch < 0) jBranch = 0;
	
		// Print branch lengths
		//std::cout << iBranch << std::endl;
		//std::cout << jBranch << std::endl;
		//std::cout << std::endl;
		
		// Build tree
		TVertexDescriptor internalVertex = addVertex(g);
		addEdge(g, internalVertex, connector[mini], (TCargo) iBranch);
		addEdge(g, internalVertex, connector[minj], (TCargo) jBranch);

		// Remember the average branch length for the new combined node
		// Must be subtracted from all branches that include this node
		if(dmin < 0) dmin = 0;
		av[mini] = dmin / 2;


		// Re-initialisation
		// mini becomes the new combined node, minj is killed
		--fnseqs;
		connector[minj] = nilVertex;
		connector[mini] = internalVertex;

		for(TSize j=0; j<nseq; ++j) {
			if( connector[j] != nilVertex ) {
				// Use upper triangle
				if((TSize) mini < j) mat[mini*nseq+j] = (TValue) ((mat[mini*nseq+j] + mat[minj*nseq+j]) / 2);
				if((TSize) mini > j) mat[j*nseq+mini] = (TValue) ((mat[mini*nseq+j] + mat[minj*nseq+j]) / 2);
			}
		}
		for(TSize j=0; j<nseq; ++j)
			mat[j*nseq+minj] = mat[minj*nseq+j] = 0;
	}

	// Only three nodes left

	// Find the remaining nodes
	String<TSize> l;
	resize(l,3);
	TSize count = 0;
	for(TSize i=0; i<nseq; ++i) {
		if(connector[i] != nilVertex) {
			l[count] = i;
			++count;
		}
	}

	// Remaining nodes
	//std::cout << l[0] << std::endl;
	//std::cout << l[1] << std::endl;
	//std::cout << l[2] << std::endl;
	//std::cout << std::endl;

	String<TCargo> branch;
	resize(branch, 3);
	branch[0] = (mat[l[0]*nseq+l[1]] + mat[l[0]*nseq+l[2]] - mat[l[1]*nseq+l[2]]) / 2;
	branch[1] = (mat[l[1]*nseq+l[2]] + mat[l[0]*nseq+l[1]] - mat[l[0]*nseq+l[2]]) / 2;
	branch[2] = (mat[l[1]*nseq+l[2]] + mat[l[0]*nseq+l[2]] - mat[l[0]*nseq+l[1]]) / 2;
    
	branch[0] -= av[l[0]];
	branch[1] -= av[l[1]];
	branch[2] -= av[l[2]];

	// Print branch lengths
	//std::cout << branch[0] << std::endl;
	//std::cout << branch[1] << std::endl;
	//std::cout << branch[2] << std::endl;
	//std::cout << std::endl;
    
	// Reset negative branch lengths to zero
	if( branch[0] < 0) branch[0] = 0;
	if( branch[1] < 0) branch[1] = 0;
	if( branch[2] < 0) branch[2] = 0;
    
	// Build tree
	TVertexDescriptor internalVertex = addVertex(g);
	addEdge(g, internalVertex, getValue(connector, l[0]), (TCargo) branch[0]);
	addEdge(g, internalVertex, getValue(connector, l[1]), (TCargo) branch[1]);
	TVertexDescriptor the_root = addVertex(g);
	addEdge(g, the_root, getValue(connector, l[2]), (TCargo) branch[2] / 2);
	addEdge(g, the_root, internalVertex, (TCargo) branch[2] / 2);
	g.data_root = the_root;
}






//////////////////////////////////////////////////////////////////////////////
// Unweighted Pair Group Mean Average (UPGMA)
//////////////////////////////////////////////////////////////////////////////


/**
.Tag.Upgma Configurator:
..summary:A tag to configure the guide tree construction.
*/

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Upgma Configurator.value.UpgmaMin:
	Uses the min operation in the upgma algorithm
*/

struct UpgmaMin_;
typedef Tag<UpgmaMin_> const UpgmaMin;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Upgma Configurator.value.UpgmaMax:
	Uses the max operation in the upgma algorithm
*/

struct UpgmaMax_;
typedef Tag<UpgmaMax_> const UpgmaMax;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Upgma Configurator.value.UpgmaAvg:
	Uses the average operation in the upgma algorithm
*/

struct UpgmaAvg_;
typedef Tag<UpgmaAvg_> const UpgmaAvg;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Upgma Configurator.value.UpgmaWeightAvg:
	Uses the weighted average operation in the upgma algorithm
*/

struct UpgmaWeightAvg_;
typedef Tag<UpgmaWeightAvg_> const UpgmaWeightAvg;



//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix, typename TActive, typename TSize>
inline void
_upgmaTreeMerge(TMatrix& mat, 
				TActive& active,
				TSize index_i,
				TSize index_j,
				TSize nseq,
				UpgmaWeightAvg) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TMatrix>::Type TValue;
	// Average
	for(TSize i=0;i<nseq;++i) {
		if ((i != index_i) && (i != index_j) && (active[i] != 0)) {
			if (index_i < i) {
				mat[index_i*nseq + i] = ((TValue) active[index_i] / (TValue) (active[index_i] + active[index_j])) * mat[index_i * nseq + i];
				if (index_j < i) mat[index_i*nseq + i] += ((TValue) active[index_j] / (TValue) (active[index_i] + active[index_j])) * mat[index_j * nseq + i];
				else mat[index_i*nseq + i] += ((TValue) active[index_j] / (TValue) (active[index_i] + active[index_j])) * mat[i * nseq + index_j];
			} else {
				mat[i*nseq + index_i] = ((TValue) active[index_i] / (TValue) (active[index_i] + active[index_j])) * mat[i * nseq + index_i];
				if (index_j < i) value(mat, i*nseq + index_i) += ((TValue) active[index_j] / (TValue) (active[index_i] + active[index_j])) * mat[index_j * nseq + i];
				else mat[i*nseq + index_i] += ((TValue) active[index_j] / (TValue) (active[index_i] + active[index_j])) * mat[i * nseq + index_j];
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////


template<typename TMatrix, typename TActive, typename TSize>
inline void
_upgmaTreeMerge(TMatrix& mat, 
				TActive& active,
				TSize index_i,
				TSize index_j,
				TSize nseq,
				UpgmaAvg) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TMatrix>::Type TValue;

	// Minimum
	for(TSize i=0;i<nseq;++i) {
		if ((i != index_i) && (i != index_j) && (active[i] != 0)) {
			TValue val1 = (index_i < i) ? mat[index_i * nseq + i] : mat[i * nseq + index_i];
			TValue val2 = (index_j < i) ? mat[index_j * nseq + i] : mat[i * nseq + index_j];
			if (index_i < i) mat[index_i * nseq + i] = (val1 + val2) / 2;
			else mat[i * nseq + index_i] = (val1 + val2) / 2;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////


template<typename TMatrix, typename TActive, typename TSize>
inline void
_upgmaTreeMerge(TMatrix& mat, 
				TActive& active,
				TSize index_i,
				TSize index_j,
				TSize nseq,
				UpgmaMin) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TMatrix>::Type TValue;

	// Minimum
	for(TSize i=0;i<nseq;++i) {
		if ((i != index_i) && (i != index_j) && (active[i] != 0)) {
			TValue newDist = (index_i < i) ? mat[index_i * nseq + i] : mat[i * nseq + index_i];
			TValue newDist2 = (index_j < i) ? mat[index_j * nseq + i] : mat[i * nseq + index_j];
			if (index_i < i) mat[index_i * nseq + i] = _min(newDist, newDist2);
			else mat[i * nseq + index_i] = _min(newDist, newDist2);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix, typename TActive, typename TSize>
inline void
_upgmaTreeMerge(TMatrix& mat, 
				TActive& active,
				TSize index_i,
				TSize index_j,
				TSize nseq,
				UpgmaMax) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TMatrix>::Type TValue;

	// Maximum
	for(TSize i=0;i<nseq;++i) {
		if ((i != index_i) && (i != index_j) && (active[i] != 0)) {
			TValue newDist = (index_i < i) ? mat[index_i * nseq + i] : mat[i * nseq + index_i];
			TValue newDist2 = (index_j < i) ? mat[index_j * nseq + i] : mat[i * nseq + index_j];
			if (index_i < i) mat[index_i * nseq + i] = _max(newDist, newDist2);
			else mat[i * nseq + index_i] = _max(newDist, newDist2);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TActive, typename TEdgeDescriptor>
inline void
_upgmaTreeMerge(Graph<Undirected<TCargo, TSpec> >& pairGraph, 
				TActive const&,
				TEdgeDescriptor best,
				UpgmaMax) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	TVertexDescriptor s = sourceVertex(pairGraph,best);
	TVertexDescriptor t = targetVertex(pairGraph,best);
	typedef String<TEdgeDescriptor> TEdgeString;
	typedef typename Iterator<TEdgeString>::Type TEdgeIter;
	TEdgeString removeEdges;
	for(TOutEdgeIterator outIt(pairGraph, s);!atEnd(outIt);goNext(outIt)) {
		if (targetVertex(outIt) == t) continue;
		TEdgeDescriptor e = findEdge(pairGraph, targetVertex(outIt), t);
		if (e == 0) appendValue(removeEdges, value(outIt));
		else {
			if (cargo(e) > cargo(value(outIt))) cargo(value(outIt)) = cargo(e);
		}
	}
	TEdgeIter eIt = begin(removeEdges);
	TEdgeIter eItEnd = end(removeEdges);
	for(;eIt != eItEnd; goNext(eIt)) removeEdge(pairGraph, value(eIt));
	removeVertex(pairGraph, t);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TActive, typename TEdgeDescriptor>
inline void
_upgmaTreeMerge(Graph<Undirected<TCargo, TSpec> >& pairGraph, 
				TActive const&,
				TEdgeDescriptor best,
				UpgmaMin) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	TVertexDescriptor s = sourceVertex(pairGraph,best);
	TVertexDescriptor t = targetVertex(pairGraph,best);
	for(TOutEdgeIterator outIt(pairGraph, s);!atEnd(outIt);goNext(outIt)) {
		if (targetVertex(outIt) == t) continue;
		TEdgeDescriptor e = findEdge(pairGraph, targetVertex(outIt), t);
		if (e != 0) {
			if (cargo(e) < cargo(value(outIt))) cargo(value(outIt)) = cargo(e);
		}
	}
	for(TOutEdgeIterator outIt(pairGraph, t);!atEnd(outIt);goNext(outIt)) {
		if (targetVertex(outIt) == s) continue;
		TEdgeDescriptor e = findEdge(pairGraph, targetVertex(outIt), s);
		if (e == 0) addEdge(pairGraph, s, targetVertex(outIt), cargo(value(outIt)));
	}
	removeVertex(pairGraph, t);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TActive, typename TEdgeDescriptor>
inline void
_upgmaTreeMerge(Graph<Undirected<TCargo, TSpec> >& pairGraph, 
				TActive const& active,
				TEdgeDescriptor best,
				UpgmaWeightAvg) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	TCargo infCargo = _getInfinity<TCargo>();
	TVertexDescriptor s = sourceVertex(pairGraph,best);
	TVertexDescriptor t = targetVertex(pairGraph,best);
	
	for(TOutEdgeIterator outIt(pairGraph, s);!atEnd(outIt);goNext(outIt)) {
		if (targetVertex(outIt) == t) continue;
		TEdgeDescriptor e1 = value(outIt);
		TEdgeDescriptor e2 = findEdge(pairGraph, targetVertex(outIt), t);
		if (e2 != 0) cargo(e1) = ((TCargo) value(active,s) / (TCargo) (value(active,s) + value(active,t))) * cargo(e1) + ((TCargo) value(active,t) / (TCargo) (value(active,s) + value(active,t))) * cargo(e2);
		else cargo(e1) = ((TCargo) value(active,s) / (TCargo) (value(active,s) + value(active,t))) * cargo(e1) + ((TCargo) value(active,t) / (TCargo) (value(active,s) + value(active,t))) * infCargo;
	}
	for(TOutEdgeIterator outIt(pairGraph, t);!atEnd(outIt);goNext(outIt)) {
		if (targetVertex(outIt) == s) continue;
		TEdgeDescriptor e = findEdge(pairGraph, targetVertex(outIt), s);
		TCargo c = ((TCargo) value(active,s) / (TCargo) (value(active,s) + value(active,t))) * infCargo + ((TCargo) value(active,t) / (TCargo) (value(active,s) + value(active,t))) * cargo(value(outIt));
		if (e == 0)  addEdge(pairGraph, s, targetVertex(outIt), c);
	}
	removeVertex(pairGraph, t);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TActive, typename TEdgeDescriptor>
inline void
_upgmaTreeMerge(Graph<Undirected<TCargo, TSpec> >& pairGraph, 
				TActive const&,
				TEdgeDescriptor best,
				UpgmaAvg) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	TCargo infCargo = _getInfinity<TCargo>();
	TVertexDescriptor s = sourceVertex(pairGraph,best);
	TVertexDescriptor t = targetVertex(pairGraph,best);
	
	for(TOutEdgeIterator outIt(pairGraph, s);!atEnd(outIt);goNext(outIt)) {
		if (targetVertex(outIt) == t) continue;
		TEdgeDescriptor e1 = value(outIt);
		TEdgeDescriptor e2 = findEdge(pairGraph, targetVertex(outIt), t);
		cargo(e1) = (e2 != 0) ? (cargo(e1) + cargo(e2)) / 2 : (cargo(e1) + infCargo) / 2;
	}
	for(TOutEdgeIterator outIt(pairGraph, t);!atEnd(outIt);goNext(outIt)) {
		if (targetVertex(outIt) == s) continue;
		TEdgeDescriptor e = findEdge(pairGraph, targetVertex(outIt), s);
		if (e == 0) {
			TCargo c = (infCargo + cargo(value(outIt))) / 2;
			addEdge(pairGraph, s, targetVertex(outIt), c);
		}
	}
	removeVertex(pairGraph, t);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringValue, typename TStringSpec, typename TCargo, typename TSpec, typename TTag>
inline void
upgmaTree(String<TStringValue, TStringSpec>& mat, 
		  Graph<Tree<TCargo, TSpec> >& g,
		  TTag) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	typedef String<TStringValue, TStringSpec> TMatrix;
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename Iterator<TMatrix>::Type TMatrixIter;
	
	
	// First initialization
	TSize nseq = (TSize) std::sqrt((double)length(mat));
	clearVertices(g);

	// Is it possible to make a guide tree?
	if (nseq == 1) {
		g.data_root = addVertex(g);
		return;
	} else if (nseq == 2) {
		TVertexDescriptor v1 = addVertex(g);
		TVertexDescriptor v2 = addVertex(g);
		TVertexDescriptor internalVertex = addVertex(g);
		addEdge(g, internalVertex, v1, (TCargo) 1);
		addEdge(g, internalVertex, v2, (TCargo) 1);
		g.data_root = internalVertex;
		return;
	}

	// Which entries in the matrix are still active and how many members belong to this group
	String<TSize> active;
	fill(active, nseq, 1);
	// Vertex descriptor that represents that entry
	String<TVertexDescriptor> nodes;
	reserve(nodes, nseq);


	// Find the minimal value
	bool notFound = true;
	TValue minVal = 0;
	TSize index_i = 0;
	TSize index_j = 1;
	for(TSize row=0;row<nseq;++row) {
		for(TSize col=row+1;col<nseq;++col) {
			if ((notFound) || (minVal > mat[row*nseq + col])) {
				notFound = false;
				minVal = mat[row*nseq + col];
				index_i = row;
				index_j = col;
			}
		}
		appendValue(nodes, addVertex(g));	// For each sequence one vertex
	}

	// Property map for sum of weights for each node
	String<TCargo> weights;
	fill(weights, nseq, (TCargo) 0);
	reserve(weights, 2*nseq - 1);

	// Merge groups
	TSize m = nseq;
	while (m>1) {
		// Merge nodes
		TVertexDescriptor internalNode = addVertex(g);

		//// Debug code
		//for(TSize i=0;i<nseq;++i) {
		//	if (value(active,i)==0) continue;
		//	for(TSize j=i+1;j<nseq;++j) {
		//		if (value(active,j)==0) continue;
		//		std::cout << value(mat, i*nseq+j) << ",";
		//	}
		//	std::cout << std::endl;
		//}
		//std::cout << minVal << ',' << index_i << ',' << index_j << ',' << std::endl;
		//std::cout << nodes[index_i] << ',' << nodes[index_j] << std::endl;
		//std::cout << std::endl;

		TCargo w = (TCargo) (minVal / 2);
		addEdge(g, internalNode, nodes[index_i], w - property(weights, nodes[index_i]));
		addEdge(g, internalNode, nodes[index_j], w - property(weights, nodes[index_j]));
		appendValue(weights, w);		

		// Get the new distance values
		_upgmaTreeMerge(mat, active, index_i, index_j, nseq, TTag());

		// Inactivate one group, adjust the member count for the other one
		active[index_i] += active[index_j];
		active[index_j] = 0;
		nodes[index_i] = internalNode;
		
		// Find new minimum
		notFound = true;
		for(TSize i=0;i<nseq;++i) {
			if (active[i] == 0) continue;
			for(TSize j=i+1;j<nseq;++j) {
				if (active[j] == 0) continue;
				if ((notFound) || (minVal > mat[i*nseq + j])) {
					notFound = false;
					minVal = mat[i*nseq + j];
					index_i = i;
					index_j = j;
				}
			}
		}
		--m;
	}
	g.data_root = numVertices(g) - 1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec1, typename TCargo, typename TSpec2, typename TTag>
inline void
upgmaTree(Graph<Undirected<TValue, TSpec1> >& pairGraph, 
		  Graph<Tree<TCargo, TSpec2> >& g,
		  TTag) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TValue, TSpec1> > TPairGraph;
	typedef Graph<Tree<TCargo, TSpec2> > TGuideTree;
	typedef typename VertexDescriptor<TGuideTree>::Type TVertexDescriptor;
	typedef typename VertexDescriptor<TPairGraph>::Type TVD;
	typedef typename EdgeDescriptor<TPairGraph>::Type TED;
	typedef typename Iterator<TPairGraph, EdgeIterator>::Type TEdgeI;
	typedef typename Iterator<TPairGraph, OutEdgeIterator>::Type TEdgeOutI;
	typedef typename Iterator<TPairGraph, VertexIterator>::Type TVertexI;
	typedef typename Size<TPairGraph>::Type TSize;

	// First initialization
	TCargo const maxVal = supremumValue<TCargo>();
	TSize nseq = numVertices(pairGraph);
	TCargo infCargo = _getInfinity<TCargo>();
	clearVertices(g);

	// Is it possible to make a guide tree?
	if (nseq == 1) {
		g.data_root = addVertex(g);
		return;
	} else if (nseq == 2) {
		TVertexDescriptor v1 = addVertex(g);
		TVertexDescriptor v2 = addVertex(g);
		TVertexDescriptor internalVertex = addVertex(g);
		addEdge(g, internalVertex, v1, (TCargo) 1);
		addEdge(g, internalVertex, v2, (TCargo) 1);
		g.data_root = internalVertex;
		return;
	}

	// Which entries in the matrix are still active and how many members belong to this group
	String<TSize> active;
	fill(active, nseq, 1);
	// Vertex descriptor that represents that entry
	typedef String<TVertexDescriptor> TNodeString;
	typedef typename Iterator<TNodeString, Standard>::Type TNodeIter;
	TNodeString nodes;
	resize(nodes, nseq);
	TNodeIter nodeIt = begin(nodes, Standard() );
	TNodeIter nodeItEnd = end(nodes, Standard() );
	for(;nodeIt<nodeItEnd;goNext(nodeIt)) 
		*nodeIt = addVertex(g);	// For each sequence one vertex

	// Find the minimal value for all vertices (with respect to all greater vertices)
	typedef Pair<TValue, TVD> TWeightEdgePair;
	typedef String<TWeightEdgePair> TMinValues;
	TMinValues minValues;
	fill(minValues, nseq, TWeightEdgePair(maxVal, 0));
	TEdgeI itE(pairGraph);
	for(;!atEnd(itE);goNext(itE)) {
		TVD s = sourceVertex(itE);
		TVD t = targetVertex(itE);
		if (cargo(*itE) < (minValues[s].i1)) 
			minValues[s] = TWeightEdgePair(cargo(*itE), t);		
	}
	// Find the overall minimum
	typedef typename Iterator<TMinValues, Standard>::Type TMinIter;
	TMinIter itMin = begin(minValues, Standard() );
	TMinIter itMinEnd = end(minValues, Standard() );
	TValue minVal = maxVal;
	TVD sourceBest = 0;
	TVD targetBest = 0;
	for(TVD index=0;itMin != itMinEnd; goNext(itMin), ++index) {
		if ((*itMin).i1 < minVal) {
			minVal = (*itMin).i1;
			sourceBest = index;
			targetBest = (*itMin).i2;
		}
	}
	TED best = 0;
	if (sourceBest == targetBest) { // If none is found we have to insert a new edge
		TVertexI itV(pairGraph);
		sourceBest = value(itV); goNext(itV);
		targetBest = value(itV);
		best = addEdge(pairGraph, sourceBest, targetBest, infCargo);
	} else best = findEdge(pairGraph, sourceBest, targetBest);
	
	// Property map for sum of weights for each node
	String<TCargo> weights;
	fill(weights, nseq, (TCargo) 0);
	reserve(weights, 2*nseq - 1);

	// Merge groups
	TSize m = nseq;
	while (m>1) {
		// Merge nodes
		TVertexDescriptor internalNode = addVertex(g);
		
		// Set the weights
		TCargo w = (TCargo) (minVal / 2);
		addEdge(g, internalNode, nodes[sourceBest], w - property(weights, nodes[sourceBest]));
		addEdge(g, internalNode, nodes[targetBest], w - property(weights, nodes[targetBest]));
		appendValue(weights, w);		

		// Get the new distance values
		_upgmaTreeMerge(pairGraph, active, best, TTag());
		
		// Inactivate one group, adjust the member count for the other one
		active[sourceBest] += active[targetBest];
		active[targetBest] = 0;
		nodes[sourceBest] = internalNode;
	
		// Update the minimum values
		minValues[sourceBest] = TWeightEdgePair(maxVal, 0);
		for(TEdgeOutI itOutE(pairGraph, sourceBest);!atEnd(itOutE);goNext(itOutE)) {
			TVD localTVD = targetVertex(itOutE);
			if (sourceBest < localTVD) {
				if (cargo(value(itOutE)) < (minValues[sourceBest].i1)) minValues[sourceBest] = TWeightEdgePair(cargo(value(itOutE)), localTVD);
			}
		}
		// Find the new minimum value
		itMin = begin(minValues, Standard() );
		itMinEnd = end(minValues, Standard() );
		minVal = maxVal;
		TVD oldSourceBest = sourceBest;
		sourceBest = 0;
		targetBest = 0;
		for(TVD index= 0;itMin != itMinEnd; goNext(itMin), ++index) {
			if (active[index] == 0) continue;
			if (((*itMin).i2 == oldSourceBest) || (active[(*itMin).i2] == 0)) {
				// Update the values
				(*itMin).i1 = maxVal;
				TEdgeOutI itOutLocal(pairGraph, index);
				for(;!atEnd(itOutLocal);goNext(itOutLocal)) {
					TVD targ = targetVertex(itOutLocal);
					if (targ < index) continue;
					if (cargo(value(itOutLocal)) < (*itMin).i1) *itMin = TWeightEdgePair(cargo(value(itOutLocal)), targ);
				}
			}
			if ((*itMin).i1 < minVal) {
				minVal = (*itMin).i1;
				sourceBest = index;
				targetBest = (*itMin).i2;
			}
		}
		// If none is found we have to insert a new edge
		if ((sourceBest == targetBest) && (m>2)) {
			TVertexI itV(pairGraph);
			sourceBest = value(itV); goNext(itV);
			targetBest = value(itV);
			best = addEdge(pairGraph, sourceBest, targetBest, infCargo);
		} else best = findEdge(pairGraph, sourceBest, targetBest);
		--m;
	}
	g.data_root = numVertices(g) - 1;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.upgmaTree:
..summary:Computes a guide tree from a distance matrix.
..cat:Graph
..signature:
upgmaTree(mat, graph [,tag])
upgmaTree(sparse_mat, graph [,tag])
..param.mat:A string of pairwise distance values, representing a square matrix.
...type:Class.String
...remarks: String must use double values because the algorithm recycles this string to store intermediate possibly fractional results.
..param.sparse_mat:An undirected graph where each edge corresponds to the distance between sequence i and sequence j.
...type:Spec.Undirected graph
..param.graph:Out-parameter:The guide tree.
...type:Spec.Tree
..param.tag:Tag that indicates how to calculate cluster distances.
...remarks:Possible values are UpgmaWeightAvg, UpgmaAvg, UpgmaMax, and UpgmaMin.
...type:Tag.Upgma Configurator
...default:UpgmaWeightAvg
..returns:void
*/
template<typename TDistance, typename TCargo, typename TSpec>
inline void
upgmaTree(TDistance& dist, 
		  Graph<Tree<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	upgmaTree(dist, g, UpgmaWeightAvg());
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
