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

#ifndef SEQAN_HEADER_REPEAT_BASE_H
#define SEQAN_HEADER_REPEAT_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	template <typename TPos, typename TPeriod>
	struct Repeat {
		TPos		beginPosition;
		TPos		endPosition;
		TPeriod		period;
	};

	template <typename TPos, typename TPeriod>
	struct Value< Repeat<TPos, TPeriod> > {
		typedef TPos Type;
	};

	template <typename TPos, typename TPeriod>
	struct Size< Repeat<TPos, TPeriod> > {
		typedef TPeriod Type;
	};


	template <typename TSize>
	struct RepeatFinderParams {
		TSize minRepeatLen;
		TSize maxPeriod;
	};

	// custom TSpec for our customized wotd-Index
	struct TRepeatFinder;

	template <typename TText>
	struct Cargo<Index<TText, IndexWotd<TRepeatFinder> > > 
	{
		typedef Index<TText, IndexWotd<TRepeatFinder> >	TIndex;
		typedef typename Size<TIndex>::Type					TSize;
		typedef RepeatFinderParams<TSize>					Type;
	};


	// node predicate
	template <typename TText, typename TSpec>
	bool nodePredicate(Iter<Index<TText, IndexWotd<TRepeatFinder> >, TSpec> &it) 
	{
//		return countOccurrences(it) * nodeDepth(it) >= cargo(container(it)).minRepeatLen;
		return countOccurrences(it) * repLength(it) >= cargo(container(it)).minRepeatLen;
	}

	// monotonic hull
	template <typename TText, typename TSpec>
	bool nodeHullPredicate(Iter<Index<TText, IndexWotd<TRepeatFinder> >, TSpec> &it) 
	{
//		return nodeDepth(it) <= cargo(container(it)).maxPeriod;
		return repLength(it) <= cargo(container(it)).maxPeriod;
	}

	template <typename TPos>
	struct RepeatLess_ : public ::std::binary_function<TPos, TPos, bool>
	{
		// key less
		inline bool operator() (TPos const &a, TPos const &b) {
			return posLess(a, b);
		}
	};

	template <typename TValue>
	inline bool _repeatMaskValue(TValue) 
	{
		return false;
	}

	template <>
	inline bool _repeatMaskValue(Dna5 val) 
	{
		return val.value == 4; // 'N'
	}

	template <>
	inline bool _repeatMaskValue(Iupac val) 
	{
		static const Iupac n = 'N';
		return val.value == n.value;
	}
/*
	template <>
	inline bool _repeatMaskValue(AminoAcid val) 
	{
		return val == 'X';
	}
*/
	// period-1 optimization
	template <typename TRepeatStore, typename TString, typename TRepeatSize>
	inline void findRepeats(TRepeatStore &repString, TString const &text, TRepeatSize minRepeatLen) 
	{
		typedef typename Value<TRepeatStore>::Type	TRepeat;
		typedef typename Iterator<TString>::Type	TIterator;
		typedef typename Value<TString>::Type		TValue;
		typedef typename Size<TString>::Type		TSize;

		TRepeat rep;
		rep.period = 1;
		clear(repString);

		TIterator it = begin(text, Standard());
		TIterator itEnd = end(text, Standard());
		if (it == itEnd) return;

		TValue last = *it;
		TSize repLeft = 0;
		TSize repRight = 1;

		for (++it; it != itEnd; ++it, ++repRight) 
		{
			if (last != *it)
			{
				if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen)
				{
					// insert repeat
					rep.beginPosition = repLeft;
					rep.endPosition = repRight;
//					::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
					appendValue(repString, rep);
				}
				repLeft = repRight;
				last = *it;
			}
		}
		if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen)
		{
			// insert repeat
			rep.beginPosition = repLeft;
			rep.endPosition = repRight;
//			::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
			appendValue(repString, rep);
		}
	}

	template <typename TRepeatStore, typename TString, typename TSpec, typename TRepeatSize>
	inline void findRepeats(TRepeatStore &repString, StringSet<TString, TSpec> const &text, TRepeatSize minRepeatLen) 
	{
		typedef typename Value<TRepeatStore>::Type	TRepeat;
		typedef typename Iterator<TString>::Type	TIterator;
		typedef typename Value<TString>::Type		TValue;
		typedef typename Size<TString>::Type		TSize;

		TRepeat rep;
		rep.period = 1;
		clear(repString);

		for (unsigned i = 0; i < length(text); ++i)
		{
			TIterator it = begin(text[i], Standard());
			TIterator itEnd = end(text[i], Standard());
			if (it == itEnd) continue;

			TValue last = *it;
			TSize repLeft = 0;
			TSize repRight = 1;
			rep.beginPosition.i1 = i;
			rep.endPosition.i1 = i;

			for (++it; it != itEnd; ++it, ++repRight) 
			{
				if (last != *it)
				{
					if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen)
					{
						// insert repeat
						rep.beginPosition.i2 = repLeft;
						rep.endPosition.i2 = repRight;
//						::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
						appendValue(repString, rep);
					}
					repLeft = repRight;
					last = *it;
				}
			}
			if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen)
			{
				// insert repeat
				rep.beginPosition.i2 = repLeft;
				rep.endPosition.i2 = repRight;
//				::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
				appendValue(repString, rep);
			}
		}
	}

	// main function
	template <typename TRepeatStore, typename TText, typename TRepeatSize, typename TPeriodSize>
	void findRepeats(TRepeatStore &repString, TText const &text, TRepeatSize minRepeatLen, TPeriodSize maxPeriod) 
	{
		typedef Index<TText, IndexWotd<TRepeatFinder> >					TIndex;
		typedef typename Size<TIndex>::Type									TSize;
		typedef typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type	TNodeIterator;
		typedef typename Fibre<TIndex, FibreSA>::Type const				TSA;
		typedef typename Infix<TSA>::Type									TOccString;
		typedef typename Iterator<TOccString>::Type							TOccIterator;

		typedef typename Value<TRepeatStore>::Type							TRepeat;
		typedef typename Value<TOccString>::Type							TOcc;

		typedef ::std::map<TOcc,TRepeat,RepeatLess_<TOcc> >					TRepeatList;

		if (maxPeriod < 1) return;
		if (maxPeriod == 1) 
		{
			findRepeats(repString, text, minRepeatLen);
			return;
		}

		TIndex		index(text);
		TRepeatList list;

		// set repeat finder parameters
		cargo(index).minRepeatLen = minRepeatLen;
		cargo(index).maxPeriod = maxPeriod;

		TNodeIterator nodeIt(index);
		TOccIterator itA, itB, itRepBegin, itEnd;
		TRepeat rep;
		for (; !atEnd(nodeIt); goNext(nodeIt))
		{
			if (isRoot(nodeIt)) continue;

			// get occurrences
			TOccString occ = getOccurrences(nodeIt);
			itA = begin(occ, Standard());
			itEnd = end(occ, Standard());
			itRepBegin = itB = itA;

			TSize repLen = repLength(nodeIt);		// representative length
			if ((TSize)minRepeatLen <= repLen) continue;

			TSize diff, period = 0;					// period of current repeat
			TSize repeatLen = 0;					// overall length of current repeat
			TSize minLen = minRepeatLen - repLen;	// minimum repeat length minus length of representative

			for (++itB; itB != itEnd; ++itB)
			{
				diff = posSub(*itB, *itA);
				if (diff != period || getSeqNo(*itA) != getSeqNo(*itB))
				{
					// is the repeat long enough?
					if (repeatLen >= minLen)
						// is the repeat self overlapping or connected?
						if (parentRepLength(nodeIt) < period && period <= repLen)
						{
							// insert repeat
							rep.beginPosition = *itRepBegin;
							rep.endPosition = posAdd(*itA, period);
							rep.period = period;
//							::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
							list.insert(::std::pair<TOcc,TRepeat>(rep.beginPosition, rep));
						}
					itRepBegin = itA;
					period = diff;
					repeatLen = 0;
				}
				repeatLen += period;
				itA = itB;
			}

			// is the last repeat long enough?
			if (repeatLen >= minLen)
				// is the repeat self overlapping or connected?
				if (parentRepLength(nodeIt) < period && period <= repLen)
				{
					// insert repeat
					rep.beginPosition = *itRepBegin;
					rep.endPosition = posAdd(*itA, period);
					rep.period = period;
//					::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
					list.insert(::std::pair<TOcc,TRepeat>(rep.beginPosition, rep));
				}
		}

		// copy low-complex regions to result string
		resize(repString, list.size());
		typename TRepeatList::const_iterator lit = list.begin();
		typename TRepeatList::const_iterator litEnd = list.end();
		for (TSize i = 0; lit != litEnd; ++lit, ++i)
			repString[i] = (*lit).second;
	}


}	// namespace seqan

#endif
