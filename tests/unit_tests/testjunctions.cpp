/*
 *  testjunctions.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/22/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#include <iostream>
#include <cppunit/extensions/TestFactoryRegistry.h> 
#include <cppunit/ui/text/TestRunner.h> 
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>


#include "bwt_map.h"
#include "junctions.h"

using namespace std;
using namespace seqan;

bool junction_found(const SpliceSet& splices, uint32_t left, uint32_t right)
{
	pair<size_t, size_t> p = make_pair<size_t, size_t>(left, right);
	SpliceSet::iterator i = splices.find(p);
	return i != splices.end();
}

class ClosureTest : public CppUnit::TestFixture  
{
	CPPUNIT_TEST_SUITE( ClosureTest );	
	CPPUNIT_TEST( testOneHopClosure );
	CPPUNIT_TEST( testNoClosure );
	CPPUNIT_TEST( testFuzzyOneHopClosure );
	CPPUNIT_TEST( testTwoHopExactClosure );
	CPPUNIT_TEST( testBowtieOverlapClosure );
	CPPUNIT_TEST_SUITE_END();
	
public:
	typedef String< Dna5, Alloc<> > Reference;
	void testOneHopClosure()
	{
		SpliceSet fwd_juncs;
		//                   insert/1             insert/2
		//                   |||||-------------------|||||"
		Reference ref_str = "CCCCCGTGGGGGGGGGGGGGGGAGTTTTT";
		
		SequenceTable it;
		SequenceTable rt;
		HitTable hits1(it,rt);
		HitTable hits2(it,rt);
		 
		hits1.add_hit("insert","ref",0,5,false);
		hits2.add_hit("insert","ref",24,5,true);
		
		vector<pair<size_t, size_t> > happy_mates;
		happy_mates.push_back(make_pair<size_t,size_t>(0,0));
		
		SimpleIntronFinder<Reference> intron_finder(ref_str, "GT", "AG");
		
		JunctionFinder<Reference, SimpleIntronFinder<Reference> > finder(ref_str, intron_finder, 10, 0, 10, 5);
		
		finder.possible_junctions(fwd_juncs,
								  happy_mates,
								  *(hits1.get_hits(0)),
								  *(hits2.get_hits(0)));
		
		CPPUNIT_ASSERT(junction_found(fwd_juncs, 4, 24));
		
	}
	
	void testNoClosure()
	{
		SpliceSet fwd_juncs;
		//                   insert/1             insert/2
		//                   |||||-------------------|||||"
		Reference ref_str = "CCCCCGTGGGGGGGGGGGGGGGGGTTTTT";
		
		SequenceTable it;
		SequenceTable rt;
		HitTable hits1(it,rt);
		HitTable hits2(it,rt);
		
		hits1.add_hit("insert","ref",0,5,false);
		hits2.add_hit("insert","ref",24,5,true);
		
		vector<pair<size_t, size_t> > happy_mates;
		happy_mates.push_back(make_pair<size_t,size_t>(0,0));
		
		SimpleIntronFinder<Reference> intron_finder(ref_str, "GT", "AG");
		
		JunctionFinder<Reference, SimpleIntronFinder<Reference> > finder(ref_str, intron_finder, 10, 0, 10, 5);
		
		finder.possible_junctions(fwd_juncs,
								  happy_mates,
								  *(hits1.get_hits(0)),
								  *(hits2.get_hits(0)));
		
		CPPUNIT_ASSERT(fwd_juncs.empty());
		
	}
	
	void testFuzzyOneHopClosure()
	{
		SpliceSet fwd_juncs;
		//                   insert/1             insert/2
		//                   |||||-------------------|||||"
		Reference ref_str = "CCCCCCCGTGGGGGGGGGGGAGTTTTTTT";
		
		SequenceTable it;
		SequenceTable rt;
		HitTable hits1(it,rt);
		HitTable hits2(it,rt);
		
		hits1.add_hit("insert","ref",0,5,false);
		hits2.add_hit("insert","ref",24,5,true);
		
		vector<pair<size_t, size_t> > happy_mates;
		happy_mates.push_back(make_pair<size_t,size_t>(0,0));
		
		SimpleIntronFinder<Reference> intron_finder(ref_str, "GT", "AG");
		
		JunctionFinder<Reference, SimpleIntronFinder<Reference> > finder1(ref_str, intron_finder, 10, 5, 10, 5);
		
		finder1.possible_junctions(fwd_juncs,
								  happy_mates,
								  *(hits1.get_hits(0)),
								  *(hits2.get_hits(0)));
		
		CPPUNIT_ASSERT(junction_found(fwd_juncs, 6, 22));
		
		fwd_juncs.clear();
		
		JunctionFinder<Reference, SimpleIntronFinder<Reference> >  finder2(ref_str, intron_finder, 10, 0, 10, 5);
		
		finder2.possible_junctions(fwd_juncs,
								   happy_mates,
								   *(hits1.get_hits(0)),
								   *(hits2.get_hits(0)));
		
		CPPUNIT_ASSERT(!junction_found(fwd_juncs, 6, 22));
		
	}
	
	void testTwoHopExactClosure()
	{
		SpliceSet fwd_juncs;
		//                   insert/1    single    insert/2
		//                               |||||
		//                   |||||-------------------|||||"
		Reference ref_str = "CCCCCGTGGGAGGGGGGGTGGGAGTTTTT";
		
		SequenceTable it;
		SequenceTable rt;
		HitTable hits1(it,rt);
		HitTable hits2(it,rt);
		
		hits1.add_hit("insert","ref",0,5,false);
		hits2.add_hit("insert","ref",24,5,true);
		
		// TODO: should really have a hit in the hit table for the singleton
		// For now, we find ALL donor and acceptors.  This test will break
		// when we start filtering donors and acceptors based on their proximity
		// to a region covered by the initially mapped reads.
		
		vector<pair<size_t, size_t> > happy_mates;
		happy_mates.push_back(make_pair<size_t,size_t>(0,0));
		
		SimpleIntronFinder<Reference> intron_finder(ref_str, "GT", "AG");
		
		JunctionFinder<Reference, SimpleIntronFinder<Reference> > finder(ref_str, intron_finder, 15, 0, 6, 5);
		
		finder.possible_junctions(fwd_juncs,
								  happy_mates,
								  *(hits1.get_hits(0)),
								  *(hits2.get_hits(0)));
		
		CPPUNIT_ASSERT(junction_found(fwd_juncs, 4, 12));
		CPPUNIT_ASSERT(junction_found(fwd_juncs, 16, 24));
	}
	
	void testBowtieOverlapClosure()
	{
		SpliceSet fwd_juncs;
		//                   insert/1             insert/2
		//                   |||||-------------------|||||"
		Reference ref_str = "CCCCCGTGGGGGGGGGGGGGGGAGTTTTT";
		
		SequenceTable it;
		SequenceTable rt;
		HitTable hits1(it,rt);
		HitTable hits2(it,rt);
		
		hits1.add_hit("insert","ref",0,7,false);
		hits2.add_hit("insert","ref",22,5,true);
		
		vector<pair<size_t, size_t> > happy_mates;
		happy_mates.push_back(make_pair<size_t,size_t>(0,0));
		
		SimpleIntronFinder<Reference> intron_finder(ref_str, "GT", "AG");
		
		JunctionFinder<Reference, SimpleIntronFinder<Reference> > finder(ref_str, intron_finder, 10, 0, 10, 5, 2);
		
		finder.possible_junctions(fwd_juncs,
								  happy_mates,
								  *(hits1.get_hits(0)),
								  *(hits2.get_hits(0)));
		
		CPPUNIT_ASSERT(junction_found(fwd_juncs, 4, 24));
	}
};

class IntronFinderTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE( IntronFinderTest );	
	CPPUNIT_TEST( testEmptyMap );
	CPPUNIT_TEST( testSingleInsertMap );
	CPPUNIT_TEST( testAggressiveTolerance );
	CPPUNIT_TEST_SUITE_END();
	
	typedef String< Dna5, Alloc<> > Reference;
	
public:
	void testEmptyMap()
	{
		//                   insert/1             insert/2
		//                   |||||-------------------|||||"
		Reference ref_str = "CCCCCGTGGGGGGGGGGGGGGGAGTTTTT";
		
		SequenceTable it;
		SequenceTable rt;
		HitTable hits1(it,rt);
		
		vector<const HitList*> all_hits;

		MappedIntronFinder<Reference> intron_finder(ref_str, all_hits, "GT", "AG", 5);
		
		CPPUNIT_ASSERT(intron_finder.l_hits().empty());
		CPPUNIT_ASSERT(intron_finder.r_hits().empty());
	}
	
	void testSingleInsertMap()
	{
		//                   insert/1             insert/2
		//                   |||||-------------------|||||"
		Reference ref_str = "CCCCCGTGGGGGGGGGGGGGGGAGTTTTT";
		
		SequenceTable it;
		SequenceTable rt;
		HitTable hits1(it,rt);
		HitTable hits2(it,rt);
		
		hits1.add_hit("insert","ref",0,5,false);
		hits2.add_hit("insert","ref",24,5,true);
		
		vector<const HitList*> all_hits;
		all_hits.push_back(hits1.get_hits(0));
		all_hits.push_back(hits2.get_hits(0));
		
		MappedIntronFinder<Reference> intron_finder(ref_str, all_hits, "GT", "AG", 5);
		
		const vector<size_t> lh = intron_finder.l_hits(); 
		const vector<size_t> rh = intron_finder.r_hits(); 
		
		CPPUNIT_ASSERT(lh.size() && lh.front() == 5 && lh.back() == 23);
		CPPUNIT_ASSERT(rh.size() && rh.back() == 22);
	}
	
	void testAggressiveTolerance()
	{
		//                   insert/1             insert/2
		//                   |||||-------------------|||||"
		Reference ref_str = "CCCCCGTGGGGGGGGGGGGGGGAGTTTTT";
		
		SequenceTable it;
		SequenceTable rt;
		HitTable hits1(it,rt);
		HitTable hits2(it,rt);
		
		hits1.add_hit("insert","ref",0,5,false);
		hits2.add_hit("insert","ref",24,5,true);
		
		vector<const HitList*> all_hits;
		all_hits.push_back(hits1.get_hits(0));
		all_hits.push_back(hits2.get_hits(0));
		
		MappedIntronFinder<Reference> intron_finder(ref_str, all_hits, "GT", "AG", 50);
		
		const vector<size_t> lh = intron_finder.l_hits(); 
		const vector<size_t> rh = intron_finder.r_hits(); 
		
		CPPUNIT_ASSERT(lh.size() && lh.front() == 5 && lh.back() == 23);
		CPPUNIT_ASSERT(rh.size() && rh.back() == 22);
	}
};

/*
 * Tests the use of a MappedIntronFinder with a Junction finder.
 */
class MappedJunctionFinder : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE( MappedJunctionFinder );	
	CPPUNIT_TEST( testTwoHopMappedClosure );
	CPPUNIT_TEST( testNoMapRouteClosure );
	CPPUNIT_TEST( testReusedJunctionClosure );
	CPPUNIT_TEST_SUITE_END();
	
	typedef String< Dna5, Alloc<> > Reference;
	
public:
	void testTwoHopMappedClosure()
	{
		SpliceSet fwd_juncs;
		//                   insert/1    single    insert/2
		//                               |||||
		//                   |||||-------------------|||||"
		Reference ref_str = "CCCCCGTGGGAGGGGGGGTGGGAGTTTTT";
		
		SequenceTable it;
		SequenceTable rt;
		HitTable hits1(it,rt);
		HitTable hits2(it,rt);
		
		hits1.add_hit("insert","ref",0,5,false);
		hits1.add_hit("single","ref",12, 17, false);
		hits2.add_hit("insert","ref",24,5,true);
		
		// TODO: should really have a hit in the hit table for the singleton
		// For now, we find ALL donor and acceptors.  This test will break
		// when we start filtering donors and acceptors based on their proximity
		// to a region covered by the initially mapped reads.
		
		vector<pair<size_t, size_t> > happy_mates;
		happy_mates.push_back(make_pair<size_t,size_t>(0,0));
		
		vector<const HitList*> all_hits;
		all_hits.push_back(hits1.get_hits(0));
		all_hits.push_back(hits2.get_hits(0));
		
		MappedIntronFinder<Reference> intron_finder(ref_str, all_hits, "GT", "AG", 3);
		
		JunctionFinder<Reference, MappedIntronFinder<Reference> > finder(ref_str, intron_finder, 15, 0, 6, 5);
		
		finder.possible_junctions(fwd_juncs,
								  happy_mates,
								  *(hits1.get_hits(0)),
								  *(hits2.get_hits(0)));
		
		CPPUNIT_ASSERT(junction_found(fwd_juncs, 4, 12));
		CPPUNIT_ASSERT(junction_found(fwd_juncs, 16, 24));
	}
	
	/* The lack of an alignment in the region indicated by "single"
	 * below means there is no route from insert/1 to insert/2. This 
	 * test checks that the junction finder doesn't find any route.
	 */
	void testNoMapRouteClosure()
	{
		SpliceSet fwd_juncs;
		//                   insert/1    single    insert/2
		//                               |||||
		//                   |||||-------------------|||||"
		Reference ref_str = "CCCCCGTGGGAGGGGGGGTGGGAGTTTTT";
		
		SequenceTable it;
		SequenceTable rt;
		HitTable hits1(it,rt);
		HitTable hits2(it,rt);
		
		hits1.add_hit("insert","ref",0,5,false);
		
		// Uncommenting this line should cause a failure for this test.
		//hits1.add_hit("single","ref",12, 17, false);
		hits2.add_hit("insert","ref",24,5,true);
		
		// TODO: should really have a hit in the hit table for the singleton
		// For now, we find ALL donor and acceptors.  This test will break
		// when we start filtering donors and acceptors based on their proximity
		// to a region covered by the initially mapped reads.
		
		vector<pair<size_t, size_t> > happy_mates;
		happy_mates.push_back(make_pair<size_t,size_t>(0,0));
		
		vector<const HitList*> all_hits;
		all_hits.push_back(hits1.get_hits(0));
		all_hits.push_back(hits2.get_hits(0));
		
		MappedIntronFinder<Reference> intron_finder(ref_str, all_hits, "GT", "AG", 3);
		
		JunctionFinder<Reference, MappedIntronFinder<Reference> > finder(ref_str, intron_finder, 15, 0, 6, 5);
		
		finder.possible_junctions(fwd_juncs,
								  happy_mates,
								  *(hits1.get_hits(0)),
								  *(hits2.get_hits(0)));
		
		CPPUNIT_ASSERT(fwd_juncs.empty());
	}
	
	void testReusedJunctionClosure()
	{
		SpliceSet fwd_juncs;
		//                   insert/1    single    insert/2
		//                               |||||
		//                   |||||-------------------|||||"
		Reference ref_str = "CCCCCGTGAGAGGGGGGGTGTGAGTTTTT";
		
		SequenceTable it;
		SequenceTable rt;
		HitTable hits1(it,rt);
		HitTable hits2(it,rt);
		
		hits1.add_hit("insert","ref",0,5,false);
		hits1.add_hit("single","ref",12, 17, false);
		hits2.add_hit("insert","ref",24,5,true);
		
		// TODO: should really have a hit in the hit table for the singleton
		// For now, we find ALL donor and acceptors.  This test will break
		// when we start filtering donors and acceptors based on their proximity
		// to a region covered by the initially mapped reads.
		
		vector<pair<size_t, size_t> > happy_mates;
		happy_mates.push_back(make_pair<size_t,size_t>(0,0));
		
		vector<const HitList*> all_hits;
		all_hits.push_back(hits1.get_hits(0));
		all_hits.push_back(hits2.get_hits(0));
		
		MappedIntronFinder<Reference> intron_finder(ref_str, all_hits, "GT", "AG", 4);
		
		JunctionFinder<Reference, MappedIntronFinder<Reference> > finder(ref_str, intron_finder, 15, 4, 6, 5);
		
		finder.possible_junctions(fwd_juncs,
								  happy_mates,
								  *(hits1.get_hits(0)),
								  *(hits2.get_hits(0)));
		
		CPPUNIT_ASSERT(junction_found(fwd_juncs, 4, 12));
		CPPUNIT_ASSERT(junction_found(fwd_juncs, 16, 24));
	}
	
};

class JunctionCompare : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE( JunctionCompare );	
	CPPUNIT_TEST( testIrreflexive );
	CPPUNIT_TEST( testHeirarchical );
	CPPUNIT_TEST( testAntisymmetric );
	CPPUNIT_TEST( testTransitive );
	CPPUNIT_TEST( testEquivalent );
	CPPUNIT_TEST( testSet );
	CPPUNIT_TEST_SUITE_END();
public:
	
	void testIrreflexive()
	{
		Junction j1;
		j1.refid = 0;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		
		CPPUNIT_ASSERT(j1 < j1 == false);
		
		j1.refid = 0;
		j1.left = 0;
		j1.right = 5;
		j1.antisense = false;
		
		CPPUNIT_ASSERT(j1 < j1 == false);
		
		j1.refid = 0;
		j1.left = 5;
		j1.right = 0;
		j1.antisense = false;
		
		CPPUNIT_ASSERT(j1 < j1 == false);
		
		j1.refid = 5;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		
		CPPUNIT_ASSERT(j1 < j1 == false);
	}
	
	void testHeirarchical()
	{
		Junction j1;
		j1.refid = 0;
		j1.left = 5;
		j1.right = 0;
		j1.antisense = false;
		
		Junction j2;
		j2.refid = 3;
		j2.left = 2;
		j2.right = 0;
		j2.antisense = false;
		CPPUNIT_ASSERT( !(j2 < j1) );
		
	}
	
	
	void testEquivalent()
	{
		Junction j1;
		j1.refid = 0;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		
		Junction j2;
		j2.refid = 0;
		j2.left = 0;
		j2.right = 0;
		j2.antisense = false;
		
		CPPUNIT_ASSERT(!(j1 < j2) && !(j2 < j1));
	}

	void testAntisymmetric()
	{
		Junction j1;
		j1.refid = 0;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		
		Junction j2;
		j2.refid = 1;
		j2.left = 0;
		j2.right = 0;
		j2.antisense = false;
		
		CPPUNIT_ASSERT(j1 < j2 && !(j2 < j1));
		
		j1.refid = 0;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;

		j2.refid = 0;
		j2.left = 1;
		j2.right = 0;
		j2.antisense = false;
		
		CPPUNIT_ASSERT(j1 < j2 && !(j2 < j1));
		
		j1.refid = 0;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		
		j2.refid = 0;
		j2.left = 0;
		j2.right = 1;
		j2.antisense = false;
		
		CPPUNIT_ASSERT(j1 < j2 && !(j2 < j1));
		
		j1.refid = 0;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		
		j2.refid = 0;
		j2.left = 0;
		j2.right = 0;
		j2.antisense = true;
		
		CPPUNIT_ASSERT(j1 < j2 && !(j2 < j1));
	}
	
	void testTransitive()
	{
		Junction j1;
		j1.refid = 0;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		
		Junction j2;
		j2.refid = 1;
		j2.left = 0;
		j2.right = 0;
		j2.antisense = false;
		
		Junction j3;
		j2.refid = 3;
		j2.left = 0;
		j2.right = 0;
		j2.antisense = false;
		
		CPPUNIT_ASSERT(j1 < j2 && j2 < j3 && j1 < j3);
		

		j1.refid = 0;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		
		j2.refid = 0;
		j2.left = 1;
		j2.right = 0;
		j2.antisense = false;
		
		j2.refid = 0;
		j2.left = 2;
		j2.right = 0;
		j2.antisense = false;
		
		CPPUNIT_ASSERT(j1 < j2 && j2 < j3 && j1 < j3);
		
		j1.refid = 0;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		
		j2.refid = 0;
		j2.left = 0;
		j2.right = 1;
		j2.antisense = false;
		
		j2.refid = 0;
		j2.left = 0;
		j2.right = 2;
		j2.antisense = false;
		
		CPPUNIT_ASSERT(j1 < j2 && j2 < j3 && j1 < j3);
		
		j1.refid = 0;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		
		j2.refid = 0;
		j2.left = 0;
		j2.right = 1;
		j2.antisense = false;
		
		j2.refid = 0;
		j2.left = 0;
		j2.right = 1;
		j2.antisense = true;
		
		CPPUNIT_ASSERT(j1 < j2 && j2 < j3 && j1 < j3);
	}
	
	void testSet()
	{
		Junction j1;
		
		std::set<Junction> junctions;
		
		j1.refid = 0;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		junctions.insert(j1);
		
		CPPUNIT_ASSERT(junctions.size() == 1);

		junctions.insert(j1);
		
		CPPUNIT_ASSERT(junctions.size() == 1);
		
		j1.refid = 1;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		junctions.insert(j1);
		
		CPPUNIT_ASSERT(junctions.size() == 2);
		
		j1.refid = 1;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = true;
		junctions.insert(j1);
		
		CPPUNIT_ASSERT(junctions.size() == 3);
		
		j1.refid = 1;
		j1.left = 1;
		j1.right = 0;
		j1.antisense = true;
		junctions.insert(j1);
		
		CPPUNIT_ASSERT(junctions.size() == 4);
		
		j1.refid = 1;
		j1.left = 1;
		j1.right = 1;
		j1.antisense = true;
		junctions.insert(j1);
		
		CPPUNIT_ASSERT(junctions.size() == 5);
		
		j1.refid = 1;
		j1.left = 1;
		j1.right = 1;
		j1.antisense = false;
		junctions.insert(j1);
		
		CPPUNIT_ASSERT(junctions.size() == 6);
		
		j1.refid = 1;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = false;
		junctions.insert(j1);
		
		CPPUNIT_ASSERT(junctions.size() == 6);
		
		j1.refid = 1;
		j1.left = 0;
		j1.right = 0;
		j1.antisense = true;
		junctions.insert(j1);
		
		CPPUNIT_ASSERT(junctions.size() == 6);
		
		j1.refid = 1;
		j1.left = 1;
		j1.right = 0;
		j1.antisense = true;
		junctions.insert(j1);
		
		CPPUNIT_ASSERT(junctions.size() == 6);
		
		j1.refid = 1;
		j1.left = 1;
		j1.right = 1;
		j1.antisense = true;
		junctions.insert(j1);
		
		CPPUNIT_ASSERT(junctions.size() == 6);
		
		j1.refid = 1;
		j1.left = 1;
		j1.right = 1;
		j1.antisense = false;
		junctions.insert(j1);
		
		CPPUNIT_ASSERT(junctions.size() == 6);
	}
	
};

class JunctionsFromSplicedHits : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE( JunctionsFromSplicedHits );	
	CPPUNIT_TEST( testSplicedHit );
	CPPUNIT_TEST_SUITE_END();
	
	typedef String< Dna5, Alloc<> > Reference;
	
public:
	
	void testSplicedHit()
	{

		SequenceTable it;
		SequenceTable rt;
		HitTable hits(it,rt);
		
		
		/*                    Reference                             */
		/* ----------------------|******************|---------------*/
		/* 0                    20                   40             */
		hits.add_spliced_hit("hit1", 
							 "reference", 
							 15, 
							 45, 
							 5, 
							 5, 
							 10, 
							 false, 
							 false);
		const HitList* hitlist = hits.get_hits(0);
		
		CPPUNIT_ASSERT(hitlist);
		CPPUNIT_ASSERT(hitlist->size());
		
		
		pair<Junction, JunctionStats> junc;
		junc = junction_from_spliced_hit(hitlist->front());
		Junction& j = junc.first;
		JunctionStats s = junc.second;
		
		CPPUNIT_ASSERT(j.refid == 0xFFFFFFFF);
		j.refid = 0;
		
		CPPUNIT_ASSERT(j.left == 20);
		CPPUNIT_ASSERT(j.right == 40);
		CPPUNIT_ASSERT(j.antisense == false);
	}
	

	
};

CPPUNIT_TEST_SUITE_REGISTRATION( ClosureTest );
CPPUNIT_TEST_SUITE_REGISTRATION( IntronFinderTest );
CPPUNIT_TEST_SUITE_REGISTRATION( MappedJunctionFinder );
CPPUNIT_TEST_SUITE_REGISTRATION( JunctionsFromSplicedHits );
CPPUNIT_TEST_SUITE_REGISTRATION( JunctionCompare );

int main( int argc, char*  argv[]) 
{
	// if command line contains "-selftest" then this is the post build check
	// => the output must be in the compiler error format.
	bool selfTest = (argc > 1)  &&  (std::string("-selftest") == argv[1]);
	
	
	CppUnit::TextUi::TestRunner runner;
	
	CppUnit::TestFactoryRegistry &registry =
	CppUnit::TestFactoryRegistry::getRegistry();
	
	runner.addTest( registry.makeTest() ); 
	
	if (selfTest)
	{	
		// Change the default outputter to a compiler error format outputter
		// The test runner owns the new outputter.
		runner.setOutputter( new CppUnit::CompilerOutputter( &runner.result(),
															std::cerr ) );
	}
	
	bool wasSuccessful = runner.run( "", !selfTest );
	return !wasSuccessful;
}