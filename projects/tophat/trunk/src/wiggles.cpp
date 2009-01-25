#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*
 *  wiggles.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 12/12/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */
#include <cassert>
#include <stdio.h>
#include <vector>
#include <string>
#include "wiggles.h"

using namespace std;

void print_wiggle_header(FILE* coverage_out)
{
	fprintf(coverage_out, "track type=wiggle_0 name=\"TopHat - read coverage\"\n");
}

void print_wiggle_for_ref(FILE* coverage_out,
						  const string& ref_name,
						  const vector<short>& DoC)
{
	short last_doc = 0; // Last DoC value we wrote to the file
	size_t last_pos = 0; // Postition where the last written DoC came from
	//fprintf(coverage_out, "variableStep chrom=chr%s\n", name.c_str());
	for (size_t i = 0; i < DoC.size(); ++i)
	{
		if (last_doc != DoC[i])
		{
			size_t j = last_pos;
			while (i - j > 10000000)
			{
				fprintf(coverage_out,"%s\t%d\t%d\t%d\n",ref_name.c_str(),(int)j + 1, (int)j + 10000000 + 1, last_doc);
				j += 10000000;
			}
			fprintf(coverage_out,"%s\t%d\t%d\t%d\n",ref_name.c_str(),(int)j + 1, (int)i + 1, last_doc);
			last_pos = i;
			last_doc = DoC[i];
		}
	}
}
