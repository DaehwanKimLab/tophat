#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <getopt.h>
#include "bwt_map.h"
#include "inserts.h"
#include "reads.h"
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;
using namespace std;

int max_gap = 50;
float min_cvg_fraction = 0.90;

void driver(FILE* bwt_map,  
			ifstream& ref_stream)
{
	typedef RefSequenceTable::Sequence Reference;
	
	ReadTable it;
	RefSequenceTable rt(true);
    HitTable hits1;

	fprintf(stderr, "Reading isoform sequences from FASTA file\n");
	map<uint32_t, vector<bool> > cov_maps;
	while(ref_stream.good() && 
		  !ref_stream.eof()) 
	{
		RefSequenceTable::Sequence* ref_str = new RefSequenceTable::Sequence;
		string name;
		readMeta(ref_stream, name, Fasta());
		read(ref_stream, *ref_str, Fasta());
		
		uint32_t ref_id = rt.get_id(name, NULL);
		cov_maps[ref_id] = vector<bool>(length(*ref_str), false);
	}
	
	char bwt_buf[2048];
	uint32_t reads_extracted = 0;
	
    BowtieHitFactory hit_factory(it,rt);
	fprintf(stderr, "Reading hits from Bowtie map\n");
	while (fgets(bwt_buf, 2048, bwt_map))
	{
		// Chomp the newline
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;
		if (*bwt_buf == 0)
			continue;
		// Get a new record from the tab-delimited Bowtie map
		BowtieHit bh;
		if (hit_factory.get_hit_from_buf(bwt_buf, bh, false))
		{
			map<uint32_t, vector<bool> >::iterator itr;
			itr = cov_maps.find(bh.ref_id());
			if (itr != cov_maps.end())
			{
				vector<bool>& cov_map = itr->second;
				for (int i = bh.left(); i < bh.right(); ++i)
				{
					if (i >= 0 && i < (int)cov_map.size())
						cov_map[i] = true;
				}
			}
		}
		reads_extracted++;
	}
	
	fprintf(stderr, "Reporting expressed isoforms\n");
	for (map<uint32_t, vector<bool> >::iterator itr = cov_maps.begin();
		 itr != cov_maps.end();
		 ++itr)
	{
		bool has_gap = false;
		int last_gap_start = 99999999;
		vector<bool>& cov_map = itr->second;
		
		int bases_covered = 0;
		for (int i = 0; i < (int)cov_map.size(); ++i)
		{
			if (!cov_map[i])
			{
				if (i > 0 && cov_map[i - 1])
					last_gap_start = i;
				else if (i - last_gap_start > max_gap)
				{
					//fprintf(stderr, "found gap %d bases long\n", i - last_gap_start);
					has_gap = true;
					break;
				}	
			}
			else
			{
				bases_covered++;
			}
		}
		
		float fraction_covered = bases_covered / (float)cov_map.size();
		
		if (!has_gap && fraction_covered > min_cvg_fraction)
		{
			const char* name = rt.get_name(itr->first);
			fprintf(stdout, "%s\n", name);
		}
	}
}

void print_usage()
{
    fprintf(stderr, "Usage:   expressed_isoform <isoforms.fa> <isoform_map.bwtout>\n");
}



/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static float parseFloat(float lower, float upper, const char *errmsg) {
    float l;
    l = (float)atof(optarg);
	
    if (l < lower) {
        cerr << errmsg << endl;
        print_usage();
        exit(1);
    }
	
    if (l > upper)
    {
        cerr << errmsg << endl;
        print_usage();
        exit(1);
    }
	
    return l;
	
    cerr << errmsg << endl;
    print_usage();
    exit(1);
    return -1;
}


const char *short_options = "m:f:";

static struct option long_options[] = {
{"max-gap",       required_argument,       0,            'm'},
{"min-cvg",       required_argument,       0,            'f'},
{0, 0, 0, 0} // terminator
};

int parse_options(int argc, char** argv)
{
    int option_index = 0;
    int next_option;
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
        switch (next_option) {
			case 'm':
	   			max_gap = (uint32_t)parseInt(1, "-m/--max-gap arg must be at least 1", print_usage);
	   			break;
			case 'f':
	   			min_cvg_fraction = (uint32_t)parseFloat(0.0, 1.0, "-f/--min-cvg arg must be at between 0.0 and 1.0");
	   			break;
            case -1:     /* Done with options. */
                break;
            default:
                print_usage();
                return 1;
        }
    } while(next_option != -1);
    
    return 0;
}


int main(int argc, char** argv)
{
    int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
	
	
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string ref_fasta = argv[optind++];
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string map_name = argv[optind++];
	
	
    FILE* bwt_map = fopen(map_name.c_str(), "r");
    if (bwt_map == NULL)
    {
        fprintf(stderr, "Error: cannot open %s for reading\n",
                map_name.c_str());
        exit(1);
    }

    ifstream ref_stream(ref_fasta.c_str(), ifstream::in);
	
    driver(bwt_map, ref_stream);
	
    return 0;
}