#include "bam_merge.h"

#define USAGE "Usage: bam_merge [-Q] <out.bam> <in1.bam> <in2.bam> [...]\n"

void print_usage()
{
  fprintf(stderr, USAGE);
}

void write_bam_lines(GBamWriter& bw, vector<CBamLine>& bam_lines) {
  for (size_t i = 0; i < bam_lines.size(); ++i)
    {
      CBamLine& bam_line = bam_lines[i];
      bw.write(bam_line.b, bam_line.read_id);
      bam_line.b_free();
    }

  bam_lines.clear();
}

int main(int argc, char *argv[])
{
  int parse_ret = parse_options(argc, argv, print_usage);
  if (parse_ret)
    return parse_ret;
  if (quals) raw_merge=true; //hijack the -Q flag for this specific merging option
  char* outfname=NULL;
  if (argc-optind<3) {
    print_usage();
    if (argc>1)
      warn_msg("Error: only %d arguments given.\n", argc-1);
    return -1;
  }
  outfname=argv[optind];
  vector<string> bam_fnames;
  for (int i=optind+1;i<argc;i++)
    bam_fnames.push_back(argv[i]);

  BamMerge bamMerge(bam_fnames);
  GBamWriter bamwriter(outfname, bamMerge.get_sam_header(), index_outfile);

  vector<CBamLine> bam_lines;
  while (bamMerge.next_bam_lines(bam_lines))
    {
      write_bam_lines(bamwriter, bam_lines);
    }

  return 0;
}
