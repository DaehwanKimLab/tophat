#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <bam/bam.h>
#include <bam/sam.h>

#include "common.h"
#include "GBase.h"
#include "GList.hh"

#define USAGE "Usage: bam_merge <out.bam> <in1.bam> <in2.bam> [...]\n"

#define ERR_BAM_OPEN "Error: bam_merge failed to open BAM file %s\n"

void print_usage()
{
  fprintf(stderr, USAGE);
}

samfile_t **srcfiles; //array of SAM file handles
samfile_t *fw; //output SAM handle

class CBamLine {
 public:
   int fileno;
   long read_id;
   bam1_t* b;
   bool operator==(CBamLine& l){
     return (read_id==l.read_id);
     }
   bool operator>(CBamLine& l){
     return (read_id<l.read_id); //(last has lowest #)
     }
   bool operator<(CBamLine& l){
     return (read_id>l.read_id);
     }
   CBamLine(int fno=-1, bam1_t* br=NULL) {
     fileno=fno;
     read_id=0;
	 b=br;
     b_init();
     }
    void b_init() {
     if (b) {
       char *name  = bam1_qname(b);
       read_id=atol(name);
       if (read_id<1) {
    	  char* samline=bam_format1(srcfiles[0]->header, b);
    	  GError("Error: invalid read Id (must be numeric) for BAM record:\n%s\n",
    			  samline);
          }
       }
     }

   void b_free() {
     if (b!=NULL) {  bam_destroy1(b); }
     }
};

GList<CBamLine> lines(true,true,false);

int main(int argc, char *argv[])
{
  int parse_ret = parse_options(argc, argv, print_usage);
  if (parse_ret)
    return parse_ret;

    char* outfname=NULL;
    if (argc-optind<3) {
       print_usage();
       if (argc>1)
         fprintf(stderr, "Error: only %d arguments given.\n", argc-1);
       return -1;
       }
    outfname=argv[optind];
    int num_src=argc-optind-1;
    GMALLOC(srcfiles, (num_src*sizeof(samfile_t*)));
    for (int i=optind+1;i<argc;i++) {
       int fno=i-optind-1;
       samfile_t* fp=samopen(argv[i], "rb", 0);
       if (fp==0) {
               fprintf(stderr, ERR_BAM_OPEN, argv[i]);
               return 1;
               }
       bam1_t *b = bam_init1();
       if (samread(fp, b) > 0) {
    	  srcfiles[fno]=fp;
          lines.Add(new CBamLine(fno, b));
          }
       }
    if (lines.Count()==0) {
      GMessage("Warning: no input BAM records found.\n");
      }
    fw=samopen(outfname, "wb", srcfiles[lines[0]->fileno]->header);
    if (fw==NULL)
    	GError("Error creating output file %s\n", outfname);

    FILE* findex = NULL;
    if (!index_outfile.empty()) {
      findex = fopen(index_outfile.c_str(), "w");
      if (findex == NULL)
	err_die("Error: cannot create file %s\n", index_outfile.c_str());
    }

    int last;
    uint64_t count = 0, last_id = 0;
    while ((last=lines.Count()-1)>=0) {
    	CBamLine* from=lines[last]; //should have the smallest read_id
    	if (from->fileno<0 || from->b==NULL)
    		  GError("Invalid processTopLine() call with uninitialized value.");
	if (findex != NULL)
	  {
	    ++count;
	    if (count >= 1000 && from->read_id != last_id)
	      {
		int64_t offset = bam_tell(fw->x.bam);
		int block_offset = offset & 0xFFFF;
		
		// daehwan - this is a bug in bgzf.c in samtools
		// I'll report this bug with a temporary solution, soon.
		if (block_offset <= 0xF000)
		  {
		    fprintf(findex, "%lu\t%ld\n", from->read_id, offset);
		    count = 0;
		  }
	      }
	    last_id = from->read_id;
	  }
	samwrite(fw, from->b);

    	if (samread(srcfiles[from->fileno], from->b)>0) {
           from->b_init();
           //adjust the position in the sorted lines list
           if (last<7) {//
    		   int i=last;
    		   while (i>0 && lines[i-1]->read_id<lines[i]->read_id) {
    			   //swap
    			   CBamLine* tmp=lines[i-1];
    			   lines.Put(i-1, lines[i], false);
    			   lines.Put(i,tmp, false);
    			   i--;
    			   }
    			  }
              else { //use quick-insert
                lines.Pop();
                lines.Add(from);
              }
    	   }
    	else { //no more BAM records
           lines.Pop();
    	   from->b_free();
    	   }
       }
    samclose(fw);
    for (int i=0;i<num_src;i++)
        samclose(srcfiles[i]);
    GFREE(srcfiles);

    if (findex != NULL)
      fclose(findex);
    
    return 0;
}
