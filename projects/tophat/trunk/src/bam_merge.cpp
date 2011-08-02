#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "bam/bam.h"
#include "bam/sam.h"
#include "GBase.h"
#include "GList.hh"

#define USAGE "Usage: bam_merge <out.bam> <in1.bam> <in2.bam> [...]\n"

#define ERR_BAM_OPEN "Error: bam_merge failed to open BAM file %s\n"

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
    char* outfname=NULL;
    if (argc<4) {
       fprintf(stderr, USAGE);
       if (argc>1)
         fprintf(stderr, "Error: only %d arguments given.\n", argc-1);
       return -1;
       }
    outfname=argv[1];
    int num_src=argc-2;
    GMALLOC(srcfiles, (num_src*sizeof(samfile_t*)));
    for (int i=2;i<argc;i++) {
       int fno=i-2;
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
    int last;
    while ((last=lines.Count()-1)>=0) {
    	CBamLine* from=lines[last]; //should have the smallest read_id
    	if (from->fileno<0 || from->b==NULL)
    		  GError("Invalid processTopLine() call with uninitialized value.");
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
    return 0;
}
