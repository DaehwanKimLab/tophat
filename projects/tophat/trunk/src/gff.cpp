#include "gff.h"

//GffNames* GffReader::names=NULL;
GffNames* GffObj::names=NULL;
//global set of feature names, attribute names etc.
// -- common for all GffObjs in current application!

const uint GFF_MAX_LOCUS = 7000000; //longest known gene in human is ~2.2M, UCSC claims a gene for mouse of ~ 3.1 M
const uint GFF_MAX_EXON  =   30000; //longest known exon in human is ~11K
const uint GFF_MAX_INTRON= 6000000;
bool gff_show_warnings = false; //global setting, set by GffReader->showWarnings()
const int gff_fid_mRNA=0;
const int gff_fid_transcript=1;
const int gff_fid_exon=2;
const int gff_fid_CDS=3; //never really used in GffObj ftype_id or subftype_id
const uint gfo_flag_HAS_ERRORS       = 0x00000001;
const uint gfo_flag_CHILDREN_PROMOTED= 0x00000002;
const uint gfo_flag_IS_GENE          = 0x00000004;
const uint gfo_flag_IS_TRANSCRIPT    = 0x00000008;
const uint gfo_flag_FROM_GFF3        = 0x00000080;
const uint gfo_flag_DISCARDED        = 0x00000100;
const uint gfo_flag_LST_KEEP         = 0x00000200;
const uint gfo_flag_LEVEL_MSK        = 0x00FF0000;
const byte gfo_flagShift_LEVEL           = 16;

void gffnames_ref(GffNames* &n) {
  if (n==NULL) n=new GffNames();
  n->numrefs++;
}

void gffnames_unref(GffNames* &n) {
  if (n==NULL) GError("Error: attempt to remove reference to null GffNames object!\n");
  n->numrefs--;
  if (n->numrefs==0) { delete n; n=NULL; }
}

int gfo_cmpByLoc(const pointer p1, const pointer p2) {
 
 GffObj& g1=*((GffObj*)p1);
 GffObj& g2=*((GffObj*)p2);
 if (g1.gseq_id==g2.gseq_id) {
             if (g1.start!=g2.start)
                    return (int)(g1.start-g2.start);
               else if (g1.getLevel()!=g2.getLevel())
                        return (int)(g1.getLevel()-g2.getLevel());
                    else
                        if (g1.end!=g2.end)
                              return (int)(g1.end-g2.end);
                        else return strcmp(g1.getID(), g2.getID());
             }
             else return (int)(g1.gseq_id-g2.gseq_id);
}

char* GffLine::extractAttr(const char* pre, bool caseStrict, bool enforce_GTF2) {
 //parse a key attribute and remove it from the info string
 //(only works for attributes that have values following them after ' ' or '=')
 static const char GTF2_ERR[]="Error parsing attribute %s ('\"' required) at GTF line:\n%s\n";
 int lpre=strlen(pre);
 char cend=pre[lpre-1];
 char* pos = (caseStrict) ? strstr(info, pre) : strifind(info, pre);
 if (pos==NULL) return NULL;
 char* findstart=info;
 //require word boundary on the left:
 while (pos!=NULL && pos!=info && *(pos-1)!=';' && *(pos-1)!=' ') {
    findstart=pos+lpre;
    pos = (caseStrict) ? strstr(findstart, pre) : strifind(findstart, pre);
    }
 if (pos==NULL) return NULL;
 if (cend!=' ' && cend!='=') {
    //require word boundary on the right:
    while (pos!=NULL && *(pos+lpre)!=' ' && *(pos+lpre)!='=') {
       findstart=pos+lpre;
       pos = (caseStrict) ? strstr(findstart, pre) : strifind(findstart, pre);
       }
    }
 if (pos==NULL) return NULL;
 char* vp=pos+lpre;
 while (*vp!=';' && *vp!=0 && (*vp==' ' || *vp=='"')) vp++;
 if (*vp==';' || *vp==0)
      GError("Error parsing value of GFF attribute \"%s\", line:\n%s\n", pre, dupline);
 if (enforce_GTF2 && *(vp-1)!='"')
      GError(GTF2_ERR,pre, dupline);
 char* vend=vp;
 while (*vend!=';' && *vend!=0 && *vend!='"') vend++;
 if (enforce_GTF2 && *vend!='"')
     GError(GTF2_ERR, pre, dupline);
 char *r=Gstrdup(vp, vend-1);
 //-- now remove this attribute from the info string
 //while (pos>info && *pos==' ') pos--;
 while (*vend!=0 && (*vend=='"' || *vend==';' || *vend==' ')) vend++;
 if (*vend==0) vend--;
 for (char *src=vend, *dest=pos;;src++,dest++) {
   *dest=*src;
   if (*src==0) break;
   }
 return r;
}

static char fnamelc[128];

GffLine::GffLine(GffReader* reader, const char* l) {
 llen=strlen(l);
 GMALLOC(line,llen+1);
 memcpy(line, l, llen+1);
 GMALLOC(dupline, llen+1);
 memcpy(dupline, l, llen+1);
 skip=true;
 gseqname=NULL;
 track=NULL;
 ftype=NULL;
 info=NULL;
 _parents=NULL;
 _parents_len=0;
 num_parents=0;
 parents=NULL;
 is_gff3=false;
 is_cds=false;
 is_transcript=false;
 is_exon=false;
 is_gene=false;
 exontype=0;
 gname=NULL;
 qstart=0;
 qend=0;
 qlen=0;
 ID=NULL;
 char* t[9];
 int i=0;
 int tidx=1;
 t[0]=line;
 
 while (line[i]!=0) {
  if (line[i]=='\t') {
   line[i]=0;
   t[tidx]=line+i+1;
   tidx++;
   if (tidx>8) break;
   }
  i++;
  }

 if (tidx<8) { // ignore non-GFF lines
  // GMessage("Warning: error parsing GFF/GTF line:\n%s\n", l);
  return;
  }
 gseqname=t[0];
 track=t[1];
 ftype=t[2];
 info=t[8];
 char* p=t[3];
 if (!parseUInt(p,fstart))
   GError("Error parsing start coordinate from GFF line:\n%s\n",l);
 p=t[4];
 if (!parseUInt(p,fend))
   GError("Error parsing end coordinate from GFF line:\n%s\n",l);
 if (fend<fstart) swap(fend,fstart); //make sure fstart>=fend, always
 p=t[5];
 if (p[0]=='.' && p[1]==0) {
  score=0;
  }
 else {
  if (!parseDouble(p,score))
       GError("Error parsing feature score from GFF line:\n%s\n",l);
  }
 strand=*t[6];
 if (strand!='+' && strand!='-' && strand!='.')
     GError("Error parsing strand (%c) from GFF line:\n%s\n",strand,l);
 phase=*t[7]; // must be '.', '0', '1' or '2'
 ID=NULL;
 // exon/CDS/mrna filter
 strncpy(fnamelc, ftype, 127);
 fnamelc[127]=0;
 strlower(fnamelc); //convert to lower case
 bool is_t_data=false;
 if (strstr(fnamelc, "utr")!=NULL) {
   exontype=exgffUTR;
   is_exon=true;
   is_t_data=true;
   }
  else if (strstr(fnamelc, "exon")!=NULL) {
   exontype=exgffExon;
   is_exon=true;
   is_t_data=true;
   }
  else if (strstr(fnamelc, "stop") && 
      (strstr(fnamelc, "codon") || strstr(fnamelc, "cds"))){
   exontype=exgffStop;
   is_cds=true; //though some place it outside the last CDS segment
   is_t_data=true;
   }
  else if (strstr(fnamelc, "start") && 
      ((strstr(fnamelc, "codon")!=NULL) || strstr(fnamelc, "cds")!=NULL)){
   exontype=exgffStart;
   is_cds=true;
   is_t_data=true;
   }
 else if (strcmp(fnamelc, "cds")==0) {
   exontype=exgffCDS;
   is_cds=true;
   is_t_data=true;
   }
 else if (endsWith(fnamelc, "gene") || startsWith(fnamelc, "gene")) {
   is_gene=true;
   is_t_data=true; //because its name will be attached to parented transcripts
   }
 else if (endsWith(fnamelc,"rna") || endsWith(fnamelc,"transcript")) {
   is_transcript=true;
   is_t_data=true;
   }

if (reader->transcriptsOnly && !is_t_data) {
        char* id=extractAttr("ID=");
        if (id==NULL) id=extractAttr("transcript_id");
        //GMessage("Discarding non-transcript line:\n%s\n",l);
        if (id!=NULL) {
          reader->discarded_ids.Add(id, new int(1));
          GFREE(id);
          }
        return; //skip this line, unwanted feature name
        }
 ID=extractAttr("ID=");
 char* Parent=extractAttr("Parent=");
 is_gff3=(ID!=NULL || Parent!=NULL);
 if (is_gff3) {
   //parse as GFF3
    if (ID!=NULL) {
       //has ID attr so it's likely to be a parent feature
       //look for explicit gene name
       gname=extractAttr("gene_name=",false);
       if (gname==NULL) {
           gname=extractAttr("gene_id=",false);
           if (gname==NULL) {
               gname=extractAttr("gene=",false);
               if (gname==NULL) {
                   gname=extractAttr("geneid=",false);
                   }
               }
           }
       if (gname==NULL && is_gene) {
          //special case: take Name= attribute as the gene name, if it exists
          gname=extractAttr("Name=");
          if (gname==NULL) {//no Name attribute, use the ID
             gname=Gstrdup(ID);
             }
          skip=false;
          return;
          }
       }// has GFF3 ID
   if (Parent!=NULL) {
        //keep Parent attr
         //parse multiple parents
         num_parents=1;
         p=Parent;
         int last_delim_pos=-1;
         while (*p!=';' && *p!=0) {
             if (*p==',' && *(p+1)!=0 && *(p+1)!=';') {
                 num_parents++;
                 last_delim_pos=(p-Parent);
                 }
             p++;
             }
         _parents_len=p-Parent+1;
         _parents=Parent;
         GMALLOC(parents, num_parents*sizeof(char*));
         parents[0]=_parents;
         int i=1;
         if (last_delim_pos>0) {
           for (p=_parents+1;p<=_parents+last_delim_pos;p++) {
              if (*p==',') {
                 char* ep=p-1;
                 while (*ep==' ' && ep>_parents) ep--;
                 *(ep+1)=0; //end the string there
                 parents[i]=p+1;
                 i++;
                 }
              }
           }
         }
   } //GFF3
  else { // GTF-like expected
   Parent=extractAttr("transcript_id");
   if (Parent!=NULL) { //GTF2 format detected
     if (is_transcript) {
         // atypical GTF with a parent transcript line declared
         ID=Parent;
         Parent=NULL;
         }
     //check for gene_id
     gname=extractAttr("gene_id");
     //prepare for parseAttr by adding '=' character instead of spaces for all attributes
     //after the attribute name
     p=info;
     bool noed=true; //not edited after the last delim
     bool nsp=false; //non-space found after last delim
     while (*p!=0) {
       if (*p==' ') {
          if (nsp && noed) {
             *p='=';
             noed=false;
             p++;
             continue;
             }
           }
         else nsp=true; //non-space
       if (*p==';') { noed=true; nsp=false; }
       p++;
       }
     } //GTF2 detected (no parent line)
    else {// Parent is NULL, check for jigsaw format or other pre-GTF2 format
     //char* fexon=strstr(fnamelc, "exon");
     //if (fexon!=NULL) {
     if (exontype==exgffExon) {
       if (startsWith(track,"jigsaw")) {
          is_cds=true;
          strcpy(track,"jigsaw");
          p=strchr(info,';');
          if (p==NULL) { Parent=Gstrdup(info); info=NULL; }
           else { Parent=Gstrdup(info,p-1);
                  info=p+1;
                }
          }
        } //exon feature?
        if (Parent==NULL && exontype>=exgffCDS &&
               (i=strcspn(info,"; \t\n\r"))<=(int)(strlen(info)+1)) {
          //one word ID ? really desperate attempt to parse it here
          Parent=Gstrdup(info,info+i-1);
          info=NULL; //discard anything else on the line
          }
     }
   if (Parent!=NULL) { //GTF transcript_id for exon/CDS feature
      _parents=Parent;
      GMALLOC(parents,sizeof(char*));
      num_parents=1;
      parents[0]=_parents;
      }
   } //GTF-like

 //parse other potentially useful features
 if (is_gff3) {
   if ((p=strstr(info,"Target="))!=NULL) { //has Target attr
      p+=7;
      while (*p!=';' && *p!=0 && *p!=' ') p++;
      if (*p!=' ') {
         GError("Error parsing target coordinates from GFF line:\n%s\n",l);
         }
      if (!parseUInt(p,qstart))
         GError("Error parsing target start coordinate from GFF line:\n%s\n",l);
      if (*p!=' ') {
         GError("Error parsing next target coordinate from GFF line:\n%s\n",l);
         }
      p++;
      if (!parseUInt(p,qend))
         GError("Error parsing target end coordinate from GFF line:\n%s\n",l);
      }
   if ((p=strifind(info,"Qreg="))!=NULL) { //has Qreg attr
       p+=5;
       if (!parseUInt(p,qstart))
         GError("Error parsing target start coordinate from GFF line:\n%s\n",l);
       if (*p!='-') {
          GError("Error parsing next target coordinate from GFF line:\n%s\n",l);
          }
       p++;
       if (!parseUInt(p,qend))
         GError("Error parsing target end coordinate from GFF line:\n%s\n",l);
       if (*p=='|' || *p==':') {
         p++;
         if (!parseUInt(p,qlen))
           GError("Error parsing target length from GFF Qreg|: \n%s\n",l);
         }
       }//has Qreg attr
   if (qlen==0 && (p=strifind(info,"Qlen="))!=NULL) {
     p+=5;
     if (!parseUInt(p,qlen))
         GError("Error parsing target length from GFF Qlen:\n%s\n",l);
     }
   }//parsing some useful attributes in GFF3 records
 if (ID==NULL && parents==NULL) {
      if (reader->gff_warns)
          GMessage("Warning: could not parse ID or Parent from GFF line:\n%s\n",dupline);
      return; //skip
      }
 skip=false;
}

int GffObj::addExon(GffReader* reader, GffLine* gl, bool keepAttr, bool noExonAttr) {
  //this will make sure we have the right subftype_id!
  int subf_id=-1;
  //if (ftype_id==gff_fid_mRNA) { //for mRNAs only parse known subfeatures!
  if (isTranscript()) {
     if (exon_ftype_id<0) {//exon_ftype_id=gff_fid_exon;
          if (gl->exontype>0) exon_ftype_id=gff_fid_exon;
                         else exon_ftype_id=names->feats.addName(gl->ftype);
          }
     //any recognized mRNA segment gets the generic "exon" type (also applies to CDS)
     if (gl->exontype==0 && !gl->is_transcript) {
          //extraneous mRNA feature, discard
          if (reader->gff_warns)
            GMessage("Warning: discarding unrecognized transcript subfeature %s of %s\n", 
                gl->ftype, gffID);
          return -1;
          }
     }
  else { //non-mRNA parent feature, check this subf type
    subf_id=names->feats.addName(gl->ftype);
    if (exon_ftype_id<0 || exons.Count()==0) //never assigned a subfeature type before (e.g. first exon being added)
       exon_ftype_id=subf_id;
     else {
       if (exon_ftype_id!=subf_id) {
         //if (subftype_id==ftype_id && exons.Count()==1 && exons[0]->start==start && exons[0]->end==end) {
         if (exon_ftype_id==ftype_id && exons.Count()==1 && exons[0]->start==start && exons[0]->end==end) {
            //the existing exon was just a dummy one created by default, discard it
            exons.Clear();
            covlen=0;
            exon_ftype_id=subf_id; //allow the new subfeature to completely takeover
            }
         else { //multiple subfeatures, prefer those with 
             if (reader->gff_warns)
               GMessage("GFF Warning: multiple subfeatures (%s and %s) found for %s, discarding ", 
                  names->feats.getName(subf_id), names->feats.getName(exon_ftype_id),gffID);
            if (gl->exontype!=0) { //new feature is an exon, discard previously parsed subfeatures
               if (reader->gff_warns) GMessage("%s.\n", names->feats.getName(exon_ftype_id));
               exon_ftype_id=subf_id;
               exons.Clear();
               covlen=0;
               }
              else { //discard new feature
               if (reader->gff_warns) GMessage("%s.\n", names->feats.getName(subf_id));
               return -1; //skip this 2nd subfeature type for this parent!
               }
            }
         } //incoming subfeature is of different type
       } //new subfeature type
    } //non-mRNA parent
  int eidx=addExon(gl->fstart, gl->fend, gl->score, gl->phase,
         gl->qstart,gl->qend, gl->is_cds, gl->exontype);
  if (eidx<0) return eidx; //this should never happen
  if (keepAttr) {
     if (noExonAttr) { 
         if (attrs==NULL) //place the parsed attributes directly at transcript level
           parseAttrs(attrs, gl->info, noExonAttr);
         }
       else { //need all exon-level attributes
         parseAttrs(exons[eidx]->attrs, gl->info, noExonAttr);
         }
      }
  return eidx;
}


int GffObj::addExon(uint segstart, uint segend, double sc, char fr, int qs, int qe, bool iscds, char exontype) {
  if (exons.Count()==0) {
      if (iscds) isCDS=true; //for now, assume CDS only if first "exon" given is a CDS
      if (exon_ftype_id<0) {
         exon_ftype_id = isTranscript() ? gff_fid_exon : ftype_id;
         }
      }
  //special treatment of start/stop codon features, they might be broken/split between exons 
  //and in that case some providers will still give the wrong end coordinate as start+2 (e.g. UCSC)
  //so we should not trust the end coordinate for such features
  if (exontype==exgffStart || exontype==exgffStop) {
     if (strand=='-') segstart=segend;
                else  segend=segstart;
     if (exontype==exgffStart) {
           if (CDstart==0 || segstart<CDstart) CDstart=segstart;
           }
         else {
           if (segstart>CDend) CDend=segstart;
           }
     }
    else if (iscds) { //update CDS anchors:
     if (CDstart==0 || segstart<CDstart)  {
           CDstart=segstart;
           if (exontype==exgffCDS && strand=='+') CDphase=fr;
           }
     if (segend>CDend) {
           if (exontype==exgffCDS && strand=='-') CDphase=fr;
           CDend=segend;
           }
     }
   else { // not a CDS/start/stop 
     isCDS=false;
     }
  if (qs || qe) {
    if (qs>qe) swap(qs,qe);
    if (qs==0) qs=1;
    }
  if (exontype>0) { //check for overlaps between exon-type segments
      int ovlen=0;
      int oi=exonOverlapIdx(segstart, segend, &ovlen);
      if (oi>=0) { //overlap existing segment
         if (ovlen==0) {
              //adjacent segments will be merged
              if ((exons[oi]->exontype==exgffUTR && exontype==exgffCDS) ||
                  (exons[oi]->exontype==exgffCDS && exontype==exgffUTR)) {
                    expandExon(oi, segstart, segend, exgffCDSUTR, sc, fr, qs, qe);
                    return oi;
                    }
             }
         //only allow this for CDS within exon, stop_codon within exon, stop_codon within UTR,
         //                   start_codon within CDS or stop_codon within CDS
        if (exons[oi]->exontype>exontype && 
             exons[oi]->start<=segstart && exons[oi]->end>=segend &&
             !(exons[oi]->exontype==exgffUTR && exontype==exgffCDS)) {
              //larger segment given first, now the smaller included one is redundant
              return oi; //only used to store attributes from current GffLine
              }
        if (exontype>exons[oi]->exontype && 
             segstart<=exons[oi]->start && segend>=exons[oi]->end &&
             !(exontype==exgffUTR && exons[oi]->exontype==exgffCDS)) {
               //smaller segment given first, so we have to enlarge it
              expandExon(oi, segstart, segend, exontype, sc, fr, qs, qe); 
                //this should also check for overlapping next exon (oi+1) ?
              return oi;
              }
        //there is also the special case of "ribosomal slippage exception" (programmed frameshift)
        //where two CDS segments may actually overlap for 1 or 2 bases, but there should be only one encompassing exon
        //if (ovlen>2 || exons[oi]->exontype!=exgffCDS || exontype!=exgffCDS) {
        // --> had to relax this because of some weird UCSC annotations with exons partially overlapping the CDS segments
        if (ovlen>2) {
           //important structural warning, will always print:
           if (gff_show_warnings) 
               GMessage("GFF Warning: discarding overlapping feature segment (%d-%d) (vs %d-%d (%s)) for GFF ID %s on %s\n", 
               segstart, segend, exons[oi]->start, exons[oi]->end, getSubfName(), gffID, getGSeqName());
           hasErrors(true);
           return -1; //segment NOT added
           }
          // else add the segment if the overlap is small and between two CDS segments
          //we might want to add an attribute here with the slippage coordinate and size?
        }//overlap of existing segment
       } //check for overlap
   // --- no overlap, or accepted micro-overlap (ribosomal slippage)
   // create & add the new segment
   GffExon* enew=new GffExon(segstart, segend, sc, fr, qs, qe, exontype);
   int eidx=exons.Add(enew);
   if (eidx<0) {
    //this would actually be acceptable if the object is a "Gene" and "exons" are in fact isoforms
     if (gff_show_warnings) 
       GMessage("GFF Warning: failed adding segment %d-%d for %s (discarded)!\n",
            segstart, segend, gffID);
     delete enew;
     hasErrors(true);
     return -1;            
     }
   covlen+=(int)(exons[eidx]->end-exons[eidx]->start)+1;
   start=exons.First()->start;
   end=exons.Last()->end;
   if (uptr!=NULL) { //collect stats about the underlying genomic sequence
       GSeqStat* gsd=(GSeqStat*)uptr;
       if (start<gsd->mincoord) gsd->mincoord=start;
       if (end>gsd->maxcoord) gsd->maxcoord=end;
       if (this->len()>gsd->maxfeat_len) {
          gsd->maxfeat_len=this->len();
          gsd->maxfeat=this;
          }
       }
   return eidx;
}

void GffObj::expandExon(int oi, uint segstart, uint segend, char exontype, double sc, char fr, int qs, int qe) {
  //oi is the index of the *first* overlapping segment found that must be enlarged
  covlen-=exons[oi]->len();
  if (segstart<exons[oi]->start)
    exons[oi]->start=segstart;
  if (qs && qs<exons[oi]->qstart) exons[oi]->qstart=qs;
  if (segend>exons[oi]->end)
    exons[oi]->end=segend;
  if (qe && qe>exons[oi]->qend) exons[oi]->qend=qe;
  //warning: score cannot be properly adjusted! e.g. if it's a p-value it's just going to get worse
  if (sc!=0) exons[oi]->score=sc;
  covlen+=exons[oi]->len();
  //if (exons[oi]->exontype< exontype) -- always true
  exons[oi]->exontype = exontype;
  if (exontype==exgffCDS) exons[oi]->phase=fr;
  //we must check if any more exons are also overlapping this 
  int ni=oi+1; //next exon index after oi
  while (ni<exons.Count() && segend>=exons[ni]->start) { // next segment overlaps new enlarged segment
     //only allow this if next segment is fully included, and a subordinate
     if (exons[ni]->exontype<exontype && exons[ni]->end<=segend) {
/* I guess we have to relax this due to stupid UCSC hg18 files having a start_codon sticking out
chr1	hg18_knownGene	start_codon	69806911	69806913	0.000000	+	.
chr1	hg18_knownGene	CDS	69806911	69806912	0.000000	+	0
chr1	hg18_knownGene	exon	69805456	69806912	0.000000	+	.     
*/
         if (exons[ni]->qstart<exons[oi]->qstart) exons[oi]->qstart=exons[ni]->qstart;
         if (exons[ni]->qend>exons[oi]->qend) exons[oi]->qend=exons[ni]->qend;
         exons.Delete(ni);
         }
      else {
         if (gff_show_warnings) GMessage("GFF Warning: overlapping existing exon(%d-%d) while expanding to %d-%d for GFF ID %s\n",
                exons[ni]->start, exons[ni]->end, segstart, segend, gffID);
         //hasErrors(true);
         break;
         }
     }
  // -- make sure any other related boundaries are updated:
  start=exons.First()->start;
  end=exons.Last()->end;
  if (uptr!=NULL) { //collect stats about the underlying genomic sequence
    GSeqStat* gsd=(GSeqStat*)uptr;
    if (start<gsd->mincoord) gsd->mincoord=start;
    if (end>gsd->maxcoord) gsd->maxcoord=end;
    if (this->len()>gsd->maxfeat_len) {
        gsd->maxfeat_len=this->len();
        gsd->maxfeat=this;
        }
    }
}

void GffObj::removeExon(int idx) {
  /*
   if (idx==0 && segs[0].start==gstart)
                  gstart=segs[1].start;
   if (idx==segcount && segs[segcount].end==gend)
                  gend=segs[segcount-1].end;
  */
  if (idx<0 || idx>=exons.Count()) return;
  int segstart=exons[idx]->start;
  int segend=exons[idx]->end;
  exons.Delete(idx);
  covlen -= (int)(segend-segstart)+1;
  start=exons.First()->start;
  end=exons.Last()->end;
  if (isCDS) { CDstart=start; CDend=end; }
}

void GffObj::removeExon(GffExon* p) {
  for (int idx=0;idx<exons.Count();idx++) {
     if (exons[idx]==p) {
        int segstart=exons[idx]->start;
        int segend=exons[idx]->end;
        exons.Delete(idx);
        covlen -= (int)(segend-segstart)+1;
        start=exons.First()->start;
        end=exons.Last()->end;
        if (isCDS) { CDstart=start; CDend=end; }
        return;
        }
     }
}



GffObj::GffObj(GffReader *gfrd, GffLine* gffline, bool keepAttr, bool noExonAttr):
     GSeg(0,0), exons(true,true,false), children(1,false) {
  xstart=0;
  xend=0;
  xstatus=0;
  partial=false;
  isCDS=false;
  uptr=NULL;
  ulink=NULL;
  parent=NULL;
  udata=0;
  flags=0;
  CDstart=0;
  CDend=0;
  CDphase=0;
  gname=NULL;
  attrs=NULL;
  gffID=NULL;
  track_id=-1;
  gseq_id=-1;
  ftype_id=-1;
  exon_ftype_id=-1;
  strand='.';
  if (gfrd==NULL)
    GError("Cannot use this GffObj constructor with a NULL GffReader!\n");
  gffnames_ref(names);
  if (gfrd->names==NULL) gfrd->names=names;
  //qlen=0;qstart=0;qend=0;
  gscore=0;
  uscore=0;
  covlen=0;
  qcov=0;
  start=gffline->fstart;
  end=gffline->fend;
  gseq_id=names->gseqs.addName(gffline->gseqname);
  track_id=names->tracks.addName(gffline->track);
  strand=gffline->strand;
  qlen=gffline->qlen;
  qstart=gffline->qstart;
  qend=gffline->qend;
  //setup flags from gffline
  isCDS=gffline->is_cds; //for now
  isGene(gffline->is_gene);
  isTranscript(gffline->is_transcript || gffline->exontype!=0);
  fromGff3(gffline->is_gff3);

  if (gffline->parents!=NULL) {
    //GTF style -- create a GffObj directly by subfeature
    //(also possible orphan GFF3 exon line)
    if (gffline->exontype!=0) { //recognized exon-like feature
       ftype_id=gff_fid_transcript; //so this is some sort of transcript
       exon_ftype_id=gff_fid_exon; //subfeatures MUST be exons
       }
     else {//unrecognized subfeatures
       //make this GffObj of the same feature type
       ftype_id=names->feats.addName(gffline->ftype);
       }
    if (gffline->ID==NULL) { //typical GTF
        gffID=Gstrdup(gffline->parents[0]);
        //this is likely the first exon/segment of the feature
        addExon(gfrd, gffline, keepAttr, noExonAttr);
        }
      else { //a parented feature with an ID -- probably an orphan GFF3 line
        //just save the attributes but don't add itself as exon
        gffID=Gstrdup(gffline->ID);
        if (keepAttr) this->parseAttrs(attrs, gffline->info, noExonAttr);
        }
    } //subfeature given directly
  else { //gffline->parents==NULL
    //create a parent feature in its own right
    gscore=gffline->score;
    if (gffline->ID==NULL || gffline->ID[0]==0)
      GError("Error: no ID found for GFF record start\n");
    gffID=Gstrdup(gffline->ID); //there must be an ID here
    //if (gffline->is_transcript) ftype_id=gff_fid_mRNA;
      //else
    ftype_id=names->feats.addName(gffline->ftype);
    if (gffline->is_transcript)
      exon_ftype_id=gff_fid_exon;

    if (keepAttr) this->parseAttrs(attrs, gffline->info, noExonAttr);
    }//no parent

  if (gffline->gname!=NULL) {
     gname=Gstrdup(gffline->gname);
     }

  GSeqStat* gsd=gfrd->gseqstats.AddIfNew(new GSeqStat(gseq_id,names->gseqs.lastNameUsed()),true);
  uptr=gsd;
  if (start<gsd->mincoord) gsd->mincoord=start;
  if (end>gsd->maxcoord) gsd->maxcoord=end;
    if (this->len()>gsd->maxfeat_len) {
        gsd->maxfeat_len=this->len();
        gsd->maxfeat=this;
        }
}

GffLine* GffReader::nextGffLine() {
 if (gffline!=NULL) return gffline; //caller should free gffline after processing
 while (gffline==NULL) {
    int llen=0;
    buflen=GFF_LINELEN-1;
    char* l=fgetline(linebuf, buflen, fh, &fpos, &llen);
    if (l==NULL) {
         return NULL; //end of file
         }
    int ns=0; //first nonspace position
    while (l[ns]!=0 && isspace(l[ns])) ns++;
    if (l[ns]=='#' || llen<10) continue;
    gffline=new GffLine(this, l);
    if (gffline->skip) {
       delete gffline;
       gffline=NULL;
       continue;
       }
    if (gffline->ID==NULL && gffline->parents==NULL)  { //it must have an ID
        //this might not be needed, already checked in the GffLine constructor
        if (gff_warns)
            GMessage("Warning: malformed GFF line, no parent or record Id (kipping\n");
        delete gffline;
        gffline=NULL;
        //continue;
        }
    }
return gffline;
}

char* GffReader::gfoBuildId(const char* id, const char* ctg) {
//caller must free the returned pointer
 char* buf=NULL;
 int idlen=strlen(id);
 GMALLOC(buf, idlen+strlen(ctg)+2);
 strcpy(buf, id);
 buf[idlen]='~';
 strcpy(buf+idlen+1, ctg);
 return buf;
}

void GffReader::gfoRemove(const char* id, const char* ctg) {
 char* buf=gfoBuildId(id,ctg);
 phash.Remove(buf);
 GFREE(buf);
}

//Warning: if gflst gets altered, idx becomes obsolete
GfoHolder* GffReader::gfoAdd(const char* id, const char* ctg, GffObj* gfo, int idx) {
 char* buf=gfoBuildId(id,ctg);
 GfoHolder* r=new GfoHolder(gfo,idx);
 phash.Add(buf, r);
 GFREE(buf);
 return r;
}

GfoHolder* GffReader::gfoFind(const char* id, const char* ctg) {
 char* buf=gfoBuildId(id,ctg);
 GfoHolder* r=phash.Find(buf);
 GFREE(buf);
 return r;
}

GfoHolder* GffReader::replaceGffRec(GffLine* gffline, bool keepAttr, bool noExonAttr, int replaceidx) {
  GffObj* newgfo=new GffObj(this, gffline, keepAttr, noExonAttr);
  GfoHolder* r=NULL;
  if (replaceidx>=0) {
     gflst.Put(replaceidx,newgfo);
     r=gfoAdd(newgfo->gffID, gffline->gseqname, newgfo, replaceidx);
     }
   else {
     int gfoidx=gflst.Add(newgfo);
     r=gfoAdd(newgfo->gffID, gffline->gseqname, newgfo, gfoidx);
     }
  if (gff_warns) {
    int* pcount=tids.Find(newgfo->gffID);
    if (pcount!=NULL) {
       if (gff_warns) GMessage("Warning: duplicate GFF ID: %s\n", newgfo->gffID);
       (*pcount)++;
       }
     else {
       tids.Add(newgfo->gffID,new int(1));
       }
    }
  return r;
}


GfoHolder* GffReader::newGffRec(GffLine* gffline, bool keepAttr, bool noExonAttr,
                          GffObj* parent, GffExon* pexon) {
  GffObj* newgfo=new GffObj(this, gffline, keepAttr, noExonAttr);
  GfoHolder* r=NULL;
  int gfoidx=gflst.Add(newgfo);
  r=gfoAdd(newgfo->gffID, gffline->gseqname, newgfo, gfoidx);
  if (parent!=NULL) {
    parent->children.Add(newgfo);
    newgfo->parent=parent;
    newgfo->setLevel(parent->getLevel()+1);
    if (parent->isGene() && parent->gname!=NULL && newgfo->gname==NULL)
       newgfo->gname=Gstrdup(parent->gname);
    if (pexon!=NULL) parent->removeExon(pexon);
    }
  if (gff_warns) {
    int* pcount=tids.Find(newgfo->gffID);
    if (pcount!=NULL) {
       if (gff_warns) GMessage("Warning: duplicate GFF ID: %s\n", newgfo->gffID);
       (*pcount)++;
       }
     else {
       tids.Add(newgfo->gffID,new int(1));
       }
    }
  return r;
}

bool GffReader::addExonFeature(GfoHolder* prevgfo, GffLine* gffline, GHash<CNonExon>& pex, bool noExonAttr) {
  bool r=true;
  if (gffline->strand!=prevgfo->gffobj->strand) {
     GMessage("GFF Error: duplicate GFF ID '%s' (exons found on different strands of %s)\n",
        prevgfo->gffobj->gffID, prevgfo->gffobj->getGSeqName());
      r=false;
     }
  int gdist=(gffline->fstart>prevgfo->gffobj->end) ? gffline->fstart-prevgfo->gffobj->end :
                      ((gffline->fend<prevgfo->gffobj->start)? prevgfo->gffobj->start-gffline->fend :
                         0 );
  if (gdist>(int)GFF_MAX_LOCUS) { //too far apart, most likely this is a duplicate ID
    GMessage("Error: duplicate GFF ID '%s' (or exons too far apart)!\n",prevgfo->gffobj->gffID);
    //validation_errors = true;
    r=false;
    if (!gff_warns) exit(1);
    }
  int eidx=prevgfo->gffobj->addExon(this, gffline, !noExonAttr, noExonAttr);
  if (eidx>=0 && gffline->ID!=NULL && gffline->exontype==0)
      subfPoolAdd(pex, prevgfo);
  return r;
}

CNonExon* GffReader::subfPoolCheck(GffLine* gffline, GHash<CNonExon>& pex, char*& subp_name) {
  CNonExon* subp=NULL;
  subp_name=NULL;
  for (int i=0;i<gffline->num_parents;i++) {
    if (transcriptsOnly && discarded_ids.Find(gffline->parents[i])!=NULL)
        continue;
    subp_name=gfoBuildId(gffline->parents[i], gffline->gseqname); //e.g. mRNA name
    subp=pex.Find(subp_name);
    if (subp!=NULL)
       return subp;
    GFREE(subp_name);
    }
  return NULL;
}

void GffReader::subfPoolAdd(GHash<CNonExon>& pex, GfoHolder* newgfo) {
//this might become a parent feature later
if (newgfo->gffobj->exons.Count()>0) {
   char* xbuf=gfoBuildId(gffline->ID, gffline->gseqname);
   pex.Add(xbuf, new CNonExon(newgfo->idx, newgfo->gffobj,
       newgfo->gffobj->exons[0], gffline));
   GFREE(xbuf);
   }
}

GfoHolder* GffReader::promoteFeature(CNonExon* subp, char*& subp_name, GHash<CNonExon>& pex,
    bool keepAttr, bool noExonAttr) {
  GffObj* prevp=subp->parent; //grandparent of gffline (e.g. gene)
  if (prevp!=gflst[subp->idx])
    GError("Error promoting subfeature %s, gflst index mismatch?!\n", subp->gffline->ID);
  subp->gffline->discardParent();
  GfoHolder* gfoh=newGffRec(subp->gffline, keepAttr, noExonAttr, prevp, subp->exon);
  pex.Remove(subp_name); //no longer a potential parent, moved it to phash already
  prevp->promotedChildren(true);
  return gfoh; //returns the holder of newly promoted feature
}

//have to parse the whole file because exons can be scattered all over
void GffReader::readAll(bool keepAttr, bool mergeCloseExons, bool noExonAttr) {
  bool validation_errors = false;
  //loc_debug=false;
  GHash<CNonExon> pex; //keep track of any "exon"-like features that have an ID
                     //and thus could become promoted to parent features
  while (nextGffLine()!=NULL) {
    if (gffline->parents==NULL) {//start GFF3-like record with no parent
       //had this gff ID before?
       GfoHolder* f=gfoFind(gffline->ID, gffline->gseqname);
       if (f!=NULL) {
            GMessage("Error: duplicate GFF ID '%s' encountered!\n",gffline->ID);
            validation_errors = true;
            if (gff_warns) { delete gffline; gffline=NULL; continue; }
                        else exit(1);
            }
       newGffRec(gffline, keepAttr, noExonAttr);
       }
    else { //--- it's a parented feature:
       bool found_parent=false;
       GfoHolder* newgfo=NULL;
       for (int i=0;i<gffline->num_parents;i++) {
            if (transcriptsOnly && discarded_ids.Find(gffline->parents[i])!=NULL)
               continue;
            GfoHolder* prevgfo=gfoFind(gffline->parents[i], gffline->gseqname);
            if (prevgfo!=NULL) { //parent GffObj
                   found_parent=true;
                   if (prevgfo->gffobj->isGene() && gffline->is_transcript
                                   && gffline->exontype==0) {
                       //not an exon, but a transcript parented by a gene
                       if (newgfo!=NULL) {
                           prevgfo->gffobj->children.Add(newgfo->gffobj);
                           //newgfo->gffobj->parent already set
                           }
                         else {
                           newgfo=newGffRec(gffline, keepAttr, noExonAttr, prevgfo->gffobj);                           }
                           }
                   else { //potential exon subfeature
                       if (!addExonFeature(prevgfo, gffline, pex, noExonAttr))
                         validation_errors=true;
                       }
                   }
            } //for each parsed parent Id
       if (!found_parent) { //new GTF-like record starting here with a subfeature directly
             //check if this feature isn't parented by a previously stored "exon" subfeature
            char* subp_name=NULL;
            CNonExon* subp=subfPoolCheck(gffline, pex, subp_name);
            if (subp!=NULL) { //found a subfeature that is the parent of this gffline
               //promote that subfeature to a full GffObj
               GfoHolder* gfoh=promoteFeature(subp, subp_name, pex, keepAttr, noExonAttr);
               //add current gffline as an exon of the newly promoted subfeature
               if (!addExonFeature(gfoh, gffline, pex, noExonAttr))
                      validation_errors=true;
               }
              else { //no parent seen before, create one directly with this exon
               //loc_debug=true;
               GfoHolder* newgfo=newGffRec(gffline, keepAttr, noExonAttr);
               if (gffline->ID!=NULL && gffline->exontype==0)
                     subfPoolAdd(pex, newgfo);
               //even those with errors will be added here!
               }
            GFREE(subp_name);
            } //no previous parent found
       } //parented feature
        //--
      delete gffline;
      gffline=NULL;
      }//while gff lines
  gflst.finalize(this, mergeCloseExons); //force sorting by locus if so constructed
 // all gff records are now loaded in GList gflst
 // so we can free the hash
  phash.Clear();
  tids.Clear();
  if (validation_errors) {
    exit(1);
    }
}

GffObj* GffObj::finalize(GffReader* gfr, bool mergeCloseExons) {
 //merge
 //always merge adjacent or overlapping segments
 //but if mergeCloseExons then merge even when distance is up to 5 bases
 udata=0;
 uptr=NULL;
 if (gfr->transcriptsOnly && !(isTranscript() || (isGene() && children.Count()==0))) {
       isDiscarded(true);
       }
 if (ftype_id==gff_fid_transcript && CDstart>0) {
    ftype_id=gff_fid_mRNA;
    //exon_ftype_id=gff_fid_exon;
    }
 //if (ftype_id==gff_fid_mRNA || exon_ftype_id==gff_fid_exon || mergeCloseExons) {
 if (isTranscript() || exon_ftype_id==gff_fid_exon || mergeCloseExons) {
   int mindist=mergeCloseExons ? 5:1;
   for (int i=0;i<exons.Count()-1;i++) {
     int ni=i+1;
     uint mend=exons[i]->end;
     while (ni<exons.Count()) {
       int dist=(int)(exons[ni]->start-mend);
       if (dist>mindist) break; //no merging with next segment
       if (gfr!=NULL && gfr->gff_warns) {
          GMessage("GFF warning: merging adjacent/overlapping segments of %s on %s (%d-%d, %d-%d)\n",
               gffID, getGSeqName(), exons[i]->start, exons[i]->end,exons[ni]->start, exons[ni]->end);
          }
       mend=exons[ni]->end;
       covlen-=exons[i]->len();
       exons[i]->end=mend;
       covlen+=exons[i]->len();
       covlen-=exons[ni]->len();
       if (exons[ni]->attrs!=NULL && (exons[i]->attrs==NULL || 
            exons[i]->attrs->Count()<exons[ni]->attrs->Count())) {
              //use the other exon attributes, if more
              delete(exons[i]->attrs);
              exons[i]->attrs=exons[ni]->attrs;
              exons[ni]->attrs=NULL;
              }
       exons.Delete(ni);
       } //check for merge with next exon
     } //for each exon
   }
 return this;
}

void GffObj::parseAttrs(GffAttrs*& atrlist, char* info, bool noExonAttr) {
  if (names==NULL)
     GError(ERR_NULL_GFNAMES, "parseAttrs()");
  if (atrlist==NULL)
      atrlist=new GffAttrs();
  char* endinfo=info+strlen(info);
  char* start=info;
  char* pch=start;
  while (start<endinfo) {
    //skip spaces
    while (*start==' ' && start<endinfo) start++;
    pch=strchr(start, ';');
    if (pch==NULL) pch=endinfo;
       else {
            *pch='\0';
            pch++;
            }
    char* ech=strchr(start,'=');
    if (ech!=NULL) { // attr=value format found
       *ech='\0';
       /*
       if (strcmp(start, "Target")==0 || strcmp(start,"Parent")==0 ||
           strcmp(start,"transcript_id")==0 || strcmp(start,"gene_id")==0)
            { start=pch; continue; } //skip these already parsed attributes
       */
       if (noExonAttr && (strcmp(start, "exon_number")==0 || strcmp(start, "exon")==0)) { start=pch; continue; }
       ech++;
       while (*ech==' ' && ech<endinfo) ech++;//skip extra spaces after the '='
       atrlist->Add(new GffAttr(names->attrs.addName(start),ech));
       }
      /*
      else { //not an attr=value format
        atrlist->Add(new GffAttr(names->attrs.addName(start),"1"));
        }
      */
    start=pch;
    }
  if (atrlist->Count()==0) { delete atrlist; atrlist=NULL; }
}

void GffObj::addAttr(const char* attrname, char* attrvalue) {
  if (this->attrs==NULL)
      this->attrs=new GffAttrs();
  this->attrs->Add(new GffAttr(names->attrs.addName(attrname),attrvalue));
}

int GffObj::removeAttr(const char* attrname, const char* attrval) {
  if (this->attrs==NULL || attrname==NULL || attrname[0]==0) return 0;
  int aid=this->names->attrs.getId(attrname);
  if (aid<0) return 0;
  int delcount=0;  //could be more than one ?
  for (int i=0;i<this->attrs->Count();i++) {
     if (aid==this->attrs->Get(i)->attr_id) {
       if (attrval==NULL || 
          strcmp(attrval, this->attrs->Get(i)->attr_val)==0) {
             delcount++;
             this->attrs->freeItem(i);
             }
       }
     }
  if (delcount>0) this->attrs->Pack(); 
  return delcount;
}

void GffObj::getCDS_ends(uint& cds_start, uint& cds_end) {
  cds_start=0;
  cds_end=0;
  if (CDstart==0 || CDend==0) return; //no CDS info
  int cdsadj=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsadj=CDphase-'0';
      }
  cds_start=CDstart;
  cds_end=CDend;
  if (strand=='-') cds_end-=cdsadj;
              else cds_start+=cdsadj;
  }

void GffObj::mRNA_CDS_coords(uint& cds_mstart, uint& cds_mend) {
  //sets cds_start and cds_end to the CDS start,end coordinates on the spliced mRNA transcript
  cds_mstart=0;
  cds_mend=0;
  if (CDstart==0 || CDend==0) return; //no CDS info
  //restore normal coordinates, just in case
  unxcoord();
  int cdsadj=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsadj=CDphase-'0';
      }
  /*
   uint seqstart=CDstart;
   uint seqend=CDend;
  */
  uint seqstart=exons.First()->start;
  uint seqend=exons.Last()->end;
  int s=0; //resulting nucleotide counter
  if (strand=='-') {
    for (int x=exons.Count()-1;x>=0;x--) {
       uint sgstart=exons[x]->start;
       uint sgend=exons[x]->end;
       if (seqend<sgstart || seqstart>sgend) continue;
       if (seqstart>=sgstart && seqstart<=sgend)
             sgstart=seqstart; //seqstart within this segment
       if (seqend>=sgstart && seqend<=sgend)
             sgend=seqend; //seqend within this segment
       s+=(int)(sgend-sgstart)+1;
       if (CDstart>=sgstart && CDstart<=sgend) {
             //CDstart in this segment
             //and we are getting the whole transcript
             cds_mend=s-(int)(CDstart-sgstart);
             }
       if (CDend>=sgstart && CDend<=sgend) {
             //CDstart in this segment
             //and we are getting the whole transcript
             cds_mstart=s-(int)(CDend-cdsadj-sgstart);
             }
      } //for each exon
    } // - strand
   else { // + strand
    for (int x=0;x<exons.Count();x++) {
      uint sgstart=exons[x]->start;
      uint sgend=exons[x]->end;
      if (seqend<sgstart || seqstart>sgend) continue;
      if (seqstart>=sgstart && seqstart<=sgend)
            sgstart=seqstart; //seqstart within this segment
      if (seqend>=sgstart && seqend<=sgend)
            sgend=seqend; //seqend within this segment
      s+=(int)(sgend-sgstart)+1;
      /* for (uint i=sgstart;i<=sgend;i++) {
          spliced[s]=gsubseq[i-gstart];
          s++;
          }//for each nt
          */
      if (CDstart>=sgstart && CDstart<=sgend) {
            //CDstart in this segment
            cds_mstart=s-(int)(sgend-CDstart-cdsadj);
            }
      if (CDend>=sgstart && CDend<=sgend) {
            //CDend in this segment
            cds_mend=s-(int)(sgend-CDend);
            }
      } //for each exon
    } // + strand
  //spliced[s]=0;
  //if (rlen!=NULL) *rlen=s;
  //return spliced;
}

char* GffObj::getUnspliced(GFaSeqGet* faseq, int* rlen, GList<GSeg>* seglst) 
{
    if (faseq==NULL) { GMessage("Warning: getUnspliced(NULL,.. ) called!\n");
        return NULL;
    }
    //restore normal coordinates:
    unxcoord();
    if (exons.Count()==0) return NULL;
    int fspan=end-start+1;
    const char* gsubseq=faseq->subseq(start, fspan);
    if (gsubseq==NULL) {
        GError("Error getting subseq for %s (%d..%d)!\n", gffID, start, end);
    }
    char* unspliced=NULL;

    int seqstart=exons.First()->start;
    int seqend=exons.Last()->end;
    
    int unsplicedlen = 0;

    unsplicedlen += seqend - seqstart + 1;

    GMALLOC(unspliced, unsplicedlen+1); //allocate more here
    //uint seqstart, seqend;

    int s = 0; //resulting nucleotide counter
    if (strand=='-') 
    {
        if (seglst!=NULL)
            seglst->Add(new GSeg(s+1,s+1+seqend-seqstart));
        for (int i=seqend;i>=seqstart;i--) 
        {
            unspliced[s] = ntComplement(gsubseq[i-start]);
            s++;
        }//for each nt
    } // - strand
    else 
    { // + strand
        if (seglst!=NULL)
            seglst->Add(new GSeg(s+1,s+1+seqend-seqstart));
        for (int i=seqstart;i<=seqend;i++) 
        {
            unspliced[s]=gsubseq[i-start];
            s++;
        }//for each nt
    } // + strand
    //assert(s <= unsplicedlen);
    unspliced[s]=0;
    if (rlen!=NULL) *rlen=s;
    return unspliced;
}

char* GffObj::getSpliced(GFaSeqGet* faseq, bool CDSonly, int* rlen, uint* cds_start, uint* cds_end,
          GList<GSeg>* seglst) {
  if (CDSonly && CDstart==0) return NULL;
  if (faseq==NULL) { GMessage("Warning: getSpliced(NULL,.. ) called!\n");
              return NULL;
              }
  //restore normal coordinates:
  unxcoord();
  if (exons.Count()==0) return NULL;
  int fspan=end-start+1;
  const char* gsubseq=faseq->subseq(start, fspan);
  if (gsubseq==NULL) {
        GError("Error getting subseq for %s (%d..%d)!\n", gffID, start, end);
        }
  if (fspan<(int)(end-start+1)) { //special case: stop coordinate was extended past the gseq length, must adjust
     int endadj=end-start+1-fspan;
     uint prevend=end;
     end-=endadj;
     if (CDend>end) CDend=end;
     if (exons.Last()->end>end) {
         exons.Last()->end=end; //this could get us into trouble if exon start is also > end
         if (exons.Last()->start>exons.Last()->end) {
            GError("GffObj::getSpliced() error: improper genomic coordinate %d on %s for %s\n",
                  prevend,getGSeqName(), getID());
            }
         covlen-=endadj;
         }
     }
  char* spliced=NULL;
  GMALLOC(spliced, covlen+1); //allocate more here
  uint seqstart, seqend;
  int cdsadj=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsadj=CDphase-'0';
      }
  if (CDSonly) {
     seqstart=CDstart;
     seqend=CDend;
     if (strand=='-') seqend-=cdsadj;
           else seqstart+=cdsadj;
     }
   else {
     seqstart=exons.First()->start;
     seqend=exons.Last()->end;
     }
  int s=0; //resulting nucleotide counter
  if (strand=='-') {
    for (int x=exons.Count()-1;x>=0;x--) {
       uint sgstart=exons[x]->start;
       uint sgend=exons[x]->end;
       if (seqend<sgstart || seqstart>sgend) continue;
       if (seqstart>=sgstart && seqstart<=sgend)
             sgstart=seqstart; //seqstart within this segment
       if (seqend>=sgstart && seqend<=sgend)
             sgend=seqend; //seqend within this segment
       if (seglst!=NULL)
           seglst->Add(new GSeg(s+1,s+1+sgend-sgstart));
       for (uint i=sgend;i>=sgstart;i--) {
            spliced[s] = ntComplement(gsubseq[i-start]);
            s++;
            }//for each nt

       if (!CDSonly && cds_start!=NULL && CDstart>0) {
          if (CDstart>=sgstart && CDstart<=sgend) {
             //CDstart in this segment
             //and we are getting the whole transcript
             *cds_end=s-(CDstart-sgstart);
             }
          if (CDend>=sgstart && CDend<=sgend) {
             //CDstart in this segment
             //and we are getting the whole transcript
             *cds_start=s-(CDend-cdsadj-sgstart);
             }
         }//update local CDS coordinates
      } //for each exon
    } // - strand
   else { // + strand
    for (int x=0;x<exons.Count();x++) {
      uint sgstart=exons[x]->start;
      uint sgend=exons[x]->end;
      if (seqend<sgstart || seqstart>sgend) continue;
      if (seqstart>=sgstart && seqstart<=sgend)
            sgstart=seqstart; //seqstart within this segment
      if (seqend>=sgstart && seqend<=sgend)
            sgend=seqend; //seqend within this segment
      if (seglst!=NULL)
          seglst->Add(new GSeg(s+1,s+1+sgend-sgstart));
      for (uint i=sgstart;i<=sgend;i++) {
          spliced[s]=gsubseq[i-start];
          s++;
          }//for each nt
      if (!CDSonly && cds_start!=NULL && CDstart>0) {
         if (CDstart>=sgstart && CDstart<=sgend) {
            //CDstart in this segment
            //and we are getting the whole transcript
            *cds_start=s-(sgend-CDstart-cdsadj);
            }
         if (CDend>=sgstart && CDend<=sgend) {
            //CDstart in this segment
            //and we are getting the whole transcript
            *cds_end=s-(sgend-CDend);
            }
        }//update local CDS coordinates
      } //for each exon
    } // + strand
  spliced[s]=0;
  if (rlen!=NULL) *rlen=s;
  return spliced;
}

char* GffObj::getSplicedTr(GFaSeqGet* faseq, bool CDSonly, int* rlen) {
  if (CDSonly && CDstart==0) return NULL;
  //restore normal coordinates:
  unxcoord();
  if (exons.Count()==0) return NULL;
  int fspan=end-start+1;
  const char* gsubseq=faseq->subseq(start, fspan);
  if (gsubseq==NULL) {
    GError("Error getting subseq for %s (%d..%d)!\n", gffID, start, end);
    }

  char* translation=NULL;
  GMALLOC(translation, (int)(covlen/3)+1);
  uint seqstart, seqend;
  int cdsadj=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsadj=CDphase-'0';
      }
  if (CDSonly) {
     seqstart=CDstart;
     seqend=CDend;
     if (strand=='-') seqend-=cdsadj;
           else seqstart+=cdsadj;
     }
   else {
     seqstart=exons.First()->start;
     seqend=exons.Last()->end;
     }
  Codon codon;
  int nt=0; //codon nucleotide counter (0..2)
  int aa=0; //aminoacid count
  if (strand=='-') {
    for (int x=exons.Count()-1;x>=0;x--) {
       uint sgstart=exons[x]->start;
       uint sgend=exons[x]->end;
       if (seqend<sgstart || seqstart>sgend) continue;
       if (seqstart>=sgstart && seqstart<=sgend)
             sgstart=seqstart; //seqstart within this segment
       if (seqend>=sgstart && seqend<=sgend) {
             sgend=seqend; //seqend within this segment
             }
       for (uint i=sgend;i>=sgstart;i--) {
            codon.nuc[nt]=ntComplement(gsubseq[i-start]);
            nt++;
            if (nt==3) {
               nt=0;
               translation[aa]=codon.translate();
               aa++;
               }
            }//for each nt
      } //for each exon
    } // - strand
   else { // + strand
    for (int x=0;x<exons.Count();x++) {
      uint sgstart=exons[x]->start;
      uint sgend=exons[x]->end;
      if (seqend<sgstart || seqstart>sgend) continue;
      if (seqstart>=sgstart && seqstart<=sgend)
            sgstart=seqstart; //seqstart within this segment
      if (seqend>=sgstart && seqend<=sgend)
            sgend=seqend; //seqend within this segment
      for (uint i=sgstart;i<=sgend;i++) {
          codon.nuc[nt]=gsubseq[i-start];
          nt++;
          if (nt==3) {
             nt=0;
             translation[aa]=codon.translate();
             aa++;
             }
          }//for each nt
        } //for each exon
    } // + strand
 translation[aa]=0;
 if (rlen!=NULL) *rlen=aa;
 return translation;
}

void GffObj::printSummary(FILE* fout) {
 if (fout==NULL) fout=stdout;
 fprintf(fout, "%s\t%c\t%d\t%d\t%4.2f\t%4.1f\n", gffID,
          strand, start, end, gscore, (float)qcov/10.0);
}

void GffObj::printGxfLine(FILE* fout, char* tlabel, char* gseqname, bool iscds,
                             uint segstart, uint segend, int exidx, char phase, bool gff3) {
  static char scorestr[14];
  strcpy(scorestr,".");
  GffAttrs* xattrs=NULL;
  if (exidx>=0) {
     if (exons[exidx]->score) sprintf(scorestr,"%.2f", exons[exidx]->score);
     xattrs=exons[exidx]->attrs;
  }
  if (phase==0 || !iscds) phase='.';
  const char* ftype=iscds ? "CDS" : getSubfName();
  if (gff3) {
    fprintf(fout,
      "%s\t%s\t%s\t%d\t%d\t%s\t%c\t%c\tParent=%s",
      gseqname, tlabel, ftype, segstart, segend, scorestr, strand,
      phase, gffID);
    if (xattrs!=NULL) {
      for (int i=0;i<xattrs->Count();i++)
         fprintf(fout, ";%s=%s",names->attrs.getName(xattrs->Get(i)->attr_id),
                           xattrs->Get(i)->attr_val);
         }
    fprintf(fout, "\n");
    } //GFF
  else {//for GTF -- we print only transcripts here
    char* geneid=(gname!=NULL)? gname : gffID;
    fprintf(fout, "%s\t%s\t%s\t%d\t%d\t%s\t%c\t%c\t",
           gseqname, tlabel, ftype, segstart, segend, scorestr, strand, phase);
    //if (isValidTranscript())
    fprintf(fout,"gene_id \"%s\"; transcript_id \"%s\";", geneid, gffID);
    bool gene_name_attr=false;
    if (xattrs!=NULL) {
          for (int i=0;i<xattrs->Count();i++) {
            if (xattrs->Get(i)->attr_val==NULL) continue;
            const char* attrname=names->attrs.getName(xattrs->Get(i)->attr_id);
            fprintf(fout, " %s ",attrname);
            if (strcmp(attrname, "gene_name")==0) gene_name_attr=true;
            if (xattrs->Get(i)->attr_val[0]=='"')
                     fprintf(fout, "%s;",xattrs->Get(i)->attr_val);
                else fprintf(fout, "\"%s\";",xattrs->Get(i)->attr_val);
             }
          }
    if (gname!=NULL && !gene_name_attr) {
       fprintf(fout, " gene_name ");
       if (gname[0]=='"') fprintf (fout, "%s;",gname);
                     else fprintf(fout, "\"%s\";",gname);
       }
    fprintf(fout, "\n");
    }//GTF
}

void GffObj::printGxf(FILE* fout, GffPrintMode gffp, char* tlabel) {
 static char tmpstr[255];
 if (tlabel==NULL) {
    tlabel=track_id>=0 ? names->tracks.Get(track_id)->name :
         (char*)"gffobj" ;
    }
 unxcoord();
 //if (exons.Count()==0) return;
 char* gseqname=names->gseqs.Get(gseq_id)->name;
 bool gff3 = (gffp>=pgffAny);
 bool showCDS = (gffp==pgtfAny || gffp==pgtfCDS || gffp==pgffCDS || gffp==pgffAny || gffp==pgffBoth);
 bool showExon = (gffp<=pgtfExon || gffp==pgffAny || gffp==pgffExon || gffp==pgffBoth);
 if (gff3) {
   //print GFF3 mRNA line:
   if (gscore>0.0) sprintf(tmpstr,"%.2f", gscore);
          else strcpy(tmpstr,".");
   uint pstart, pend;
   if (gffp==pgffCDS) {
      pstart=CDstart;
      pend=CDend;
      }
   else { pstart=start;pend=end; }
   //const char* ftype=isTranscript() ? "mRNA" : getFeatureName();
   const char* ftype=getFeatureName();
   fprintf(fout,
     "%s\t%s\t%s\t%d\t%d\t%s\t%c\t.\tID=%s",
     gseqname, tlabel, ftype, pstart, pend, tmpstr, strand, gffID);
   if (CDstart>0 && !showCDS && !isCDS) fprintf(fout,";CDS=%d-%d",CDstart,CDend);
   if (parent!=NULL && !parent->isDiscarded())
       fprintf(fout, ";Parent=%s",parent->getID());
   bool gene_name_attr=false;
   if (attrs!=NULL) {
      for (int i=0;i<attrs->Count();i++) {
        const char* attrname=names->attrs.getName(attrs->Get(i)->attr_id);
        if (strcmp("gene_name",attrname)==0) gene_name_attr=true;
        fprintf(fout,";%s=%s", attrname,
               attrs->Get(i)->attr_val);
        }
      }
    if (gname!=NULL && !gene_name_attr)
       fprintf(fout, ";gene_name=%s",gname);
    fprintf(fout,"\n");
   }// gff3 mRNA line
 if (showExon) {
   //print exons
    if (isCDS && exons.Count()>0 && 
        ((strand=='-' && exons.Last()->phase<'0') || (strand=='+' && exons.Last()->phase<'0')))
         updateExonPhase();

    for (int i=0;i<exons.Count();i++) {
      printGxfLine(fout, tlabel, gseqname, isCDS, exons[i]->start, exons[i]->end, i, exons[i]->phase, gff3);
      }
    }//printing exons
 if (showCDS && !isCDS && CDstart>0) {
    GArray<GffCDSeg> cds(true,true);
    getCDSegs(cds);
    for (int i=0;i<cds.Count();i++) {
      printGxfLine(fout, tlabel, gseqname, true, cds[i].start, cds[i].end, -1, cds[i].phase, gff3);
      }
  } //showCDS
}

void GffObj::updateExonPhase() {
  if (!isCDS) return;
  int cdsacc=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsacc+= 3-(CDphase-'0');
      }
  if (strand=='-') { //reverse strand
     for (int i=exons.Count()-1;i>=0;i--) {
         exons[i]->phase='0'+ (3-cdsacc%3)%3;
         cdsacc+=exons[i]->end-exons[i]->start+1;
         }
     }
    else { //forward strand
     for (int i=0;i<exons.Count();i++) {
         exons[i]->phase='0'+ (3-cdsacc%3)%3;
         cdsacc+=exons[i]->end-exons[i]->start+1;
         }
     }
}


void GffObj::getCDSegs(GArray<GffCDSeg>& cds) {
  GffCDSeg cdseg;
  int cdsacc=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsacc+= 3-(CDphase-'0');
      }
  if (strand=='-') {
     for (int x=exons.Count()-1;x>=0;x--) {
        uint sgstart=exons[x]->start;
        uint sgend=exons[x]->end;
        if (CDend<sgstart || CDstart>sgend) continue;
        if (CDstart>=sgstart && CDstart<=sgend)
              sgstart=CDstart; //cdstart within this segment
        if (CDend>=sgstart && CDend<=sgend)
              sgend=CDend; //cdend within this segment
        cdseg.start=sgstart;
        cdseg.end=sgend;
        cdseg.exonidx=x;
        //cdseg.phase='0'+(cdsacc>0 ? (3-cdsacc%3)%3 : 0);
        cdseg.phase='0'+ (3-cdsacc%3)%3;
        cdsacc+=sgend-sgstart+1;
        cds.Add(cdseg);
       } //for each exon
     } // - strand
    else { // + strand
     for (int x=0;x<exons.Count();x++) {
       uint sgstart=exons[x]->start;
       uint sgend=exons[x]->end;
       if (CDend<sgstart || CDstart>sgend) continue;
       if (CDstart>=sgstart && CDstart<=sgend)
             sgstart=CDstart; //seqstart within this segment
       if (CDend>=sgstart && CDend<=sgend)
             sgend=CDend; //seqend within this segment
       cdseg.start=sgstart;
       cdseg.end=sgend;
       cdseg.exonidx=x;
       //cdseg.phase='0'+(cdsacc>0 ? (3-cdsacc%3)%3 : 0);
       cdseg.phase='0' + (3-cdsacc%3)%3 ;
       cdsacc+=sgend-sgstart+1;
       cds.Add(cdseg);
       } //for each exon
   } // + strand
}
