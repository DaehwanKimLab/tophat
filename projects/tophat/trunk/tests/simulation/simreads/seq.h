#ifndef STDSEQ_H_
#define STDSEQ_H_

#include <stdio.h>

#define SEQ_MAX_NAME_LEN 255

#define INIT_SEQ(seq) (seq).s = 0; (seq).l = (seq).m = 0

#define CHAR2QUAL(c) \
	((isdigit(c))? ((c)-'0') : ((islower(c))? ((c)-'a'+10) : ((isupper(c))? ((c)-'A'+36) : 0)))
#define QUAL2CHAR(q) \
	(((q)<10)? ((q)+'0') : (((q)<36)? ((q)-10+'a') : (((q)<62)? ((q)-36+'A') : 'Z')))

typedef struct
{
	int l, m; /* length and maximum buffer size */
	unsigned char *s; /* sequence */
} seq_t;

#ifdef __cplusplus
extern "C" {
#endif

	void seq_set_block_size(int size);
	int seq_read_fasta(FILE*, seq_t*, char*, char*);
	int seq_read_fastq(FILE*, seq_t*, seq_t *, char*);

#ifdef __cplusplus
}
#endif

#endif /* STDSEQ_H_ */
