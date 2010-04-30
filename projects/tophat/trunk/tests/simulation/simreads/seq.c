#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "seq.h"

static int SEQ_BLOCK_SIZE = 512;

void seq_set_block_size(int size)
{
	SEQ_BLOCK_SIZE = size;
}

/* Read sequences from file "fp" in FASTA format. Sequence will be saved
 * in "seq", sequence ID in "locus", and comment saved in "comment",
 * provided "comment != 0". Sequence length will be returned. If -1 is
 * returned, no sequence is left in the file. */
int seq_read_fasta(FILE *fp, seq_t *seq, char *locus, char *comment)
{
	int c, l, max;
	char *p;
	
	c = 0;
	while (!feof(fp) && fgetc(fp) != '>');
	if (feof(fp)) return -1;
	p = locus;
	while (!feof(fp) && (c = fgetc(fp)) != ' ' && c != '\t' && c != '\n')
		if (c != '\r') *p++ = c;
	*p = '\0';
	if (comment) {
		p = comment;
		if (c != '\n') {
			while (!feof(fp) && ((c = fgetc(fp)) == ' ' || c == '\t'));
			if (c != '\n') {
				*p++ = c;
				while (!feof(fp) && (c = fgetc(fp)) != '\n')
					if (c != '\r') *p++ = c;
			}
		}
		*p = '\0';
	} else if (c != '\n') while (!feof(fp) && fgetc(fp) != '\n');
	l = 0; max = seq->m;
	while (!feof(fp) && (c = fgetc(fp)) != '>') {
		if (isalpha(c) || c == '-' || c == '.') {
			if (l + 1 >= max) {
				max += SEQ_BLOCK_SIZE;
				seq->s = (unsigned char*)realloc(seq->s, sizeof(char) * max);
			}
			seq->s[l++] = (unsigned char)c;
		}
	}
	if (c == '>') ungetc(c,fp);
	seq->s[l] = 0;
	seq->m = max; seq->l = l;
	return l;
}
int seq_read_fastq(FILE *fp, seq_t *seq, seq_t *qual, char *name)
{
	int c, l, max;
	char *p;
	
	c = 0;
	/* first read the sequence */
	while (!feof(fp) && fgetc(fp) != '@');
	if (feof(fp)) return -1;
	p = name;
	while (!feof(fp) && (c = fgetc(fp)) != ' ' && c != '\t' && c != '\n')
		if (c != '\r') *p++ = c;
	*p = '\0';
	if (c != '\n') while (!feof(fp) && fgetc(fp) != '\n');
	l = 0; max = seq->m;
	while (!feof(fp) && (c = fgetc(fp)) != '+') {
		if (isalpha(c) || c == '-' || c == '.') {
			if (l + 1 >= max) {
				max += SEQ_BLOCK_SIZE;
				seq->s = (unsigned char*)realloc(seq->s, sizeof(char) * max);
			}
			seq->s[l++] = (unsigned char)c;
		}
	}
	if (c == '+') ungetc(c, fp);
	seq->s[l] = 0;
	seq->m = max; seq->l = l;
	/* now come to the qualities */
	while (!feof(fp) && fgetc(fp) != '+');
	p = name;
	while (!feof(fp) && (c = fgetc(fp)) != ' ' && c != '\t' && c != '\n')
		if (c != '\r' && *p++ != c) {
			fprintf(stderr, "[seq_read_fastq] Inconsistent sequence name: %s. Continue anyway.\n", name);
			return seq->l;
		}
	l = 0; max = qual->m;
	while (!feof(fp) && l < seq->l) {
		c = fgetc(fp);
		if (c >= 33 && c < 127) {
			if (l + 1 >= max) {
				max += SEQ_BLOCK_SIZE;
				qual->s = (unsigned char*)realloc(qual->s, sizeof(char) * max);
			}
			qual->s[l++] = (unsigned char)(c - 33);
		}
	}
	c = fgetc(fp);
	if (c == '@') ungetc(c, fp);
	qual->s[l] = 0;
	qual->m = max; qual->l = l;
	if (seq->l != qual->l)
		fprintf(stderr, "[seq_read_fastq] Inconsistent length: %d!=%d. Continue anyway.", seq->l, qual->l);
	return seq->l;
}
