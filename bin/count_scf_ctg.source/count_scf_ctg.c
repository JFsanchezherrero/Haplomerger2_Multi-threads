/***************** count_scf_ctg.c ************************
Program: count_scf_ctg.c
Function: count scaffold and contig in a fasta file.
Input: fasta file
Output: table file
**********************************************************/

#ifdef NOLIB
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <getopt.h>
#include <stdarg.h>
#else
#include <mydef.h>
#include <stdio.h>
#include <regex.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <seq.h>
#include <myio.h>
#endif

#define LEN 256 

#ifdef NOLIB

#define seqNB0 4
#define seqNA0 20
#define seqNB 18
#define seqNA 28
#define seqNCODON 64

typedef
enum seqSeqType {
	seqSEQ_ERROR = 0,
	seqSEQ_DNA = 1,
	seqSEQ_RNA = 2,
	seqSEQ_PROTEIN = 3
} seqSeqType;

typedef 
enum seqCodonType {
	seqCODON_ERROR = 0,
	seqCODON_GEN = 1,
	seqCODON_MAM_MITO = 2, // mammalian mitochondria codon
	seqCODON_INV_MITO = 3, // mammalian mitochondria codon
	seqCODON_YEAST_MITO = 4, // yeast mitochondria codon
	seqCODON_OTHER=10
} seqCodonType;

typedef 
struct seqSeq {
	seqSeqType type;
	long lname;
	char* name;
	long lseq;
	char* seq;
} seqSeq;


char seqBase[seqNB+1]="TCAGUYRMKSWHBVDN?-";
short seqNBase[seqNB]={1,1,1,1, 1,2,2,2,2,2,2, 3,3,3,3, 4,4,0};
char seqEquateBase[seqNB][seqNB0+1]={"T","C","A","G", "T", "TC","AG","CA","TG","CG","TA", "TCA","TCG","CAG","TAG", "TCAG","TCAG","TCAG"};
char seqAA1[seqNA+1]="ARNDCQEGHILKMFPSTWYVBJZ*-?XU";
char seqAA3[seqNA][4]={"Ala","Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Asx", "Xle", "Glx", "***", "---", "???", "Xxx", "Unk"};
short seqBaseIndex[128];
short seqAAIndex[128];
char seqCodon[seqNB0][seqNB0][seqNB0];
static short is_set_index=0;

void mydef_perrorf(char* format, ...);
void mydef_pwarnf(char* format, ...);
long myio_sgets(char *src, long num, char *delimit, char *str);
long seq_num_fasta(FILE* pffas);
long seq_get_fasta(FILE* pffas, seqSeqType type, long n, seqSeq* seq);

#endif
static struct option lopt[] = {
	{"verbose", no_argument, 0, 'h'},
	{"in", optional_argument, 0, 'i'},
	{"ingz", required_argument, 0, 1},
	{"iz", required_argument, 0, 1},
	{"out", required_argument, 0, 'o'},
	{"outgz", required_argument, 0, 2},
	{"oz", required_argument, 0, 2},
	{0, 0, 0, 0}
};

typedef struct {
	long lname;
	char *name;
	long lseq;
	long nctg;
	long *ctg;
} Scf;

typedef struct {
	long lseq;
} Ctg;

int main(int argc, char *argv[])
{
	void usage()
	{
		printf("Usage:\n");
		printf("  %s\t-i [infile] -o [outfile]\n", argv[0]);
		printf("Option:\n");
		printf(" --iz [infile] : using gzip infile, instead of '-i',\n");
		printf("                 if no '-i' or '--iz', use standard input.\n");
		printf(" --oz [outfile] : using gzip outfile, instead of '-o',\n");
		printf("                  if no '-o' or '--oz', use standard output.\n");
		printf(" -h : print this usage.\n");
		printf(" --verbose : print this usage.\n");
		printf("Function: count scaffold and contig in a fasta file, and calculate N50 and L50.\n");
	}
	int opt_idx;
	long i, j, k, l;
	short is_in_gz=0, is_out_gz=0;
	char file[LEN]="", in[LEN]="", out[LEN]="";
	char c;
	FILE *pfr=stdin, *pfw=stdout, *pftmp=NULL;

	while (1) {
		i=getopt_long(argc, argv, "hi:o:", lopt, &opt_idx);
		if (i==-1) break;
		switch(i) {
		case 'h':
			usage();
			return 0;
		case 'i':
			strcpy(in, optarg);
			break;
		case 1:
			strcpy(in, optarg);
			is_in_gz=1;
			break;
		case 'o':
			strcpy(out, optarg);
			break;
		case 2:
			strcpy(out, optarg);
			is_out_gz=1;
			break;
		default:
			break;
		}
	}
	
	// open in and out file
	// {{{
	if (in[0]=='\0') {
		is_in_gz=0;
		pfr=tmpfile();
		do {
			c=fgetc(stdin);
			fputc(c, pfr);
		} while(c!=EOF);
		rewind(pfr);
	} else {
		if (is_in_gz) {
			sprintf(file, "gunzip -c %s", in);
			pfr=popen(file, "r");
			if (pfr==NULL) {
				mydef_perrorf("Can't open input fasta gzip file %s. (%x)", in, errno);
				return errno;
			}
		} else {
			pfr=fopen(in, "r");
			if (pfr==NULL) {
				mydef_perrorf("Can't open input fasta file %s. (%x)", in, errno);
				return errno;
			}
		}
	}
	//
	if (out[0]=='\0') {
		pfw=stdout;
		is_out_gz=0;
	} else {
		if (is_out_gz) {
			sprintf(file, "gzip -c > %s", out);
			pfw=popen(file, "w");
			if (pfw==NULL) {
				mydef_perrorf("Can't open output gzip file %s. (%x)", out, errno);
				return errno;
			}
		} else {
			pfw=fopen(out, "w");
			if (pfw==NULL) {
				mydef_perrorf("Can't open output file %s. (%x)", out, errno);
				return errno;
			}
		}
	}
	// }}}

	seqSeq seq;
	long nscf;
	long Nctg=128, nctg=0;
	Ctg *ctg=malloc(Nctg*sizeof(Ctg));
	nscf=seq_num_fasta(pfr);
	if (is_in_gz) {
		pclose(pfr);
		sprintf(file, "gunzip -c %s", in);
		pfr=popen(file, "r");
		if (pfr==NULL) {
			mydef_perrorf("Can't open input fasta gzip file %s. (%x)", in, errno);
			return errno;
		}
	}
	Scf *scf=malloc(nscf*sizeof(Scf));
	for (i=0; i<nscf; i++) {
		seq_get_fasta(pfr, seqSEQ_DNA, 1, &seq);
		scf[i].lname=seq.lname;
		scf[i].name=seq.name;
		scf[i].lseq=seq.lseq;
		l=0;
		k=0;
		do {
			while (k<seq.lseq && (seq.seq[k]=='N' || seq.seq[k]=='n')) k++;
			j=myio_sgets(seq.seq+k, seq.lseq, "Nn", NULL);
			if (nctg>=Nctg) {
				Nctg<<=1;
				ctg=realloc(ctg, Nctg*sizeof(Ctg));
			}
			ctg[nctg].lseq=j-1;
			nctg++;
			l++;
			k+=j;
			while (k<seq.lseq && (seq.seq[k]=='N' || seq.seq[k]=='n')) k++;
		} while (k<seq.lseq);
		scf[i].nctg=l;
		scf[i].ctg=malloc(l*sizeof(long));
		for (j=l; j>0; j--) {
			scf[i].ctg[l-j]=nctg-j;
		}
		free(seq.seq);
	}

	long *idx_scf=malloc(nscf*sizeof(long));
	long total_scf_len=0;
	for (i=0; i<nscf; i++) {
		idx_scf[i]=i;
		total_scf_len+=scf[i].lseq;
	}
	int cmp_scf(void *a, void *b)
	{
		long i1=*(long*)a, i2=*(long*)b;
		if (scf[i1].lseq>scf[i2].lseq) return 1;
		else if (scf[i1].lseq==scf[i2].lseq) return 0;
		else return -1;
	//	return scf[i1].lseq-scf[i2].lseq;
	}
	qsort(idx_scf, nscf, sizeof(long), cmp_scf);
	//
	long *idx_ctg=malloc(nctg*sizeof(long));
	long total_ctg_len=0;
	for (i=0; i<nctg; i++) {
		idx_ctg[i]=i;
		total_ctg_len+=ctg[i].lseq;
	}
	int cmp_ctg(void *a, void *b)
	{
		long i1=*(long*)a, i2=*(long*)b;
		if (ctg[i1].lseq>ctg[i2].lseq) return 1;
		else if (ctg[i1].lseq==ctg[i2].lseq) return 0;
		else return -1;
	//	return ctg[i1].lseq-ctg[i2].lseq;
	}
	qsort(idx_ctg, nctg, sizeof(long), cmp_ctg);

	fprintf(pfw, "# scaffolds\n");
	float prob[21], nprob=20, p=0;
	for (i=1; i<=nprob; i++) {prob[i-1]=i*0.05;}
	fprintf(pfw, "# total bp size: %li\n", total_scf_len);
	fprintf(pfw, "# order by scaffold-size (descending)\n");
	fprintf(pfw, "# \%\t(accumulated_size)\tscaffold_number\tcritical_size\n");
	j=0;
	k=prob[0]*total_scf_len;
	l=0;
	for (i=0; i<nscf; i++) {
		l+=scf[idx_scf[i]].lseq;
		while (j<nprob-1 && l>k) {
			fprintf(pfw, "%.2f\t(%li)\t%li\t%li\n", prob[j], l, i, scf[idx_scf[i]].lseq);
			j++; 
			k=total_scf_len*prob[j];
		}
	}
	fprintf(pfw, "%.2f\t(%li)\t%li\t%li\n", prob[j], l, i, scf[idx_scf[nscf-1]].lseq);

	fprintf(pfw, "# contigs\n");
	fprintf(pfw, "# total bp size: %li\n", total_ctg_len);
	fprintf(pfw, "# order by contig-size (descending)\n");
	fprintf(pfw, "# \%\t(accumulated_size)\tcontig_number\tcritical_size\n");
	j=0;
	k=prob[0]*total_ctg_len;
	l=0;
	for (i=0; i<nctg; i++) {
		l+=ctg[idx_ctg[i]].lseq;
		while (j<nprob-1 && l>k) {
			fprintf(pfw, "%.2f\t(%li)\t%li\t%li\n", prob[j], l, i, ctg[idx_ctg[i]].lseq);
			j++; 
			k=prob[j]*total_ctg_len;
		}
	}
	fprintf(pfw, "%.2f\t(%li)\t%li\t%li\n", prob[j], l, i, ctg[idx_ctg[nctg-1]].lseq);

	if (is_in_gz) {
		pclose(pfr);
	} else {
		fclose(pfr);
	}
	if (is_out_gz) {
		pclose(pfw);
	} else {
		fclose(pfw);
	}
}

#ifdef NOLIB
void mydef_perrorf(char* format, ...)
// {{{
{
	char err_msg0[256];
	va_list args;
	va_start(args, format);
	vsprintf(err_msg0, format, args);
	va_end(args);
	printf("Error: %s\n", err_msg0);
}
// }}}

long myio_sgets(char *src, long num, char *delimit, char *str)
// {{{
{
	long i=0, j;
	char c;
	while (i<num && src[i]!='\0') {
		c=src[i];
		j=0;
		while (delimit[j]!='\0' && c!=delimit[j]) j++;
		if (str!=NULL) str[i]=c;
		i++;
		if (delimit[j]!='\0') {
			break;
		}
	}
	if (str!=NULL) str[i]='\0';
	return i;
}
// }}}

long seq_get_fasta(FILE* pffas, seqSeqType type, long n, seqSeq* seq)
// {{{1
{	
	{
	short i;
	char c;
	for (i=0; i<128; i++) seqBaseIndex[i]=-1;
	i=0;
	while(seqBase[i]!='\0') {
		c=seqBase[i];
		if (isalpha(c)) {
			seqBaseIndex[c|0x20]=i;
			seqBaseIndex[c&~0x20]=i;
		} else {
			seqBaseIndex[c]=i;
		}
		i++;
	}

	i=0;
	for (i=0; i<128; i++) seqAAIndex[i]=-1;
	while(seqAA1[i]!='\0') {
		c=seqAA1[i];
		seqAAIndex[c]=i;
		if (isalpha(c)) {
			seqAAIndex[c|0x20]=i;
		}
		i++;
	}
	}

	char c, *p;
	long i, j;
	long L=512, dL=L, l;

	short *tag=NULL;
	j=0;
	if (type==seqSEQ_DNA || type==seqSEQ_RNA ) {
		tag=seqBaseIndex;
	} else if (type==seqSEQ_PROTEIN) {
		tag=seqAAIndex;
	}

	c=fgetc(pffas);
	while (c!=EOF && c!='>') c=fgetc(pffas);	
	for (i=0; i<n; i++) {
		if (c==EOF) {break;}
		if (c=='>') {
			/* get name */
			l=0;
			p=(char*)malloc((L+1)*sizeof(char));
			if (p==NULL) return -1;
			c=fgetc(pffas);
			while (c!=EOF && c!='\n' && c!='\r') {
				if (l>=L) {
					L+=dL;
					dL<<=1;
					p=(char*)realloc(p, (L+1)*sizeof(char));
					if (p==NULL) return -1;
				}
				p[l++]=c;
				c=fgetc(pffas);
			}
			c=fgetc(pffas);
			while (isspace(c)) c=fgetc(pffas);
			if (c==EOF) {free(p); break;}
			if (l<L) p=(char*)realloc(p, l+1);
			p[l]='\0';
			(seq+i)->lname=l;
			(seq+i)->name=p;

			// get sequence 
			L=1024; dL=L;
			l=0;
			p=(char*)malloc((L+1)*sizeof(char));
			if (p==NULL) {
				mydef_perrorf("can't malloc for sequence.");
				return -1;
			}
			while(c!=EOF && c!='>') {
				if (tag!=NULL && tag[c]>=0) {
					if (l>=L) {
						L+=dL;
						dL<<=1;
						p=(char*)realloc(p, (L+1)*sizeof(char));
						if (p==NULL) {free((seq+i)->name); return -1;}
					}
					if (type==seqSEQ_DNA) {
						if (c=='u') c='t';
						else if (c=='U') c='T';
					} else if (type==seqSEQ_RNA) {
						if (c=='t') c='u';
						else if (c=='T') c='U';
					}
					p[l++]=c;
				} else if (!isspace(c)) {
					mydef_pwarnf("%c is invalid character", c);
				}
				c=fgetc(pffas);
			}
			if (l<L) p=(char*)realloc(p, (l+1)*sizeof(char));
			p[l]='\0';
			(seq+i)->lseq=l;
			(seq+i)->seq=p;
			(seq+i)->type=type;
		}
	}
	if (c=='>') ungetc(c, pffas);
	return i;
}
// 1}}}

long seq_num_fasta(FILE* pffas)
// {{{2
{
	fpos_t oldfpos;
	fgetpos(pffas, &oldfpos);
	long n=0;
	char c;
	c=fgetc(pffas);
	while(c!=EOF && c!='>') c=fgetc(pffas);
	while(c=='>') {
		n++;
		c=fgetc(pffas);
		while(c!=EOF && c!='>') c=fgetc(pffas);
	}
	clearerr(pffas);
	fsetpos(pffas, &oldfpos);
	return n;
}
// 2}}}

void mydef_pwarnf(char* format, ...)
// {{{
{
	fprintf(stderr, "Warnning: ");
	va_list args;
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	fprintf(stderr, "\n");
}
// }}}
#endif
