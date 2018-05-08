/*  bali_score.c  - version 3.01 					*/
/*                   							*/
/*  Changes since version 3.0						*/
/*  1. Bug fixed in full length alignment scores using MSF reference file  */
/*                   							*/
/*  Changes since version 2.0						*/
/*  1. The BAliBASE alignments are now stored in XML format. The expat  */
/*  XML parser has therefore been incorporated in version 3.0. The 	*/
/*  expat parser is available from http://expat.sourceforge.net		*/
/*                   							*/
/*  Changes since version 1.0						*/
/*  1. Bug fixed that scored less than 1, for columns containing gaps	*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include "score.h"

#define XML 1
#define GCG 2

int SeqGCGCheckSum(char *seq, int len);

ALN read_msf(FILE *fin,int nseqs);
int score_ref(int *refseq_col,int maxlen);
int checkref(FILE *fin);
int countmsf(FILE *fin);
void columnscore(ALNPTR test_aln,int nseqs,int seqlength,int *refseq_col,int **seq_code,Boolean verbose,float *pc_res,int *tc_score);
void code_refseq(ALNPTR mult_aln,int maxlen,int **refseq_code);
void code_seq(ALNPTR ref_aln,ALNPTR test_aln,int *seq_xref,int *refseq_col,int **refseq_code,int **seq_code,int refseqlength,int seqlength);
void aln_seq(int i);
void ref_gaps(ALNPTR ref_aln,int maxlen,int cutoff,int *refseq_col);
int get_coreblocks(ALNPTR mult_aln,int *refseq_col);

Boolean verbose;

int refscore;
/* cutoff for number of gaps allowed in a column in the reference alignment 
for column to be included in alignment score */
int cutoff;

/* core block data read in from annotation file */
int nblocks;
typedef struct Block {
	int s;
	int e;
	int nseqs;
} Block;

Block blocks[100];


int main(int argc, char **argv)
{
	FILE *ifd,*tfd,*afd;
	int  err,i,j,ires,iseq=0;
	int ix,n;
	int format;
	int nseqs,refnseqs;
	int maxlen,refmaxlen;
	char t[MAXLINE+1];
	char seq[MAXLINE+1];
	char clen[MAXLINE+1];
	char cvar[MAXLINE+1];
	char cprog[MAXLINE+1];
	char method;
	Boolean eof,found;
	int *refseq_col;
	int *seq_xref;
	int **seq_code,**refseq_code;
	int tc_score;
	float maxpc_res;
	float pc_res;
	ALN ref_aln;
	ALN test_aln;

	if(argc<3) {
		fprintf(stderr,"Usage: %s ref_aln test_aln [-v]\n",argv[0]);
		fprintf(stderr,"                where ref_aln       reference alignment in xml/msf format \n");
		fprintf(stderr,"                      test_aln      test alignment in msf format \n");
		fprintf(stderr,"                      -v            verbose mode\n");
		return 1;
	}
	if(argc==4) verbose=TRUE;
	else verbose=FALSE;

/* open the reference aln file */

        if((ifd=fopen(argv[1],"r"))==NULL) {
            fprintf(stderr,"Cannot open reference aln file [%s]",argv[1]);
            return 1;
        }
/* open the test aln file */

        if((tfd=fopen(argv[2],"r"))==NULL) {
            fprintf(stderr,"Cannot open test aln file [%s]",argv[2]);
            return 1;
        }

/* set method="C" if total columns score, or 'B' for just core blocks */
	method='C';

/* read the reference alignment into refaln */
	format=checkref(ifd);

	if(format==GCG) {
		refnseqs=countmsf(ifd);
	}
	else if(format==XML) {
		refnseqs=count_xml_seqs(ifd);
	}

        fseek(ifd,0,0);
	if(refnseqs==0) {
		fprintf(stderr,"Error: no sequences in %s\n",argv[1]);
		return 1;
	}
	if(format==GCG) {
		ref_aln=read_msf(ifd,refnseqs);
	}
	else if(format==XML) {
		ref_aln=read_xml(ifd,0);
		ref_aln.nseqs=refnseqs;
		method='B';
	}
	refmaxlen=0;
	for(i=0;i<ref_aln.nseqs;i++)
		if(refmaxlen<ref_aln.seqs[i].len) refmaxlen=ref_aln.seqs[i].len;

/* read the test alignment into names, seq_array, seqlength */
	nseqs = countmsf(tfd);
	if(nseqs==0) {
		fprintf(stderr,"Error: no sequences in %s\n",argv[2]);
		return 1;
	}
	test_aln=read_msf(tfd,nseqs);

	if(nseqs != refnseqs) {
		fprintf(stderr,"Error: %d sequences in %s and %d in %s",refnseqs,argv[1],nseqs,argv[2]);
		return 1;
	}

	maxlen=0;
	for(i=0;i<test_aln.nseqs;i++)
		if(maxlen<test_aln.seqs[i].len) maxlen=test_aln.seqs[i].len;



/* cross-reference the sequence refnames, in case they're not in the same order
in the reference and test alignment files */
	seq_xref=(int *)ckalloc((refnseqs+2)*sizeof(int));
	for(i=0;i<refnseqs;i++)
		seq_xref[i]=-1;
	for(i=0;i<refnseqs;i++) {
		found=FALSE;
		for(j=0;j<test_aln.nseqs;j++) {
			if(strcasecmp(test_aln.seqs[j].name,ref_aln.seqs[i].name)==0)
			{
				found=TRUE;
				seq_xref[j]=i;
				break;
			}
		}
		if(found==FALSE) {
			fprintf(stderr,"Error: sequence %s not found in test aln %s\n",ref_aln.seqs[i].name,argv[2]);
			return 1;
		}
	}

        fprintf(stdout,"\nComparing test alignment in %s\nwith reference alignment in %s\n",argv[2],argv[1]);

/* get the core blocks */
	refseq_col=(int *)ckalloc((refmaxlen+2)*sizeof(int));
        if(method=='B') {
                n=get_coreblocks(&ref_aln,refseq_col);
		if(n<0) method='C';
        }

/* if we're using all columns, only consider columns with less than 'cutoff' 
gaps. Here we use 20% of the number of sequences */
	if(method=='C')
	{
        	cutoff=(float)refnseqs*20.0/100.0;
        	if(cutoff<1) cutoff=1;
        	ref_gaps(&ref_aln,refmaxlen,cutoff,refseq_col);
	}
	if(method=='B') fprintf(stdout,"\nUsing core blocks defined in %s\n",argv[1]);

/* code the reference alignment - assign to each residue the number of the column it's in 
   gap positions are coded 0 */

	refseq_code=(int **)ckalloc((ref_aln.nseqs+2)*sizeof(int *));
	for(i=0;i<ref_aln.nseqs;i++)
		refseq_code[i]=(int *)ckalloc((refmaxlen+2)*sizeof(int));
	code_refseq(&ref_aln,refmaxlen,refseq_code);


/* calculate the max score possible ie the score for the reference alignment */
	maxpc_res=score_ref(refseq_col,refmaxlen);
	if(maxpc_res<=0) {
		fprintf(stdout,"Error in reference alignment\n");
		exit(1);
	}

/* code the test alignment - look up each residue from the test alignment in the reference
alignment and assign the reference column number */
	seq_code=(int **)ckalloc((test_aln.nseqs+2)*sizeof(int *));
	for(i=0;i<test_aln.nseqs;i++)
		seq_code[i]=(int *)ckalloc((maxlen+2)*sizeof(int));
	code_seq(&ref_aln,&test_aln,seq_xref,refseq_col,refseq_code,seq_code,refmaxlen,maxlen);

/* calculate the scores */
	columnscore(&test_aln,test_aln.nseqs,maxlen,refseq_col,seq_code,verbose,&pc_res,&tc_score);
	pc_res/=maxpc_res;

        fprintf(stdout,"\n\tSP score= %.3f\n",pc_res);
        fprintf(stdout,"\n\tTC score= %.3f\n",(float)tc_score/100.0);
fprintf(stdout,"auto %s %.3f %.3f\n",argv[2],pc_res,(float)tc_score/100.0);

	exit(0);

	
}

int get_coreblocks(ALNPTR mult_aln,int *refseq_col)
{
	int i,n;

	n=(-1);
	for(i=0;i<mult_aln->ncol_scores;i++) {
		if(strcmp(mult_aln->col_score[i].name,"coreblock")==0) {
			n=i;
			break;
		}
	}
	if(n<0) {
		fprintf(stdout,"Warning: no core blocks in XML reference alignment, using all columns\n");
		return n;
	}
	
	for(i=0; i<mult_aln->col_score[n].length;i++) {
		if(mult_aln->col_score[n].data[i]==1) 
			refseq_col[i]=mult_aln->nseqs;
		else refseq_col[i]=0;
	}
	return n;
}

void code_refseq(ALNPTR mult_aln,int maxlen,int **refseq_code)
{
	int seq,i;

/* assign column no. in reference alignment to each residue */

	for(seq=0;seq<mult_aln->nseqs;seq++) {
		for(i=0;i<maxlen;i++) {
			if(i>=mult_aln->seqs[seq].len) {
				refseq_code[seq][i]=0;
			}
               		else if(mult_aln->seqs[seq].data[i]=='-') {
				refseq_code[seq][i]=0;
			}
			else {
				refseq_code[seq][i]=i+1;
			}
		}
	}
}

void ref_gaps(ALNPTR ref_aln,int maxlen,int cutoff,int *refseq_col)
{
	int i,j;
	int *gap;			 /* number of gaps in a column in the reference alignment */

/* find columns with gaps in the reference sequnce - set refseq_col[]=0
if gaps, =nseqs otherwise */

	gap=(int *)ckalloc((maxlen+2)*sizeof(int));
	for(i=0;i<maxlen;i++)
	{
                gap[i]=0;
                for(j=0;j<ref_aln->nseqs;j++)
                        if(ref_aln->seqs[j].data[i]=='-')
                        {
                                gap[i]++;
                        }
	}

	for(i=0;i<maxlen;i++)
	{
		if(gap[i]>=cutoff) refseq_col[i]=0;
		else refseq_col[i]=ref_aln->nseqs-gap[i];
	}
	ckfree(gap);
}

int score_ref(int *refseq_col,int maxlen)
{
	int i,j;
	int maxpc_res;

/* calculates the maximum score possible for an alignment */
	maxpc_res=0;	
	for(i=0;i<maxlen;i++)
	{
		if(refseq_col[i]>1) 
		{
			maxpc_res+=refseq_col[i]*(refseq_col[i]-1)/2.0;
		}
	}
	return maxpc_res;
}

void code_seq(ALNPTR ref_aln,ALNPTR test_aln,int *seq_xref,int *refseq_col,int **refseq_code,int **seq_code,int refseqlength,int seqlength)
{
	int i,j,ix,seq;

	for(seq=0;seq<test_aln->nseqs;seq++) {
		if ( seq_xref[seq] == -1) continue;
/* find the first residue in the reference sequence */
		ix=0;
		for(j=0;j<refseqlength;j++)
			if(ref_aln->seqs[seq_xref[seq]].data[j]!='-') {
				ix=refseq_code[seq_xref[seq]][j];
				break;
			}

		for(i=0;i<seqlength;i++) {
			if(test_aln->seqs[seq].data[i]=='-') {
				seq_code[seq_xref[seq]][i]=0;
			}
			else {
				if (refseq_col[ix-1]>0 && ref_aln->seqs[seq_xref[seq]].data[ix-1]!='-') seq_code[seq_xref[seq]][i]=ix;
				for(j+=1;j<refseqlength;j++)
					if(ref_aln->seqs[seq_xref[seq]].data[j]!='-') {
						ix=refseq_code[seq_xref[seq]][j];
						break;
					}
			}
		}
	}
}

void columnscore(ALNPTR test_aln,int nseqs,int seqlength,int *refseq_col,int **seq_code,Boolean verbose,float *pc_res,int *tc_score)
{
	int i,j,k;
	int iseq,ncols;
	int *scores;
	int *index;
	int n;
	int *colscore1;
	int pc;
	int sp,tc;
	int nblocks;
	int blen=50;
	Boolean found;

	scores=(int *)ckalloc((nseqs+1)*sizeof(int));
	index=(int *)ckalloc((nseqs+1)*sizeof(int));
	colscore1=(int *)ckalloc((seqlength+1)*sizeof(int));

	tc=ncols=0;
	sp=0;
	for(j=0;j<seqlength;j++)
		colscore1[j]=0;
	

	for(i=0;i<seqlength;i++)
	{
		for(j=0;j<nseqs;j++)
			scores[j]=0;
		n=0;
		for(j=0;j<nseqs;j++)
		{
			if(seq_code[j][i]!=0)
			{
				found=FALSE;
				for(k=0;k<n;k++)
				{
					if(seq_code[j][i]==index[k])
					{
						scores[k]++;
						found=TRUE;
						break;
					}
				}
				if(found==FALSE)
				{
					scores[n]=1;
					index[n]=seq_code[j][i];
					n++;
				}
			}
		}
		for(j=0;j<nseqs;j++)
		{
			if(scores[j]>1) {
			/*	if(scores[j]>refseq_col[seq_code[0][i]-1]) scores[j]=refseq_col[seq_code[0][i]-1];*/
				sp+=scores[j]*(scores[j]-1)/2.0;
			}
		}
/* count 1 for each column */
		if(seq_code[0][i]>0 && scores[0]>=refseq_col[seq_code[0][i]-1])
			colscore1[i]=1;
		if (seq_code[0][i]!=0) ncols++;

		tc+=colscore1[i];
	}

	if(verbose)
	{
		nblocks=seqlength/blen;
                fprintf(stdout,"\n\n");
		for(k=0;k<=nblocks;k++) {
        		for(i=0;i<nseqs;i++) {
                		fprintf(stdout,"%20s   ",test_aln->seqs[i].name);
               			for(j=0;j<blen && k*blen+j<seqlength;j++) {
                       			fprintf(stdout,"%c",test_aln->seqs[i].data[k*blen+j]);
                		}
                		fprintf(stdout,"\n");
			}
                	fprintf(stdout,"\n");
                	fprintf(stdout,"\n                       ");
               		for(j=0;j<blen && k*blen+j<seqlength;j++) {
				if(seq_code[0][k*blen+j]!=0)
                       			fprintf(stdout,"%d",colscore1[k*blen+j]);
				else
                       			fprintf(stdout,".");
                	}
                	fprintf(stdout,"\n\n");
        	}

	}

/* count 1 for each column */
	if(ncols>0) tc=100*(float)tc/(float)ncols;

	(*pc_res)=sp;
	(*tc_score)=tc;
}


ALN read_msf(FILE *fin,int nseqs)
{
        static char line[MAXLINE+1];
	char name[MAXNAMES+1];
        char *seq = NULL;
        int len,seqno,i,j,k;
        unsigned char c;
	ALN mult_aln;

	alloc_aln(nseqs,&mult_aln);
	mult_aln.nseqs=nseqs;

	for(seqno=0;seqno<nseqs;seqno++) {

        	fseek(fin,0,0);                 /* start at the beginning */

        	len=0;                         /* initialise length to zero */
        	for(i=0;;i++) {
                	if(fgets(line,MAXLINE+1,fin)==NULL) return mult_aln; /* read the title*/
                	if(line[0]=='/' && line[1]=='/') break;
        	}

        	while (fgets(line,MAXLINE+1,fin) != NULL) {
                	if(!blankline(line)) {

                        	for(i=0;i<seqno;i++) fgets(line,MAXLINE+1,fin);
                        	for(j=0;j<=strlen(line);j++) if(line[j] != ' ') break;
                        	for(k=j;k<=strlen(line);k++) if(line[k] == ' ') break;
                        	strncpy(name,line+j,k-j);
                        	name[k-j]='\0';
                        	name[MAXNAMES]='\0';
	
				if(seq==NULL)
					seq=(char *)ckalloc((MAXLINE+2)*sizeof(char));
				else
					seq=(char *)ckrealloc(seq,(len+MAXLINE+2)*sizeof(char));

                        	for(i=k;i<=MAXLINE;i++) {
                                	c=line[i];
                                	if(c == '.' || c == '~' ) c = '-';
                                	if(c == '*') c = 'X';
                                	if(c == '\n' || c == EOS) break; /* EOL */
                                	if(isalpha(c) || c=='-') seq[len++]=c;
                        	}
	
                        	for(i=0;;i++) {
                                	if(fgets(line,MAXLINE+1,fin)==NULL) break;
                                	if(blankline(line)) break;
                        	}
                	}
        	}
		alloc_seq(&mult_aln.seqs[seqno],len);
                strcpy(mult_aln.seqs[seqno].name,name);
		strcpy(mult_aln.seqs[seqno].data,seq);
		mult_aln.seqs[seqno].len=len;
		seq=ckfree(seq);
	}

	return mult_aln;
}

int checkref(FILE *fin)
{

        char line[MAXLINE+1];

	int format;

	format=0;

        while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(strncmp(line,"<?xml",5)==0) {
			fprintf(stdout,"Using XML format reference alignment\n");
			format=XML;
			break;
		}
                if(line[0]=='/' && line[1]=='/') {
			fprintf(stdout,"Using GCG format reference alignment\n");
			format=GCG;
			break;
		}
        }
        fseek(fin,0,0);

	return format;
}

int countmsf(FILE *fin)
{
/* count the number of sequences in a PILEUP alignment file */

        char line[MAXLINE+1];
        int  lnseqs;

        while (fgets(line,MAXLINE+1,fin) != NULL) {
                if(line[0]=='/' && line[1]=='/') break;
        }

        while (fgets(line,MAXLINE+1,fin) != NULL) {
                if(!blankline(line)) break;             /* Look for next non- */
        }                                               /* blank line */
        lnseqs = 1;

        while (fgets(line,MAXLINE+1,fin) != NULL) {
                if(blankline(line)) return lnseqs;
                lnseqs++;
        }

        return 0; /* if you got to here-funny format/no seqs.*/
}

void fatal( char *msg,...)
{
        va_list ap;

        va_start(ap,msg);
        fprintf(stdout,"\n\nFATAL ERROR: ");
        vfprintf(stdout,msg,ap);
        fprintf(stdout,"\n\n");
        va_end(ap);
        exit(1);
}

void error( char *msg,...)
{
        va_list ap;

        va_start(ap,msg);
        fprintf(stdout,"\n\nERROR: ");
        vfprintf(stdout,msg,ap);
        fprintf(stdout,"\n\n");
        va_end(ap);
}

void warning( char *msg,...)
{
        va_list ap;

        va_start(ap,msg);
        fprintf(stdout,"\n\nWARNING: ");
        vfprintf(stdout,msg,ap);
        fprintf(stdout,"\n\n");
        va_end(ap);
}

void info( char *msg,...)
{
        va_list ap;

        va_start(ap,msg);
        fprintf(stdout,"\n");
        vfprintf(stdout,msg,ap);
        va_end(ap);
}

