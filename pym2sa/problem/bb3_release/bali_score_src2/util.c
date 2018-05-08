#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include "clustalw.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

static char revision_level[10]="";
static char help_file_name[FILENAMELEN]="";
static Boolean usemenu=FALSE;

void set_usemenu(Boolean flag)
{
	usemenu=flag;
}

Boolean get_usemenu(void)
{
	return usemenu;
}


void set_revision_level(char *rev)
{
	strcpy(revision_level,rev);
}

char * get_revision_level(void)
{
	char *rev;

	rev=(char *)ckalloc((strlen(revision_level)+2)*sizeof(char));
	strcpy(rev,revision_level);

	return rev;
}

void set_help_file_name(char *filename)
{
	strcpy(help_file_name,filename);
}

char * get_help_file_name(void)
{
        char *file_name;

        file_name=ckalloc((strlen(help_file_name)+2)*sizeof(char));
        strcpy(file_name,help_file_name);

        return file_name;
}

/*
*	ckalloc()
*
*	Tries to allocate "bytes" bytes of memory. Exits program if failed.
*	Return value:
*		Generic pointer to the newly allocated memory.
*/

void *ckalloc(size_t bytes)
{
	register void *ret;
	
	if( (ret = calloc(bytes, sizeof(char))) == NULL)
		fatal("Out of memory\n");
	else
		return ret;	

	return ret;	
}

/*
*	ckrealloc()
*
*	Tries to reallocate "bytes" bytes of memory. Exits program if failed.
*	Return value:
*		Generic pointer to the re-allocated memory.
*/

void *ckrealloc(void *ptr, size_t bytes)
{
	register void *ret=NULL;

	if (ptr == NULL)	
		fatal("Bad call to ckrealloc\n");
	else if( (ret = realloc(ptr, bytes)) == NULL)
		fatal("Out of memory\n");
	else
		return ret;	

	return ret;	
}

/*
*	ckfree()
*
*	Tries to free memory allocated by ckalloc.
*	Return value:
*		None.
*/

void *ckfree(void *ptr)
{
	if (ptr == NULL)
		warning("Bad call to ckfree\n");
	else {
	 	free(ptr);
		ptr = NULL;
	}
	return ptr;
}


/*
*	rtrim()
*
*	Removes trailing blanks from a string
*
*	Return values:
*		Pointer to the processed string
*/

char * rtrim(char *str)
{
	register int p;

	p = strlen(str) - 1;
	
	while ( isspace(str[p]) )
		p--;
		
	str[p + 1] = EOS;
	
	return str;
}


/*
*	blank_to_()
*
*	Replace blanks in a string with underscores
*
*       Also replaces , ; : ( or ) with _
*
*	Return value:
*		Pointer to the processed string
*/

char * blank_to_(char *str)
{
	int i,p;

	p = strlen(str);
	
	for(i=0;i<p;i++) 
		if(
                     (str[i]==' ') ||
                     (str[i]==';') ||
                     (str[i]==',') ||
                     (str[i]=='(') ||
                     (str[i]==')') ||
                     (str[i]==':')
                  )
                      str[i] = '_';
	
	return str;
}


/*
*	upstr()
*
*	Converts string str to uppercase.
*	Return values:
*		Pointer to the converted string.
*/

char * upstr(char *str)
{
	register char *s = str;
	
	while( (*s = toupper(*s)) )
		s++;
		
	return str;
}

/*
*	lowstr()
*
*	Converts string str to lower case.
*	Return values:
*		Pointer to the converted string.
*/

char * lowstr(char *str)
{
	sint i,len;

	len=strlen(str);
	for(i=0;i<len;i++) {
		if(isupper(str[i])) str[i] = tolower(str[i]);
	}
		
	return str;
}

void getstr(char *instr,char *outstr)
{	
	fprintf(stdout,"%s: ",instr);
	gets(outstr);
}

double getreal(char *instr,double minx,double maxx,double def)
{
	int status;
	float ret;
	char line[MAXLINE];	
	
	while(TRUE) {
		fprintf(stdout,"%s (%.1f-%.1f)   [%.1f]: ",instr,minx,maxx,def);
		gets(line);
		status=sscanf(line,"%f",&ret);
		if(status == EOF) return def;
		if(ret>maxx) {
			fprintf(stdout,"ERROR: Max. value=%.1f\n\n",maxx);
			continue;
		}
		if(ret<minx) {
			fprintf(stdout,"ERROR: Min. value=%.1f\n\n",minx);
			continue;
		}
		break;
	}
	return (double)ret;
}


int getint(char *instr,int minx,int maxx, int def)
{
	int ret,status;
	char line[MAXLINE];	

	while(TRUE) {
		fprintf(stdout,"%s (%d..%d)    [%d]: ",
		instr,(pint)minx,(pint)maxx,(pint)def);
		gets(line);
		status=sscanf(line,"%d",&ret);
		if(status == EOF) return def;
		if(ret>maxx) {
			fprintf(stdout,"ERROR: Max. value=%d\n\n",(pint)maxx);
			continue;
		}
		if(ret<minx) {
			fprintf(stdout,"ERROR: Min. value=%d\n\n",(pint)minx);
			continue;
		}
		break;
	}
	return ret;
}

void do_system(void)
{
	char line[MAXLINE];
	
	getstr("\n\nEnter system command",line);
	if(*line != EOS)
		system(line);
	fprintf(stdout,"\n\n");
}


Boolean linetype(char *line,char *code)
{
	return( strncmp(line,code,strlen(code)) == 0 );
}

Boolean keyword(char *line,char *code)
{
	int i,j;
	char key[MAXLINE];

	for(i=0;isspace(line[i]) && line[i]!=EOS;i++);
	for(j=0;!isspace(line[i]) && line[i]!=EOS;i++)
		key[j++]=line[i];
	key[j]=EOS;
	return( strcmp(key,code) == 0 );
}

Boolean blankline(char *line)
{
	int i;

	for(i=0;line[i]!='\n' && line[i]!=EOS;i++) {
		if( isdigit(line[i]) ||
		    isspace(line[i]) ||
		    (line[i] == '*') ||
		    (line[i] == ':') ||
                    (line[i] == '.')) 
			;
		else
			return FALSE;
	}
	return TRUE;
}


void get_path(char *str,char *path)
{
	register int i;
	
	strcpy(path,str);
	for(i=strlen(path)-1;i>-1;--i) {
		if(str[i]==DIRDELIM) {
			i = -1;
			break;
		}
		if(str[i]=='.') break;
	}
	if(i<0)
		strcat(path,".");
	else
		path[i+1]=EOS;
}

void alloc_aln(sint nseqs,ALNPTR mult_aln)
{
	sint i,j;

	mult_aln->seqs = (SEQ *)ckalloc((nseqs+1) * sizeof (SEQ));
	mult_aln->ft = (FT *)ckalloc((nseqs+1) * sizeof (FT));
	mult_aln->repeat = (REP *)ckalloc((nseqs+1) * sizeof (REP));
	mult_aln->go = (GO *)ckalloc((nseqs+1) * sizeof (GO));
	for(i=0;i<nseqs;i++) {
		mult_aln->seqs[i].name = (char *)ckalloc((MAXNAMES+1) * sizeof (char));
		mult_aln->seqs[i].access = (char *)ckalloc((MAXNAMES+1) * sizeof (char));
		mult_aln->seqs[i].nid = (char *)ckalloc((MAXNAMES+1) * sizeof (char));
		mult_aln->seqs[i].title = (char *)ckalloc((MAXTITLES+1) * sizeof (char));
		mult_aln->seqs[i].org = (char *)ckalloc((MAXORGANISMS+1) * sizeof (char));
                mult_aln->seqs[i].data=NULL;
                mult_aln->seqs[i].mask=NULL;
		mult_aln->seqs[i].weight=100;
		mult_aln->seqs[i].len=0;
		mult_aln->seqs[i].simgroup=0;
		for(j=0;j<MAXFTTYPE;j++)
			mult_aln->ft[i].nentries[j]=0;
		mult_aln->repeat[i].nrepeats=0;
		mult_aln->go[i].ngorefs=0;
	}

	mult_aln->motifs = NULL;
	mult_aln->groups.ngroups = 0;

        mult_aln->nseqs=0;
        mult_aln->ncol_scores=0;
        mult_aln->nanchors=0;
        mult_aln->dnaflag=FALSE;
        strcpy(mult_aln->alphabet,"ABCDEFGHIJKLMNOPQRSTUVWXYZ");
        mult_aln->prf1.nseqs=0;
        mult_aln->prf2.nseqs=0;
        mult_aln->treename[0]='\0';
        mult_aln->prf1.treename[0]='\0';
        mult_aln->prf2.treename[0]='\0';

}

void realloc_aln(sint first_seq,sint nseqs,ALNPTR mult_aln)
{
	sint i,j;

	mult_aln->seqs = (SEQ *)ckrealloc(mult_aln->seqs,(first_seq+nseqs+1) * sizeof (SEQ));
	mult_aln->ft = (FT *)ckrealloc(mult_aln->ft,(first_seq+nseqs+1) * sizeof (FT));
	mult_aln->repeat = (REP *)ckrealloc(mult_aln->repeat,(first_seq+nseqs+1) * sizeof (REP));
	mult_aln->go = (GO *)ckrealloc(mult_aln->go,(first_seq+nseqs+1) * sizeof (GO));
	for(i=first_seq;i<first_seq+nseqs;i++) {
                mult_aln->seqs[i].name = (char *)ckalloc((MAXNAMES+1) * sizeof (char));
		mult_aln->seqs[i].access = (char *)ckalloc((MAXNAMES+1) * sizeof (char));
		mult_aln->seqs[i].nid = (char *)ckalloc((MAXNAMES+1) * sizeof (char));
                mult_aln->seqs[i].title = (char *)ckalloc((MAXTITLES+1) * sizeof (char));
                mult_aln->seqs[i].org = (char *)ckalloc((MAXORGANISMS+1) * sizeof (char));
                mult_aln->seqs[i].data=NULL;
                mult_aln->seqs[i].mask=NULL;
                mult_aln->seqs[i].weight=100;
                mult_aln->seqs[i].len=0;
                mult_aln->seqs[i].simgroup=0;
		for(j=0;j<MAXFTTYPE;j++)
			mult_aln->ft[i].nentries[j]=0;
		mult_aln->repeat[i].nrepeats=0;
		mult_aln->go[i].ngorefs=0;
        }
}

void free_aln(ALNPTR mult_aln)
{
	sint i,j,k;

	if(mult_aln->nseqs<=0) return;

        for(i=0;i<mult_aln->nseqs;i++) {
		ckfree(mult_aln->seqs[i].name);
		ckfree(mult_aln->seqs[i].title);
		ckfree(mult_aln->seqs[i].org);
		ckfree(mult_aln->seqs[i].data);
		ckfree(mult_aln->seqs[i].mask);
        	for(j=0;j<MAXFTTYPE;j++) {
        		for(k=0;k<mult_aln->ft[i].nentries[j];k++) {
				ckfree(mult_aln->ft[i].data[j][k].type);
				ckfree(mult_aln->ft[i].data[j][k].name);
			}
		}
        	for(k=0;k<mult_aln->go[i].ngorefs;k++) {
			ckfree(mult_aln->go[i].goref[k].id);
			ckfree(mult_aln->go[i].goref[k].desc);
		}
	}
	ckfree(mult_aln->seqs);
	ckfree(mult_aln->ft);
	ckfree(mult_aln->repeat);
	ckfree(mult_aln->go);
}

void alloc_seq(SEQ *seq,sint length)
{
	seq->data = (char *)ckalloc((length+2) * sizeof (char));
	seq->mask = (char *)ckalloc((length+2) * sizeof (char));
}

void realloc_seq(SEQ *seq,sint length)
{
	seq->data = (char *)ckrealloc(seq->data, (length+2) * sizeof (char));
	seq->mask = (char *)ckrealloc(seq->mask, (length+2) * sizeof (char));
}

void alloc_ft_entry(FT_ENTRY *data)
{
	data->type = (char *)ckalloc((10+2) * sizeof (char));
	data->name = (char *)ckalloc((100+2) * sizeof (char));
}

void alloc_go_entry(GOREF *goref)
{
	goref->id = (char *)ckalloc((10+2) * sizeof (char));
	goref->desc = (char *)ckalloc((100+2) * sizeof (char));
}

int getargs(char *inline1,char *args[],int max)
{

        char    *inptr;
/*
#ifndef MAC
        char    *strtok(char *s1, const char *s2);
#endif
*/
        int     i;

        inptr=inline1;
        for (i=0;i<=max;i++)
        {
                if ((args[i]=strtok(inptr," \t\n"))==NULL)
                        break;
                inptr=NULL;
        }

        return(i);
}

int getintargs(char *inline1,sint *args,int max)
{

        char    *inptr;
	char	*tstring;
/*
#ifndef MAC
        char    *strtok(char *s1, const char *s2);
#endif
*/
        int     i;

        inptr=inline1;
        for (i=0;i<=max;i++)
        {
                if ((tstring=strtok(inptr," \t\n"))==NULL)
                        break;
		args[i]=atoi(tstring);
                inptr=NULL;
        }

        return(i);
}


/* 
   count the number of identities between two sequences
	seq1	first sequence
	seq2	second sequence
*/
float countid(SEQ seq1,SEQ seq2)
{
   char c1,c2;
   sint i;
   sint count,total;
   float score;

   count = total = 0;
   for (i=0;i<seq1.len && i<seq2.len;i++) {
     c1 = seq1.data[i];
     c2 = seq2.data[i];
     if (isalpha(c1) && isalpha(c2)) {
       total++;
       if (c1 == c2) count++;
     }

   }

   if(total==0) score=0;
   else
   score = 100.0 * (float)count / (float)total;
   return(score);

}

/*
   count the number of identities between two sequences and divide by length of sequence
        seq1    first sequence
        seq2    second sequence
*/

float countid1(SEQ seq1,SEQ seq2)
{
   char c1,c2;
   sint i;
   sint count,total;
   sint len;
   float score;

   len = count = total = 0;
   for (i=0;i<seq1.len && i<seq2.len;i++) {
     c1 = seq1.data[i];
     c2 = seq2.data[i];
     if(isalpha(c1) || isalpha(c2)) len++;
     if (isalpha(c1) && isalpha(c2)) {
       total++;
       if (c1 == c2) count++;
     }

   }

   if(len==0) score=0;
   else
   score = (float)(100.0 * (float)count / (float)len);
   return(score);

}


void create_parameter_output(OPT opt,ALN mult_aln)
{
	char parname[FILENAMELEN+1], temp[FILENAMELEN+1];
	char path[FILENAMELEN+1];
	FILE *parout;
	Boolean usemenu;

        get_path(mult_aln.filename,path);
        strcpy(parname,path);
        strcat(parname,"par");

	usemenu=get_usemenu();
	if(usemenu) {
        	fprintf(stdout,"\nEnter a name for the parameter output file [%s]: ",
                                           parname);
               	gets(temp);
               	if(*temp != EOS)
                       	strcpy(parname,temp);
       	}

/* create a file with execute permissions first */
	remove(parname);
/*
	fd = creat(parname, 0777);
	close(fd);
*/

        if((parout = open_explicit_file(parname))==NULL) return;

        fprintf(parout,"clustalw \\\n");
	if ((mult_aln.nseqs>0) && (mult_aln.prf1.nseqs<=0)) fprintf(parout,"-infile=%s \\\n",mult_aln.filename);
	if (mult_aln.prf1.nseqs>0) fprintf(parout,"-profile1=%s\\\n",mult_aln.prf1.filename);
	if (mult_aln.prf2.nseqs>0) fprintf(parout,"-profile2=%s\\\n",mult_aln.prf2.filename);
	if (mult_aln.dnaflag == TRUE) fprintf(parout,"-type=dna \\\n");
	else                 fprintf(parout,"-type=protein \\\n");

	if (opt.pw_opt->quick_pairalign) {
		fprintf(parout,"-quicktree \\\n");
		if (!mult_aln.dnaflag) {
			fprintf(parout,"-ktuple=%d \\\n",(pint)opt.quickpw_opt->dna_ktup);
     			fprintf(parout,"-window=%d \\\n",(pint)opt.quickpw_opt->dna_window);
     			fprintf(parout,"-pairgap=%d \\\n",(pint)opt.quickpw_opt->dna_wind_gap);
     			fprintf(parout,"-topdiags=%d \\\n",(pint)opt.quickpw_opt->dna_signif);    
		}
		else {
			fprintf(parout,"-ktuple=%d \\\n",(pint)opt.quickpw_opt->prot_ktup);
     			fprintf(parout,"-window=%d \\\n",(pint)opt.quickpw_opt->prot_window);
     			fprintf(parout,"-pairgap=%d \\\n",(pint)opt.quickpw_opt->prot_wind_gap);
     			fprintf(parout,"-topdiags=%d \\\n",(pint)opt.quickpw_opt->prot_signif);    
		}
     		if (opt.quickpw_opt->percent) fprintf(parout,"/score=percent \\\n");      
     		else         fprintf(parout,"-score=absolute \\\n");      
	}
	else {
		if (!mult_aln.dnaflag) {
			fprintf(parout,"-pwmatrix=%s \\\n",opt.pw_opt->mtrxname);
			fprintf(parout,"-pwgapopen=%.2f \\\n",opt.pw_opt->prot_go_penalty);
			fprintf(parout,"-pwgapext=%.2f \\\n",opt.pw_opt->prot_ge_penalty);
		}
		else {
			fprintf(parout,"-pwgapopen=%.2f \\\n",opt.pw_opt->dna_go_penalty);
			fprintf(parout,"-pwgapext=%.2f \\\n",opt.pw_opt->dna_ge_penalty);
		}
	}

	if (!mult_aln.dnaflag) {
		fprintf(parout,"-matrix=%s \\\n",opt.mult_opt->mtrxname);
		fprintf(parout,"-gapopen=%.2f \\\n",opt.mult_opt->prot_gap_open);
		fprintf(parout,"-gapext=%.2f \\\n",opt.mult_opt->prot_gap_extend);
	}
	else {
		fprintf(parout,"-gapopen=%.2f \\\n",opt.mult_opt->dna_gap_open);
		fprintf(parout,"-gapext=%.2f \\\n",opt.mult_opt->dna_gap_extend);
	}

	fprintf(parout,"-maxdiv=%d \\\n",(pint)opt.mult_opt->divergence_cutoff);
	if (!opt.mult_opt->gap_opt->use_endgaps) fprintf(parout,"-endgaps \\\n");    

	if (!mult_aln.dnaflag) {
     		if (opt.mult_opt->neg_matrix) fprintf(parout,"-negative \\\n");   
     		if (opt.mult_opt->gap_opt->no_pref_penalties) fprintf(parout,"-nopgap \\\n");     
     		if (opt.mult_opt->gap_opt->no_hyd_penalties) fprintf(parout,"-nohgap \\\n");     
     		if (opt.mult_opt->gap_opt->no_var_penalties) fprintf(parout,"/novgap \\\n");     
    		fprintf(parout,"-hgapresidues=%s \\\n",opt.mult_opt->gap_opt->hyd_residues);
     		fprintf(parout,"-gapdist=%d \\\n",(pint)opt.mult_opt->gap_opt->gap_dist);     
	}
	else {
		fprintf(parout,"-transweight=%.2f \\\n",opt.mult_opt->transition_weight);
	}

     	if (opt.alnout_opt->output_gcg) fprintf(parout,"-output=gcg \\\n");
     	else if (opt.alnout_opt->output_gde) fprintf(parout,"-output=gde \\\n");
     	else if (opt.alnout_opt->output_nbrf) fprintf(parout,"-output=pir \\\n");
     	else if (opt.alnout_opt->output_phylip) fprintf(parout,"-output=phylip \\\n");
     	if (opt.alnout_opt->output_order==ALIGNED) fprintf(parout,"-outorder=aligned \\\n");  
     	else                      fprintf(parout,"-outorder=input \\\n");  
     	if (opt.alnout_opt->output_gde)
		if (opt.alnout_opt->lowercase) fprintf(parout,"-case=lower \\\n");
		else           fprintf(parout,"-case=upper \\\n");


        fprintf(parout,"-interactive\n");

	fclose(parout);

}




FILE *open_explicit_file(char *file_name)
{
        FILE * file_handle;

        if (*file_name == EOS) {
                error("Bad output file [%s]",file_name);
                return NULL;
        }
#ifdef VMS
        if((file_handle=fopen(file_name,"w","rat=cr","rfm=var"))==NULL) {
#else
        if((file_handle=fopen(file_name,"w"))==NULL) {
#endif
                error("Cannot open output file [%s]",file_name);
                return NULL;
        }
        return file_handle;
}

void check_fragments(ALNPTR mult_aln)
{
	sint i,j;
	sint *useqlen;
	float tmp,dscore;

/*
	useqlen = (sint *) ckalloc( (mult_aln->nseqs+1) * sizeof (sint) );
        for(i=0;i<mult_aln->nseqs;i++) {
                useqlen[i]=0;
                for(j=0;j<mult_aln->seqs[i].len;j++)
                        if(isalpha(mult_aln->seqs[i].data[j])) {
                                useqlen[i]++;
                        }
        }


        for (i=0;i<mult_aln->nseqs;i++) {
                for (j=i+1;j<mult_aln->nseqs;j++) {
                        dscore = countid(mult_aln->seqs[i],mult_aln->seqs[j]);
                        if(dscore>60) {
                                tmp=(float)useqlen[i]/(float)useqlen[j];
                                if(tmp<0.8) mult_aln->seqs[i].fragment=TRUE;
                                else if(tmp>1.25) mult_aln->seqs[j].fragment=TRUE;
                        }

                }
        }

        ckfree(useqlen);
*/

}

void remove_gap_pos(sint prf_no,ALNPTR mult_aln)
{
        int i,j,k,ngaps;
        int fseq,lseq;

        if(prf_no==0) {
                fseq=0;
                lseq=mult_aln->nseqs;
        }
        else if(prf_no==1) {
                fseq=0;
                lseq=mult_aln->prf1.nseqs;
        }
        else {
                fseq=mult_aln->prf1.nseqs;
                lseq=mult_aln->nseqs;
        }
        if(fseq>=lseq) return;

        for (i=0;i<mult_aln->seqs[fseq].len;)
        {
                ngaps=0;
                for (j=fseq;j<lseq;j++)
                        if(!isalpha(mult_aln->seqs[j].data[i])) ngaps++;
                if (ngaps==lseq-fseq)
                {
                        for (j=fseq;j<lseq;j++)
                        {
                                for(k=i+1;k<=mult_aln->seqs[j].len;k++)
                                        mult_aln->seqs[j].data[k-1]=mult_aln->seqs[j].data[k];
                                mult_aln->seqs[j].len--;
                        }
                        if(mult_aln->seqs[fseq].len<=0) break;
                }
                else i++;
        }
}

float normalise_score(float score,float n,float ntot,float ntotseq)
{
        float ret;

        if(n==0) ret=0.0;
        else
                ret=score*exp(-10.0*(float)(ntot-n)/((float)(ntot)));

        return ret;

}


void sort_scores(float *scores,int f,int l)
{
        int i,last;

        if(f>=l) return;

        swap_scores(scores,f,(f+l)/2);
        last=f;
        for(i=f+1;i<=l;i++)
        {
                if(scores[i]>scores[f])
                        swap_scores(scores,++last,i);
        }
        swap_scores(scores,f,last);
        sort_scores(scores,f,last-1);
        sort_scores(scores,last+1,l);

}

void swap_scores(float *scores,int s1, int s2)
{
        float temp;

        temp=scores[s1];
        scores[s1]=scores[s2];
        scores[s2]=temp;
}

void pos2col(char *seq,sint pstart,sint pend,sint *cstart,sint *cend)
{
        int i,ix;

        ix=0;
        if(pstart<0)
        {
                (*cstart)=-1;
                (*cend)=-1;
                return;
        }
        for(i=0;i<strlen(seq);i++)
        {
                if(isalpha(seq[i])) ix++;
                if(ix==pstart+1) break;
        }
        (*cstart)=i;

        if (pend<=pstart)
        {
                (*cend)=(*cstart);
                return;
        }

        i++;
        for(;i<strlen(seq);i++)
        {
                if(isalpha(seq[i])) ix++;
                if(ix==pend+1) break;
        }
        (*cend)=i;

}

void col2pos(char *seq,sint cstart,sint cend,sint *pstart,sint *pend)
{
        int i,ix;

        ix=0;
        if(cstart<0)
        {
                (*pstart)=-1;
                (*pend)=-1;
                return;
        }
        for(i=0;i<strlen(seq);i++)
        {
                if(isalpha(seq[i])) ix++;
                if(i==cstart) break;
        }
        (*pstart)=ix-1;

        if (cend<=cstart)
        {
                (*pend)=(*pstart);
                return;
        }

        i++;
        for(;i<strlen(seq);i++)
        {
                if(isalpha(seq[i])) ix++;
                if(i==cend) break;
        }
        (*pend)=ix-1;

}

sint overlap(sint f1,sint l1,sint f2, sint l2)
{
        sint len;

        len=MIN(l1,l2)-MAX(f1,f2)+1;
        return len;
}

sint check_ft_type(char *ft_type,char *ft_name,sint *type)
{
        sint ret=0;
        sint seq_type=0;

        if (strcmp(ft_type,"STRUCT")==0) {
                seq_type=3;
                (*type)=STRUCT;
                ret=1;
        }
        else if ((strcmp(ft_type,"HELIX")==0) ||
           (strcmp(ft_type,"STRAND")==0)) {
                strcpy(ft_name,ft_type);
                strcpy(ft_type,"STRUCT");
                seq_type=3;
                (*type)=STRUCT;
                ret=1;
        }
        else if ((strcmp(ft_type,"DNA_BIND")==0) ||
           (strcmp(ft_type,"ZN_FING")==0)) {
                strcpy(ft_name,ft_type);
                strcpy(ft_type,"DOMAIN");
                (*type)=SWDOMAIN;
                ret=1;
        }
        else if ((strcmp(ft_type,"NP_BIND")==0) ||
           (strcmp(ft_type,"CA_BIND")==0) ||
           (strcmp(ft_type,"METAL")==0) ||
           (strcmp(ft_type,"CARBOHYD")==0) ||
           (strcmp(ft_type,"DISULFID")==0) ||
           (strcmp(ft_type,"SITE")==0) ||
           (strcmp(ft_type,"ACT_SITE")==0) ||
           (strcmp(ft_type,"BINDING")==0)) { 
                strcpy(ft_name,ft_type);
                strcpy(ft_type,"SITE");
                (*type)=SITE;
                ret=1;
        }
        else if (strcmp(ft_type,"MOD_RES")==0) {
               	(*type)=MODRES;
               	ret=1;
        }
        else if (strcmp(ft_type,"REPEAT")==0) {
                (*type)=REPEAT;
                ret=1;
        }
        else if (strcmp(ft_type,"TRANSMEM")==0) {
                (*type)=TRANSMEM;
                ret=1;
        }
        else if (strcmp(ft_type,"SIGNAL")==0) {
                (*type)=SIGNAL;
                ret=1;
        }
        else if (strcmp(ft_type,"VARSPLIC")==0) {
                (*type)=VARSPLIC;
                ret=1;
        }
        else if (strcmp(ft_type,"DOMAIN")==0) {
                (*type)=SWDOMAIN;
                ret=1;
        }
        else if (strcmp(ft_type,"BLOCK")==0) {
                (*type)=COREBLOCK;
                ret=1;
        }
        else if (strcmp(ft_type,"REGION")==0) {
                (*type)=REGION;
                ret=1;
        }
        else if (strcmp(ft_type,"SEQERR")==0) {
                (*type)=SEQERRBLOCK;
                ret=1;
        }
        else if (strcmp(ft_type,"COIL")==0) {
                (*type)=COIL;
                ret=1;
        }
        else if (strcmp(ft_type,"LOWCOMP")==0) {
                (*type)=LOWC;
                ret=1;
        }
        else if (strcmp(ft_type,"PFAM-A")==0) {
                (*type)=PFAMA;
                ret=1;
        }
        else if (strcmp(ft_type,"PFAM-B")==0) {
                (*type)=PFAMB;
                ret=1;
        }
        else if (strcmp(ft_type,"PROSITE")==0) {
                (*type)=PROSITE;
                ret=1;
        }
        return ret;
}

void add_ft_entry(ALNPTR mult_aln,sint seq,sint first,sint last,sint type,sint code,float score,char *ctype,char *iname,sint is,sint ie)
{
        sint n,color;
        sint fr,lr;
        sint fc,lc;
	char name[100];

        if(mult_aln->ft[seq].nentries[type]>MAXFT) {
                fprintf(stdout,"WARNING: too many features in %s %d (%d)\n",mult_aln->seqs[seq].name,mult_aln->ft[seq].nentries[type],type);
                return;
        }
        if(last<is || first>ie) return;
        if(first<is) first=is;
        if(last>ie) last=ie;
        n=mult_aln->ft[seq].nentries[type];
        alloc_ft_entry(&mult_aln->ft[seq].data[type][n]);
        strcpy(mult_aln->ft[seq].data[type][n].type,ctype);
        col2pos(mult_aln->seqs[seq].data,first,last,&fr,&lr);
        pos2col(mult_aln->seqs[seq].data,fr,lr,&fc,&lc);
        if(fc==first)
                mult_aln->ft[seq].data[type][n].start=fr;
        else
                mult_aln->ft[seq].data[type][n].start=fr+1;
        if(lc==last)
                mult_aln->ft[seq].data[type][n].end=lr;
        else
                mult_aln->ft[seq].data[type][n].end=lr-1;

	strcpy(name,iname);
	if(strcmp(ctype,"BLOCK")==0) {
		if (strcmp(name,"SBLOCK")==0) 
			color=code;
		else if (strcmp(name,"LBLOCK")==0) 
			color=12;
	}
	else if(strcmp(ctype,"REGION")==0) {
		color=10;
	}
	else if(strcmp(ctype,"SEQERR")==0) {
		color=3;
	}
        else if(strcmp(ctype,"ANCHOR")==0) {
                color=code;
        }
        else if(strcmp(ctype,"REPEAT")==0) {
                color=code;
        }
	else if(strcmp(ctype,"STRUCT")==0) {
		if(strncmp(name,"PRED_HELIX",10)==0) color=VLRED;
		else if(strncmp(name,"PRED_STRAND",11)==0) color=VLGREEN;
		else if(strncmp(name,"PROP_HELIX",10)==0) color=LRED;
		else if(strncmp(name,"PROP_STRAND",11)==11) color=LGREEN;
		else if(strcmp(name,"HELIX")==0) color=RED;
		else if(strcmp(name,"STRAND")==0) color=GREEN;
	}
	else if(strcmp(ctype,"TRANSMEM")==0) {
		if(strncmp(name,"PRED_",5)==0) {
			sprintf(name,"%s %.2f",name,score);
			color=10;
		}
		else if(strncmp(name,"PROP_",5)==0) color=GRAY;
		else color=1;
	}
	else if(strcmp(ctype,"COIL")==0) {
		if(strncmp(name,"PRED_",5)==0) color=LGRAY;
		else if(strncmp(name,"PROP_",5)==0) color=GRAY;
		else color=6;
	}
        else if(strcmp(ctype,"LOWCOMP")==0) {
                if(strncmp(name,"PRED_",5)==0) color=LGRAY;
                else if(strncmp(name,"PROP_",5)==0) color=GRAY;
                else color=6;
        }
        else if(strcmp(ctype,"DOMAIN")==0) {
                if(strncmp(name,"PRED_",5)==0) color=LGRAY;
                else if(strncmp(name,"PROP_",5)==0) color=GRAY;
                else if(strncmp(name,"DNA_BIND",8)==0) color=2;
                else if(strncmp(name,"ZN_FING",7)==0) color=3;
                else color=1;
        }
        else if(strcmp(ctype,"PFAM-A")==0) {
                if(strncmp(name,"PRED_",5)==0) color=LGRAY;
                else if(strncmp(name,"PROP_",5)==0) color=GRAY;
                else color=code;
        }
        else if(strcmp(ctype,"PFAM-B")==0) {
                if(strncmp(name,"PRED_",5)==0) color=LGRAY;
                else if(strncmp(name,"PROP_",5)==0) color=GRAY;
                else color=code;
        }
        else if(strcmp(ctype,"PROSITE")==0) {
                if(strncmp(name,"PRED_",5)==0) color=LGRAY;
                else if(strncmp(name,"PROP_",5)==0) color=GRAY;
                else color=code;
        }
        else if(strcmp(ctype,"SIGNAL")==0) {
                if(strncmp(name,"PRED_",5)==0) {
                        sprintf(name,"%s %.2f",name,score);
                        color=LGRAY;
                }
                else if(strncmp(name,"PROP_",5)==0) color=GRAY;
                else if(strncmp(name,"NUCLEAR",7)==0) color=1;
                else if(strncmp(name,"CHLOROPLAST",11)==0) color=2;
                else if(strncmp(name,"MITOCHONDRION",13)==0) color=3;
                else if(strncmp(name,"MICROBODY",9)==0) color=4;
                else color=5;
        }
        else if(strcmp(ctype,"VARSPLIC")==0) {
                if(strncmp(name,"PRED_",5)==0) color=LGRAY;
                else if(strncmp(name,"PROP_",5)==0) color=GRAY;
                else color=5;
        }
        else if(strcmp(ctype,"MOD_RES")==0) {
                color=0;
                if(strncmp(name,"PRED_",5)==0) color=LGRAY;
                else if(strncmp(name,"PROP_",5)==0) color=GRAY;
                else if(strncmp(name,"ACETYLATION",11)==0) color=1;
                else if(strncmp(name,"ALKYLATION",10)==0) color=2;
                else if(strncmp(name,"HYDROXYLATION",13)==0) color=3;
                else if(strncmp(name,"METHYLATION",11)==0) color=4;
                else if(strncmp(name,"PHOSPHORYLATION",15)==0) color=5;
                else if(strncmp(name,"SULFATION",9)==0) color=6;
                else if(strncmp(name,"PYRROLIDONE",11)==0) color=7;
                else if(strncmp(name,"MYRISTATE",9)==0) color=8;
                else if(strncmp(name,"PALMITATE",9)==0) color=9;
                else if(strncmp(name,"FARNESYL",8)==0) color=10;
                else if(strncmp(name,"GERANYL",7)==0) color=11;
                else if(strncmp(name,"GPI-ANCHOR",3)==0) color=12;
                else if(strncmp(name,"N-ACYL",5)==0) color=13;
        }
        else if(strcmp(ctype,"SITE")==0) {
                color=0;
                if(strncmp(name,"PRED_",5)==0) color=LGRAY;
                else if(strncmp(name,"PROP_",5)==0) color=GRAY;
                else if(strncmp(name,"SITE",4)==0) color=0;
                else if(strncmp(name,"ACT_SITE",8)==0) color=1;
                else if(strncmp(name,"CARBOHYD",8)==0) color=2;
                else if(strncmp(name,"DISULFID",8)==0) color=3;
                else if(strncmp(name,"BINDING",7)==0) color=4;
                else if(strncmp(name,"NP_BIND",7)==0) color=5;
                else if(strncmp(name,"CA_BIND",7)==0) color=6;
                else if(strncmp(name,"TRANSIT",7)==0) color=7;
                else if(strncmp(name,"LIPID",5)==0) color=8;
                else if(strncmp(name,"METAL",5)==0) color=9;
        }

        strcpy(mult_aln->ft[seq].data[type][n].name,name);
        mult_aln->ft[seq].data[type][n].score=score;
        mult_aln->ft[seq].data[type][n].color=color;






        mult_aln->ft[seq].nentries[type]++;

}

