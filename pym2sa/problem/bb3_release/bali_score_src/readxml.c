#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "clustalw.h"
#include <expat.h>

#define BUFFSIZE	8192 	/* size of buffer for reading xml file */

char Buff[BUFFSIZE]; 	/* buffer for reading xml file */

int Depth;		/* depth of xml tree (not currently used!) */
char *Element;		/* element name */
char *Content;		/* element contents */
char **Attributes;	/* list of attributes */
int  NAttributes;	/* number of attributes */
int Nseqs;		/* number of sequences */
int Nftable[MAXFTTYPE];	/* number of entries in feature table */
int Ftype;		/* feature type */
int Ngorefs;		/* number of GO cross-references */
int Nkeywords;		/* number of keywords */
int Ncolscores;		/* number of column scores */
int Fseq;		/* first sequence - start index to Maln */
ALN Maln;

#define START 0
#define END 1

static void do_element(void);

/* This is the application specific stuff. When we get here, the element name is contained in
Element, any attributes are in NAttributes and Attributes, the data is in Content */

static void do_element(void)
{
	int i,j,l,n;

	if(strcmp(Element,"sequence")==0) {
/* this is the start of a new sequence, so reset the data for this sequence */
		for(i=0;i<MAXFTTYPE;i++)
			Nftable[i]=(-1);
		Ngorefs=(-1);
		Nkeywords=(-1);
		Ncolscores=(-1);
	}
	else if(strcmp(Element,"seq-name")==0) {
		strcpy(Maln.seqs[Fseq+Nseqs].name,Content);
	}
	else if(strcmp(Element,"accession")==0) {
		strcpy(Maln.seqs[Fseq+Nseqs].access,Content);
	}
	else if(strcmp(Element,"nid")==0) {
		strcpy(Maln.seqs[Fseq+Nseqs].nid,Content);
	}
	else if(strcmp(Element,"definition")==0) {
		strcpy(Maln.seqs[Fseq+Nseqs].title,Content);
	}
	else if(strcmp(Element,"organism")==0) {
		strcpy(Maln.seqs[Fseq+Nseqs].org,Content);
	}
	else if(strcmp(Element,"subgroup")==0) {
		Maln.seqs[Fseq+Nseqs].simgroup=atoi(Content);
	}
	else if(strcmp(Element,"hydrophobicity")==0) {
		Maln.seqs[Fseq+Nseqs].hydrophobicity=atof(Content);
	}
	else if(strcmp(Element,"fitem")==0) {
	}
	else if(strcmp(Element,"ftype")==0) {
		if(strcmp(Content,"DOMAIN")==0) 
			Ftype=SWDOMAIN;
		else if (strcmp(Content,"PFAM-A")==0) 
			Ftype=PFAMA;
		else if (strcmp(Content,"PFAM-B")==0) 
			Ftype=PFAMB;
		else if (strcmp(Content,"PROSITE")==0) 
			Ftype=PROSITE;
		else if (strcmp(Content,"SITE")==0) 
			Ftype=SITE;
		else if (strcmp(Content,"TRANSMEM")==0) 
			Ftype=TRANSMEM;
		else if (strcmp(Content,"REPEAT")==0) 
			Ftype=REPEAT;
		else if (strcmp(Content,"COIL")==0) 
			Ftype=COIL;
		else if (strcmp(Content,"SIGNAL")==0) 
			Ftype=SIGNAL;
		else if (strcmp(Content,"STRUCT")==0) 
			Ftype=STRUCT;
		else if (strcmp(Content,"VARSPLIC")==0) 
			Ftype=VARSPLIC;
		else if (strcmp(Content,"ANCHOR")==0) 
			Ftype=ANCHOR;
		else if (strcmp(Content,"BLOCK")==0) 
			Ftype=COREBLOCK;
		else if (strcmp(Content,"SEQERR")==0) 
			Ftype=SEQERRBLOCK;
		else if (strcmp(Content,"REGION")==0) 
			Ftype=REGION;
		else if (strcmp(Content,"MOD_RES")==0) 
			Ftype=MODRES;
		else if (strcmp(Content,"LOWCOMP")==0) 
			Ftype=LOWC;
		else if (strcmp(Content,"UNKNOWN")==0) 
			Ftype=OTHER;
		Nftable[Ftype]++;
		if(Nftable[Ftype]>=MAXFT) 
			warning("too many FTABLE elements for sequence %s",Maln.seqs[Fseq+Nseqs].name);
		else alloc_ft_entry(&Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]]);
		if(Nftable[Ftype]<MAXFT) strcpy(Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]].type,Content);
	}
	else if(strcmp(Element,"fstart")==0) {
		if(Nftable[Ftype]<MAXFT) Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]].start=atoi(Content)-1;
	}
	else if(strcmp(Element,"fstop")==0) {
		if(Nftable[Ftype]<MAXFT) Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]].end=atoi(Content)-1;
	}
	else if(strcmp(Element,"fcolor")==0) {
		if(Nftable[Ftype]<MAXFT) Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]].color=atoi(Content);
	}
	else if(strcmp(Element,"fscore")==0) {
		if(Nftable[Ftype]<MAXFT) Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]].score=atof(Content);
	}
	else if(strcmp(Element,"fnote")==0) {
		if(Nftable[Ftype]<MAXFT) if(Content!=NULL) strcpy(Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]].name,Content);
	}
	else if(strcmp(Element,"column-score")==0) {
		if(Ncolscores<MAXCSCORE) Ncolscores++;
	}
	else if(strcmp(Element,"colsco-name")==0) {
		if(Ncolscores<MAXCSCORE) if(Content!=NULL) {
			n=strlen(Content);
			Maln.col_score[Ncolscores].name=(char *)ckalloc((n+1)*sizeof(char));
			strcpy(Maln.col_score[Ncolscores].name,Content);
		}
	}
	else if(strcmp(Element,"colsco-owner")==0) {
		if(Ncolscores<MAXCSCORE) if(Content!=NULL) {
			n=strlen(Content);
			Maln.col_score[Ncolscores].owner=(char *)ckalloc((n+1)*sizeof(char));
			strcpy(Maln.col_score[Ncolscores].owner,Content);
		}
	}
	else if(strcmp(Element,"colsco-data")==0) {
		if(Ncolscores<MAXCSCORE) {
			n=0;
			if(Content!=NULL) {
				for(i=0;i<strlen(Content);i++)
					if(isspace(Content[i])) n++;
				Maln.col_score[Ncolscores].data=(sint *)ckalloc((n+12)*sizeof(sint));
				i=getintargs(Content,Maln.col_score[Ncolscores].data,n+10);
				Maln.col_score[Ncolscores].length=n;
				Maln.ncol_scores=Ncolscores+1;
			}
		}
	}
	else if(strcmp(Element,"goxref")==0) {
		Ngorefs++;
		if(Ngorefs>=MAXGOREF) 
			warning("too many GOXREF elements for sequence %s",Maln.seqs[Fseq+Nseqs].name);
	
	}
	else if(strcmp(Element,"goid")==0) {
		if(Ngorefs<MAXGOREF) {
		l=strlen(Content);
		Maln.go[Fseq+Nseqs].goref[Ngorefs].id=(void *)ckalloc((l+1)*sizeof(char));
		strcpy(Maln.go[Fseq+Nseqs].goref[Ngorefs].id,Content);
		}
	}
	else if(strcmp(Element,"goclass")==0) {
		if(Ngorefs<MAXGOREF) {
		Maln.go[Fseq+Nseqs].goref[Ngorefs].class=Content[0];
		}
	}
	else if(strcmp(Element,"godesc")==0) {
		if(Ngorefs<MAXGOREF) {
		l=strlen(Content);
		Maln.go[Fseq+Nseqs].goref[Ngorefs].desc=(void *)ckalloc((l+1)*sizeof(char));
		strcpy(Maln.go[Fseq+Nseqs].goref[Ngorefs].desc,Content);
		}
	}
	else if(strcmp(Element,"goevidence")==0) {
		if(Ngorefs<MAXGOREF) {
		strncpy(Maln.go[Fseq+Nseqs].goref[Ngorefs].evidence,Content,3);
		}
	}
	else if(strcmp(Element,"length")==0) {
		Maln.seqs[Fseq+Nseqs].len=atoi(Content);
	}
	else if(strcmp(Element,"group")==0) {
		Maln.seqs[Fseq+Nseqs].simgroup=atoi(Content);
	}
	else if(strcmp(Element,"keyword")==0) {
		Nkeywords++;
		if(Nkeywords>=MAXKEYWORDS) {
			warning("too many KEYWORD elements for sequence %s",Maln.seqs[Fseq+Nseqs].name);
		}
		else {
		l=strlen(Content);
		Maln.seqs[Fseq+Nseqs].keyword[Nkeywords]=(void *)ckalloc((l+1)*sizeof(char));
		strcpy(Maln.seqs[Fseq+Nseqs].keyword[Nkeywords],Content);
		}
	}
	else if(strcmp(Element,"seq-data")==0) {
/* then this sequence must be finished, so we can update the data for this sequence */
		Maln.seqs[Fseq+Nseqs].nkeywords=Nkeywords+1;
		Maln.go[Fseq+Nseqs].ngorefs=Ngorefs+1;
		for(i=0;i<MAXFTTYPE;i++)
			Maln.ft[Fseq+Nseqs].nentries[i]=Nftable[i]+1;
		l=strlen(Content);
		alloc_seq(&Maln.seqs[Fseq+Nseqs],l);
		for(i=0,j=0;i<l;i++)
			if(!isspace(Content[i])) Maln.seqs[Fseq+Nseqs].data[j++]=Content[i];
		Maln.seqs[Fseq+Nseqs].len=j;
		Maln.nseqs++;
		Nseqs++;
	}

	for(i=0;i<NAttributes;i++) {
		ckfree(Attributes[i*2]);
		ckfree(Attributes[i*2+1]);
	}
	if(Element!=NULL) Element=ckfree(Element);
	if(Content!=NULL) Content=ckfree(Content);
}

/* This only gets called at the end of an element that contains data. Shouldn't need to touch it again! */
static void end(void *data, const char *el)
{
	if(Element!=NULL) {
		do_element();
	}
	Depth--;
}  

/* This reads the element names and puts them in Element. Shouldn't need to touch it again! */
static void start(void *data, const char *el, const char **attr)
{
	int i,n,len;

	if(Element!=NULL) {
		do_element();
	}

	len=strlen(el);
	Element=(char *)ckalloc((len+1)*sizeof(char));
	strcpy(Element,el);

	n=0;
	for (i = 0; attr[i]; i += 2) {
		n++;
	}
	NAttributes=n;

	if(NAttributes>0) {
		Attributes=(char **)ckalloc((NAttributes*2)*sizeof(char *));
		for (i = 0; attr[i]; i += 2) {
			len=strlen(attr[i]);
			Attributes[i]=(char *)ckalloc((len+1)*sizeof(char));
			strcpy(Attributes[i],attr[i]);
			len=strlen(attr[i+1]);
			Attributes[i+1]=(char *)ckalloc((len+1)*sizeof(char));
			strcpy(Attributes[i+1],attr[i+1]);
		}
	}
/*
	for (i = 0; attr[i]; i += 2) {
    printf(" %s='%s'", attr[i], attr[i + 1]);
	}
*/

	Depth++;
}  

/* This reads the contents of each element, and puts into Content. Shouldn't need to touch it again! */
static void charhndl(void *data, const char *text, int len)
{
	int i,l=0;

	if(len<=0) return;

	if(Content==NULL)
		Content=(char *)ckalloc((len+1)*sizeof(char));
	else {
		l=strlen(Content);
		Content=(char *)ckrealloc(Content,(l+len+2)*sizeof(char));
	}
	for (i = 0; i < len; i++)
		Content[l++] = text[i];
	if(Content[l-1]=='\n') Content[l-1]='\0';
	else Content[l]='\0';
}  

/* handlers used, if we're only counting the number of sequences */
static void c_start(void *data, const char *el, const char **attr)
{
	int i,len;

	if(strcmp(el,"sequence")==0) Nseqs++;

	Depth++;
}  

static void c_end(void *data, const char *el)
{
	Depth--;
}  

sint count_xml_seqs(FILE *fin)
{
  XML_Parser p;

  p = XML_ParserCreate(NULL);

  if (! p) {
    fprintf(stderr, "Couldn't allocate memory for parser\n");
    exit(-1);
  }

  XML_SetElementHandler(p, c_start, c_end);

  Nseqs=0;
  for (;;) {
    int done;
    int len;

    len = fread(Buff, 1, BUFFSIZE, fin);
    if (ferror(fin)) {
      fprintf(stderr, "Read error\n");
      exit(-1);
    }
    done = feof(fin);

    if (! XML_Parse(p, Buff, len, done)) {
      fprintf(stderr, "Parse error at line %d:\n%s\n",
	      XML_GetCurrentLineNumber(p),
	      XML_ErrorString(XML_GetErrorCode(p)));
      exit(-1);
    }

    if (done)
      break;
  }
  return Nseqs;
}

ALN read_xml(FILE *fin,int first_seq)
{
  XML_Parser p;

  p = XML_ParserCreate(NULL);

  if (! p) {
    fprintf(stderr, "Couldn't allocate memory for parser\n");
    exit(-1);
  }
  alloc_aln(Nseqs,&Maln);

  XML_SetElementHandler(p, start, end);
  XML_SetCharacterDataHandler(p, charhndl);

  Nseqs=0;
  Fseq=first_seq;
  for (;;) {
    int done;
    int len;

    len = fread(Buff, 1, BUFFSIZE, fin);
    if (ferror(fin)) {
      fprintf(stderr, "Read error\n");
      exit(-1);
    }
    done = feof(fin);

    if (! XML_Parse(p, Buff, len, done)) {
      fprintf(stderr, "Parse error at line %d:\n%s\n",
	      XML_GetCurrentLineNumber(p),
	      XML_ErrorString(XML_GetErrorCode(p)));
      exit(-1);
    }

    if (done)
      break;
  }
  return Maln;
}
