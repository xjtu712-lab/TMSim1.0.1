
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <iostream> 
#include <fstream>

using namespace std;

int 	REF_LEN;

FILE* ingenfp;
FILE* insimfp;
FILE* outbedfp;
FILE* outbedannofp;


#define INIT_SEQ(seq) (seq).s = 0; (seq).l = (seq).max = 0
#define CLEAR_SEQ(seq) free( (seq).s )
#define INIT_MUTSEQ(mseq) (mseq).s = 0; (mseq).l = (mseq).max = 0; (mseq).gnum = (mseq).vnum = 0
#define CLEAR_MUTSEQ(mseq) free( (mseq).s )

enum muttype_t {NOCHANGE = 0, INSERT = 0x80, LONGINSERT = 0xf0000000, SUBSTITUTE = 0x40, DELETE = 0xc0, GTYPEBB=0x20, GTYPEAB=0x10};
typedef unsigned int mut_t;

static mut_t mutmsk = (mut_t)0xc0;
static int SEQ_BLOCK_SIZE = 512;

typedef struct {
	int l, max;
	char *s;
} seq_t;

typedef struct {
	int indel_no;
	int indel_beg, indel_len, indel_type, indel_num, indel_var;
	char *str;
} indel_e;

typedef struct {
	int l, max;
	indel_e *s;
} indel_t;

 int *label;
 mut_t *hap1;
 mut_t *hap2;
int sim_len;
indel_t  indeltab;
/**********************************************/
FILE *FileOpen(char *fileName, const char *mode)
{
	FILE *fp;
	fp = fopen (fileName, mode);
	if (fp == NULL)
	{
		fprintf(stdout, "Error: Cannot Open the file %s for %s !\n", fileName, mode);
		fflush(stdout);
		exit(0);
	}
	return fp;
}

/**********************************************/
char gch_rev(char orgch){
	char ch;
	orgch = _toupper(orgch);
	switch(orgch){
		case 'A': ch='T'; break;
		case 'C': ch='G'; break;
		case 'G': ch='C'; break;
		case 'T': ch='A'; break;
		default : ch='N'; break;
	}
	return ch;
}

/**********************************************/
int gch_to_num(char orgch){
	int ch;
	orgch = _toupper(orgch);
	switch(orgch){
		case 'A': ch=0; break;
		case 'C': ch=1; break;
		case 'G': ch=2; break;
		case 'T': ch=3; break;
		default : ch=4; break;
	}
	return ch;
}
/**********************************************/
char num_to_gch(int orgch){
	char ch;
	switch(orgch){
		case 0:   ch='A'; break;
		case 1:   ch='C'; break;
		case 2:   ch='G'; break;
		case 3:   ch='T'; break;
		default : ch='N'; break;
	}
	return ch;
}

///**************************************************************
int seq_read_fasta(FILE *fp, seq_t *seq, char *chrname, char *comment)
{
	int c, l, max, n;
	char *p, cht[2];
	
	c = 0;
	while (!feof(fp) && fgetc(fp) != '>');
	if (feof(fp)) return -1;
	p = chrname;
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
	n=0; cht[1] = '\0';
	l = 0; max = seq->max;
	while (!feof(fp) && (c = fgetc(fp)) ) {
		if (isalpha(c) || c == '-' || c == '.') {
			if (l + 1 >= max) {
				max += SEQ_BLOCK_SIZE;
				seq->s = (char *)realloc(seq->s, sizeof(char) * max);
			}
			cht[0] = seq->s[l++] = _toupper((char)c);
			if (strstr("ACGT", cht)) n++;
		}
	}
	seq->s[l] = 0;
	seq->max = max; 
	seq->l = l;
	
	return l;
}


/*******************************************************/
int addinsert_iterm(indel_t *indelp, int idno, int idbeg, int idlen, int gtype, int idnum, int idvar, char *idstr){
int pos;
	pos = indelp->l;
	if (pos+1 >= indelp->max){
		indelp->max += SEQ_BLOCK_SIZE;
		indelp->s = (indel_e *)realloc(indelp->s, sizeof(indel_e) * indelp->max);
	}
	indelp->s[pos].indel_no =idno;
	indelp->s[pos].indel_beg=idbeg;
	indelp->s[pos].indel_len=idlen;
	indelp->s[pos].indel_type=gtype;
	indelp->s[pos].indel_num=idnum;
	indelp->s[pos].indel_var=idvar;
	indelp->s[pos].str = (char *)calloc(strlen(idstr)+1, sizeof(char));
	strcpy(indelp->s[pos].str, idstr);
	indelp->l++;
	return 0;
}

/********************************************************************/
int read_longinsert_tab(char *insimf){
int idno, idbeg, idlen, idtype, idnum, idvar;
int n;
char *idstr, intabfile[50];
FILE *fp;

	sprintf(intabfile,"%s.idx",  insimf);
	fp = fopen( intabfile, "rt");
	if (!fp) return 0;
	idstr = new char [100000];
	n=0; 
	while ( !feof(fp) && fscanf(fp, "%d %d %d %d %d %d %s\n", &idno, &idbeg, &idlen, &idtype, &idnum, &idvar, idstr) != EOF ){
		addinsert_iterm(&indeltab, idno, idbeg, idlen, idtype, idnum, idvar, idstr);
		n++;
	} //while (reading file .....)
	printf("END read long insert table::  Long insert num=%d     Long insert table size=%d\n", n, indeltab.l);
	delete [] idstr;
	fclose(fp);
	return 0;
}

/**********************************************/

void print_bedfile( FILE *fp1, FILE *fp2 ,seq_t *seq,  char *refname){
 
    int i,j,pos,tmp;
	int mark;
	char *lstr; 
	lstr = new char [REF_LEN];
	char gt;//genotype
	char h1_h2[10];
	char type[5];//mutation  type
	char var;//variation base
	int is;//insert size


	//format: #CHROM	chromStrart		chromEnd
	//format: chrom  chromStart  chromEnd  type  ref  var  h1/h2  genotype
   
	for(i=0;i<sim_len;i++){
        pos=label[i];mark=0;
		if((hap1[i]&mutmsk)==SUBSTITUTE) mark=0;
		if((hap2[i]&mutmsk)==SUBSTITUTE) mark=1;
		if((hap1[i]&mutmsk)==INSERT && (hap1[i]&LONGINSERT)!=LONGINSERT ) mark=2;
		if((hap2[i]&mutmsk)==INSERT && (hap2[i]&LONGINSERT)!=LONGINSERT ) mark=3;
        if((hap1[i]&mutmsk)==DELETE) mark=4;
		if((hap2[i]&mutmsk)==DELETE) mark=5;
 		if( (hap1[i]&LONGINSERT)==LONGINSERT ) mark=6;
		if( (hap2[i]&LONGINSERT)==LONGINSERT ) mark=7;

		if( hap1[i]!= hap2[i] ){//AB
			gt='1';	
			switch(mark){	
			case 0://hap1 or hap2 subsitute
					strcpy(type,"SNP");
					strcpy(h1_h2,"1/0");
                    var=num_to_gch(hap1[i]&0x3);
					fprintf(fp1,"%s\t%d\t%d\n",refname,pos,pos);
					fprintf(fp2,"%s\t%d\t%d\t%s\t%c\t%c\t%s\t%c\n",refname,pos,pos,type,seq->s[pos],var,h1_h2,gt);
					break;
			case 1:
 					strcpy(type,"SNP");
					strcpy(h1_h2,"0/1");
                    var=num_to_gch(hap2[i]&0x3);
					fprintf(fp1,"%s\t%d\t%d\n",refname,pos,pos);
					fprintf(fp2,"%s\t%d\t%d\t%s\t%c\t%c\t%s\t%c\n",refname,pos,pos,type,seq->s[pos],var,h1_h2,gt);
					break;
										
			case 2://hap1 or hap2 insertion
				    lstr[0] = num_to_gch( hap2[i]&0x3 );
					strcpy(type,"INS");
					strcpy(h1_h2,"1/0");
					is = (hap1[i]>>28);
					tmp = ( (hap1[i]<<4)>>12 );
					for(j=0;j<is;j++){ lstr[is-j] = num_to_gch( tmp&0x3 ); tmp=tmp>>2; }
					lstr[is]='\0';
					fprintf(fp1,"%s\t%d\t%d\n",refname,pos,pos+1);
					fprintf(fp2,"%s\t%d\t%d\t%s\t%c\t%s\t%s\t%c\n",refname,pos,pos+1,type,lstr[0],lstr,h1_h2,gt);
					break;
			case 3:
					lstr[0] = num_to_gch( hap2[i]&0x3 );
					strcpy(type,"INS");
					strcpy(h1_h2,"0/1");
					is = (hap2[i]>>28);
					tmp = ( (hap2[i]<<4)>>12 );
					for(j=0;j<is;j++){ lstr[is-j] = num_to_gch( tmp&0x3 ); tmp=tmp>>2; }
					lstr[is]='\0';
					fprintf(fp1,"%s\t%d\t%d\n",refname,pos,pos+1);
					fprintf(fp2,"%s\t%d\t%d\t%s\t%c\t%s\t%s\t%c\n",refname,pos,pos+1,type,lstr[0],lstr,h1_h2,gt);
					break;
						
			case 4://hap1 or hap2 deletion
					strcpy(type,"DEL");
					strcpy(h1_h2,"1/0");
					lstr[0] =  seq->s[pos-1];
					lstr[1] = num_to_gch( hap1[i]&0x3 );
					tmp=2;
					for(j=i+1; j<sim_len; j++){
						if((hap1[j]&mutmsk)==DELETE && (hap2[j]&mutmsk)!=DELETE && label[j]==label[j-1]+1 ) lstr[tmp++]=num_to_gch(hap1[j]&0x3);
						else break;
					}
					lstr[tmp]='\0';
					i = j-1;
					fprintf(fp1,"%s\t%d\t%d\n",refname,pos,pos+tmp);
					fprintf(fp2,"%s\t%d\t%d\t%s\t%s\t%c\t%s\t%c\n",refname,pos,pos+tmp,type,lstr,lstr[0] ,h1_h2,gt);
					break;
			case 5:
					strcpy(type,"DEL");
					strcpy(h1_h2,"0/1");
					lstr[0] =  seq->s[pos-1];
					lstr[1] = num_to_gch( hap2[i]&0x3 );
					tmp=2;
					for(j=i+1; j<sim_len; j++){
						if((hap2[j]&mutmsk)==DELETE && (hap1[j]&mutmsk)!=DELETE && label[j]==label[j-1]+1 ) lstr[tmp++]=num_to_gch(hap2[j]&0x3);
						else break;
					}
					lstr[tmp]='\0';
					i = j-1;
					fprintf(fp1,"%s\t%d\t%d\n",refname,pos,pos+tmp);
					fprintf(fp2,"%s\t%d\t%d\t%s\t%s\t%c\t%s\t%c\n",refname,pos,pos+tmp,type,lstr,lstr[0] ,h1_h2,gt);
					break;				
			case 6://hap1 long insertion
				{
					strcpy(type,"LONGINS");
					strcpy(h1_h2,"1/0");
					int in_no, no, j, k, l,cc, loopmk;

					lstr[0] = num_to_gch( hap1[i]&0x3 );
					l = 1;
					no = ((hap1[i]<<4)>>12);
					loopmk = 0;
					for(j=0; j<indeltab.l; j++){
						if (indeltab.s[j].indel_no == no){
							in_no = indeltab.s[j].indel_len;
							//printf("no=%d ==>%d %d %d %d %d %d %d \n", no, indeltab.s[j].indel_no, indeltab.s[j].indel_beg, indeltab.s[j].indel_end, 
					//indeltab.s[j].indel_len, indeltab.s[j].indel_type, indeltab.s[j].indel_num, indeltab.s[j].indel_var);
							for (k=0; k<in_no; k++){
								lstr[l++] = indeltab.s[j].str[k];
							}
							loopmk = 1;
						}else if (loopmk) break;
					}
					lstr[l] ='\0';
					fprintf(fp1,"%s\t%d\t%d\n",refname,pos,pos+1);
					fprintf(fp2,"%s\t%d\t%d\t%s\t%c\t%s\t%s\t%c\n",refname,pos,pos+1,type,lstr[0],lstr,h1_h2,gt);
					
				}
					break;	
			case 7://hap2 long insertion
				{
					strcpy(type,"LONGINS");
					strcpy(h1_h2,"0/1");				
					int in_no, no, j, k, l,cc, loopmk;

					lstr[0] = num_to_gch( hap2[i]&0x3 );
					l = 1;
					no = ((hap2[i]<<4)>>12);
					loopmk = 0;
					for(j=0; j<indeltab.l; j++){
						if (indeltab.s[j].indel_no == no){
							in_no = indeltab.s[j].indel_len;
							for (k=0; k<in_no; k++){
								lstr[l++] = indeltab.s[j].str[k];
							}
							loopmk = 1;
						}else if (loopmk) break;
					}
					lstr[l] ='\0';
					
					fprintf(fp1,"%s\t%d\t%d\n",refname,pos,pos+1);
					fprintf(fp2,"%s\t%d\t%d\t%s\t%c\t%s\t%s\t%c\n",refname,pos,pos+1,type,lstr[0],lstr,h1_h2,gt);
					
				}
					break;		

			 default:break;
						
			 }
		}else if(hap1[i]==hap2[i]){//BB
			gt='2';					
			if((hap1[i]&mutmsk) == SUBSTITUTE){//subsitute
					strcpy(type,"SNP");
					strcpy(h1_h2,"1/1");
                    var=num_to_gch(hap1[i]&0x3);
					fprintf(fp1,"%s\t%d\t%d\n",refname,pos,pos);
					fprintf(fp2,"%s\t%d\t%d\t%s\t%c\t%c\t%s\t%c\n",refname,pos,pos,type,seq->s[pos],var,h1_h2,gt);
			}else if((hap1[i]&mutmsk)==INSERT && (hap1[i]&LONGINSERT)!=LONGINSERT ){//insertion
				    lstr[0] = num_to_gch( hap1[i]&0x3 );
					strcpy(type,"INS");
					strcpy(h1_h2,"1/1");
					is = (hap1[i]>>28);
					tmp = ( (hap1[i]<<4)>>12 );
					for(j=0;j<is;j++){ lstr[is-j] = num_to_gch( tmp&0x3 ); tmp=tmp>>2; }
					lstr[is]='\0';
					fprintf(fp1,"%s\t%d\t%d\n",refname,pos,pos+1);
					fprintf(fp2,"%s\t%d\t%d\t%s\t%c\t%s\t%s\t%c\n",refname,pos,pos+1,type,lstr[0],lstr,h1_h2,gt);
			}else if((hap1[i]&mutmsk)==DELETE){//deletion
					strcpy(type,"DEL");
					strcpy(h1_h2,"1/1");
					lstr[0] = seq->s[pos-1];
					lstr[1] = num_to_gch( hap1[i]&0x3 );
					tmp=2;
					for(j=i+1; j<sim_len; j++){
						if((hap1[j]&mutmsk)==DELETE && (hap2[j]&mutmsk)!=DELETE && label[j]==label[j-1]+1 ) lstr[tmp++]=num_to_gch(hap1[j]&0x3);
						else break;
					}
					lstr[tmp]='\0';
					i = j-1;
					fprintf(fp1,"%s\t%d\t%d\n",refname,pos,pos+tmp);
					fprintf(fp2,"%s\t%d\t%d\t%s\t%s\t%c\t%s\t%c\n",refname,pos,pos+tmp,type,lstr,lstr[0] ,h1_h2,gt);
			}else if ( (hap1[i]&LONGINSERT) == LONGINSERT ){
					strcpy(type,"LONGINS");
					strcpy(h1_h2,"1/1");				
					int in_no, no, j, k, l,cc, loopmk;
					lstr[0] = num_to_gch( hap1[i]&0x3 );
					l = 1;
					no = ((hap1[i]<<4)>>12);
					loopmk = 0;
					for(j=0; j<indeltab.l; j++){
						if (indeltab.s[j].indel_no == no){
							in_no = indeltab.s[j].indel_len;
							for (k=0; k<in_no; k++){
								lstr[l++] = indeltab.s[j].str[k];
							}
							loopmk = 1;
						}else if (loopmk) break;
					}
					lstr[l] ='\0';
					fprintf(fp1,"%s\t%d\t%d\n",refname,pos,pos+1);
					fprintf(fp2,"%s\t%d\t%d\t%s\t%c\t%s\t%s\t%c\n",refname,pos,pos+1,type,lstr[0] ,lstr,h1_h2,gt);
			}
		}	
	}
	delete [] lstr;
	return ;

}

/**********************************************/
int main(int argc, char *argv[])
{	
	int i,k;
	seq_t	 seq;
	char ref_name[50];
	int ret_eof;
	char ch0[100],fhead[50];
	char bedanno[50];

	if(argc<4){printf("Error!Missing parameters...\nusage: stv *.fa *.sim *.bed\n");return 0;}
	ingenfp  = FileOpen( argv[1],  "r"); 
	insimfp = FileOpen( argv[2],   "r");
	outbedfp = FileOpen( argv[3], "w");
	//outbedannofp = fopen(bedanno,"w");
	
	strcpy(ref_name,argv[1]);
	for (i=0; i<strlen(ref_name); i++) if (ref_name[i] == '.') break; 
	if (i<strlen(ref_name) ) { strncpy(fhead, &ref_name[0], i); fhead[i]='\0'; }
	else strcpy(fhead, ref_name); 

	strcpy(ref_name,argv[2]);
	for (i=0; i<strlen(ref_name); i++) if (ref_name[i] == '.') break; 
	if (i<strlen(ref_name) ) { strncpy(bedanno, &ref_name[0], i); bedanno[i]='\0'; }
	else strcpy(bedanno, ref_name);
	sprintf(bedanno,"%s.anno.bed",bedanno);
	outbedannofp = fopen(bedanno,"w");

	INIT_SEQ(seq);

	//* read reference fasta */	
	REF_LEN = seq_read_fasta(ingenfp, &seq, ref_name, 0);
	cout<<"seq_read_fasta() return! REF_LEN = "<<REF_LEN<<endl;
  		
    label= new int[REF_LEN];
    hap1=new mut_t[REF_LEN];
    hap2=new mut_t[REF_LEN];
	cout<<"begin sim read...."<<endl;		
			
    ret_eof=fscanf(insimfp, "%s\n", ch0);

  	k=0;
	while( ret_eof != EOF ) {
		ret_eof=fscanf(insimfp, "%d %u %u\n", &label[k] ,&hap1[k],&hap2[k]);
		k++;
	 }
	sim_len = k-1;
	label[sim_len]='\0';
	hap1[sim_len]='\0';
    hap2[sim_len]='\0';
	cout<<"sim_len = "<<sim_len<<endl;
	fclose(insimfp);

	read_longinsert_tab( argv[2] );
	
	print_bedfile( outbedfp, outbedannofp, &seq, ref_name );
	
	CLEAR_SEQ(seq);
		
    delete []label;
	delete []hap1;
	delete []hap2;
	fclose(ingenfp); 
	fclose(outbedfp);
	fclose(outbedannofp);
	return 0;
}

