
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
FILE* outvcffp;


#define INIT_SEQ(seq) (seq).s = 0; (seq).l = (seq).max = 0
#define CLEAR_SEQ(seq) free( (seq).s )
#define INIT_MUTSEQ(mseq) (mseq).s = 0; (mseq).l = (mseq).max = 0; (mseq).gnum = (mseq).vnum = 0
#define CLEAR_MUTSEQ(mseq) free( (mseq).s )

enum muttype_t {NOCHANGE = 0, INSERT = 0x80, LONGINSERT = 0xf0000000, SUBSTITUTE = 0x40, DELETE = 0xc0, GTYPEBB=0x20, GTYPEAB=0x10};
typedef unsigned int mut_t;
static mut_t mutmsk = (mut_t)0xc0;
static int SEQ_BLOCK_SIZE=512;

typedef struct {
	int l, max;
	char *s;
} seq_t;

typedef struct {
	int l, max, gnum, vnum; //gnum---valid char so as A C G T, vnum---has varied num
	unsigned int *s;
} mutseq_t;


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



/**********************************************/

void print_vcffile( FILE *fp ,seq_t *seq,  char *refname){
 
    int i,j,pos, tmp;
	char alt;
	char gt[10];
    int ac,ad1,ad2,dp;
	float af;
	int is;//insert size
	char *lstr; //[REF_LEN];
  	int n1=0,n2=0,n3=0,n4=0,n5=0;
	int mark;
    
	lstr = new char [REF_LEN];
	fprintf(fp,"##fileformat=VCFv4.1\n");
	fprintf(fp,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n");
	fprintf(fp,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n");
	fprintf(fp,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n");
	fprintf(fp,"##INFO=<ID=RU,Number=1,Type=String,Description=\"Tandem repeat unit (bases)\">\n");
	fprintf(fp,"##INFO=<ID=STR,Number=1,Type=Flag,Description=\"Variant is a short tandem repeat\">\n");
	fprintf(fp,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n");
	fprintf(fp,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(fp,"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n");
	fprintf(fp,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n");
	fprintf(fp,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n");

//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO							FORMAT	Nhap12878
//  9	    ***	.	*	*	.	  	.		AC=*;AF=*;AN=2;DP=1;RU=*;STR	GT:AD:DP	1/1:0,1:19	***	.	*	*	.	.	AC=*;AF=*;AN=2;DP=1;RU=*;STR	GT:AD:DP	1/1:0,1:1	

	for(i=0; i<sim_len; i++){
		pos=label[i];mark=0;
		if((hap1[i]&mutmsk)==SUBSTITUTE) mark=0;
		if((hap2[i]&mutmsk)==SUBSTITUTE) mark=1;
		if((hap1[i]&mutmsk)==INSERT && (hap1[i]&LONGINSERT)!=LONGINSERT ) mark=2;
		if((hap2[i]&mutmsk)==INSERT && (hap2[i]&LONGINSERT)!=LONGINSERT ) mark=3;
        if((hap1[i]&mutmsk)==DELETE) mark=4;
		if((hap2[i]&mutmsk)==DELETE) mark=5;
		if((hap1[i]&LONGINSERT)==LONGINSERT ) mark=6;
		if((hap2[i]&LONGINSERT)==LONGINSERT ) mark=7;
		if(hap1[i]!= hap2[i] ){//AB
			strcpy(gt,"0/1");
			ac=1;
			af=0.50;
			ad1=1;
			ad2=1;
 			dp=2;

			switch(mark){
		    case 0:{//hap1 subsitute
            		alt=num_to_gch(hap1[i]&0x3);
					fprintf(fp,"%s\t%d\t.\t%c\t%c\t.\t.\tAC=%d;AF=%f;AN=2\tGT:AD:DP\t%s:%d,%d:%d\n",refname,pos,seq->s[pos],
										alt,ac,af,gt,ad1,ad2,dp);
					 break;
			}
			case 1: //hap2 subsitute
					alt=num_to_gch(hap2[i]&0x3);
					fprintf(fp,"%s\t%d\t.\t%c\t%c\t.\t.\tAC=%d;AF=%f;AN=2\tGT:AD:DP\t%s:%d,%d:%d\n",refname,pos,seq->s[pos],
									alt,ac,af,gt,ad1,ad2,dp);
					 break;
			case 2: //hap1 insertion
					lstr[0] = num_to_gch( hap1[i]&0x3 );
					is = (hap1[i]>>28);
					tmp = ( (hap1[i]<<4)>>12 );
					for(j=0;j<is;j++){ lstr[is-j] = num_to_gch( tmp&0x3 ); tmp=tmp>>2; }
					lstr[is+1]='\0';
					fprintf(fp,"%s\t%d\t.\t%c\t%s\t.\t.\tAC=%d;AF=%f;AN=2\tGT:AD:DP\t%s:%d,%d:%d\n",refname,pos,lstr[0],
									lstr,ac,af,gt,ad1,ad2,dp);
					break;
			case 3: //hap2 insertion
					lstr[0] = num_to_gch( hap2[i]&0x3 );
					is = (hap2[i]>>28);
					tmp = ((hap2[i]<<4)>>12);
					for(j=0;j<is;j++){ lstr[is-j] = num_to_gch( tmp&0x3 ); tmp=tmp>>2; }
					lstr[is+1]='\0';
					fprintf(fp,"%s\t%d\t.\t%c\t%s\t.\t.\tAC=%d;AF=%f;AN=2\tGT:AD:DP\t%s:%d,%d:%d\n",refname,pos,lstr[0],
									lstr,ac,af,gt,ad1,ad2,dp);
					break;
			case 4://hap1 deletion
					lstr[0] = seq->s[pos-1];
					lstr[1] = num_to_gch( hap1[i]&0x3 );
					tmp=2;
					for(j=i+1; j<sim_len; j++){
						if((hap1[j]&mutmsk)==DELETE && (hap2[j]&mutmsk)!=DELETE && label[j]==label[j-1]+1 ) lstr[tmp++]=num_to_gch(hap1[j]&0x3);
						else break;
					}
					lstr[tmp]='\0';
					i = j-1;
					fprintf(fp,"%s\t%d\t.\t%s\t%c\t.\t.\tAC=%d;AF=%f;AN=2\tGT:AD:DP\t%s:%d,%d:%d\n",refname,pos,lstr,
									lstr[0],ac,af,gt,ad1,ad2,dp);
					break;
			case 5://hap2 deleteion
					lstr[0] = seq->s[pos-1];
					lstr[1] = num_to_gch( hap2[i]&0x3 );
					tmp=2;
					for(j=i+1; j<sim_len; j++){
						if((hap2[j]&mutmsk)==DELETE && (hap1[j]&mutmsk)!=DELETE && label[j]==label[j-1]+1 ) lstr[tmp++]=num_to_gch(hap2[j]&0x3);
						else break;
					}
					lstr[tmp]='\0';
					i = j-1;
					fprintf(fp,"%s\t%d\t.\t%s\t%c\t.\t.\tAC=%d;AF=%f;AN=2\tGT:AD:DP\t%s:%d,%d:%d\n",refname,pos,lstr,
								lstr[0],ac,af,gt,ad1,ad2,dp);
					break;
			case 6://hap1 long insertion
				{
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
					fprintf(fp,"%s\t%d\t.\t%c\t%s\t.\t.\tAC=%d;AF=%f;AN=2\tGT:AD:DP\t%s:%d,%d:%d\n",refname,pos,lstr[0],
									lstr,ac,af,gt,ad1,ad2,dp);
				}
					break;

			case 7://hap2 long insertion
				{
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
					fprintf(fp,"%s\t%d\t.\t%c\t%s\t.\t.\tAC=%d;AF=%f;AN=2\tGT:AD:DP\t%s:%d,%d:%d\n",refname,pos,lstr[0],
									lstr,ac,af,gt,ad1,ad2,dp);
				}
					break;
            default:break;
			}
		} else if(hap1[i] == hap2[i]){//BB
			strcpy(gt,"1/1");
			ac=2;
			af=1.00;
			ad1=0;
			ad2=2;
			dp=2;
			if((hap1[i]&mutmsk) == SUBSTITUTE){//subsitute
				alt=num_to_gch(hap1[i]&0x3);
				fprintf(fp,"%s\t%d\t.\t%c\t%c\t.\t.\tAC=%d;AF=%f;AN=2\tGT:AD:DP\t%s:%d,%d:%d\n",refname,pos,seq->s[pos],
								alt,ac,af,gt,ad1,ad2,dp);
			}else if((hap1[i]&mutmsk)==INSERT && (hap1[i]&LONGINSERT)!=LONGINSERT ){//insertion
				lstr[0] = num_to_gch( hap1[i]&0x3 );
				is = (hap1[i]>>28);
				tmp = ((hap1[i]<<4)>>12);
				for(j=0;j<is;j++){ lstr[is-j] = num_to_gch( tmp&0x3 ); tmp=tmp>>2; }
				lstr[is+1]='\0';
				fprintf(fp,"%s\t%d\t.\t%c\t%s\t.\t.\tAC=%d;AF=%f;AN=2\tGT:AD:DP\t%s:%d,%d:%d\n",refname,pos,lstr[0],
							lstr,ac,af,gt,ad1,ad2,dp);
			}else if((hap1[i]&mutmsk)==DELETE){//deletion
				lstr[0] = seq->s[pos-1];
				lstr[1] = num_to_gch( hap1[i]&0x3 );
				tmp=2;
				for(j=i+1; j<sim_len; j++){
					if((hap1[j]&mutmsk)==DELETE && (hap2[j]&mutmsk)==DELETE && label[j]==label[j-1]+1 ) lstr[tmp++]=num_to_gch(hap1[j]&0x3);
					else break;
				}
				lstr[tmp]='\0';
				i = j-1;
				fprintf(fp,"%s\t%d\t.\t%s\t%c\t.\t.\tAC=%d;AF=%f;AN=2\tGT:AD:DP\t%s:%d,%d:%d\n",refname,pos,lstr,
							lstr[0],ac,af,gt,ad1,ad2,dp);
			}else if ( (hap1[i]&LONGINSERT) == LONGINSERT ){
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
					fprintf(fp,"%s\t%d\t.\t%c\t%s\t.\t.\tAC=%d;AF=%f;AN=2\tGT:AD:DP\t%s:%d,%d:%d\n",refname,pos,lstr[0],
									lstr,ac,af,gt,ad1,ad2,dp);
			}
		}	
	}
	delete [] lstr;
	return ;
}

///*****************
void put_gch(seq_t *s0, seq_t *s1, int numch){
	int max;
	max = s0->max;
	if (s0->l + 1 >= max) {
		max += SEQ_BLOCK_SIZE;
		s0->s = (char *)realloc(s0->s, sizeof(char) * max);
		s1->s = (char *)realloc(s1->s, sizeof(char) * max);
		s0->max = max;
		s1->max = max;
	}
	s0->s[s0->l++] = num_to_gch(numch & 0xf);
	s1->s[s1->l++] = gch_rev( num_to_gch(numch & 0xf) );
	return;
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
	printf("END read long insert table::  Long insert num=%d\n",n);
	printf("Long insert table size=%d\n", indeltab.l);
	delete [] idstr;
	fclose(fp);
	return 0;
}

/**********************************************/
int main(int argc, char *argv[])
{	
	int i;
	seq_t	 seq;
	mutseq_t sim[2];
	char ref_name[50];
	int ret_eof;
	char ch0[100],fhead[50];
	

	if(argc<4){printf("Error!Missing parameters...\n\tusage: stv *.fa *.sim *.vcf\n");return 0;}
	ingenfp  = FileOpen( argv[1], "r"); 
	insimfp  = FileOpen( argv[2], "r");
	outvcffp = FileOpen( argv[3], "w");

	
	strcpy(ref_name,argv[1]);
	for (i=0; i<strlen(ref_name); i++) if (ref_name[i] == '.') break; 
	if (i<strlen(ref_name) ) { strncpy(fhead, &ref_name[0], i); fhead[i]='\0'; }
	else strcpy(fhead, ref_name); 

	INIT_SEQ(seq);
	INIT_SEQ(indeltab);

	//* read reference fasta */	
	REF_LEN = seq_read_fasta(ingenfp, &seq, ref_name, 0);
	cout<<"seq_read_fasta() return! REF_LEN = "<<REF_LEN<<endl;
 		
    label= new int[REF_LEN];
    hap1 =new mut_t[REF_LEN];
    hap2 =new mut_t[REF_LEN];
	cout<<"begin sim read...."<<endl;		
			
    ret_eof=fscanf(insimfp, "%s\n", ch0);
  	int k=0;
	while( fscanf(insimfp, "%d %u %u\n", &label[k], &hap1[k], &hap2[k]) != EOF ) {
		k++;
	}
	sim_len = k;
	cout<<"sim_len = "<<sim_len<<endl;
	fclose(insimfp);
	
	read_longinsert_tab( argv[2] );

	print_vcffile( outvcffp,&seq, ref_name );

	CLEAR_SEQ(seq);
	for (i=0; i<indeltab.l; i++) if (indeltab.s[i].str) free(indeltab.s[i].str);
	CLEAR_SEQ(indeltab);

    delete []label;
	delete []hap1;
	delete []hap2;
	fclose(ingenfp); 
	fclose(outvcffp);
	return 0;
}

