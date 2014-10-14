
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


/*************************************************************/
int read_line(FILE *fp, char *s){
int k;
char ch;
	k = 0;
	while(!feof(fp) && (ch=fgetc(fp)) != '\n') s[k++] = ch;
	s[k]='\0';
	if (feof(fp)) return -1;
	else return k;
}


/*******************************************************/
int splitline(char *s, char *format[20], char splitmark){
int j, k, n, slen;
	j=0; n = 0;
	slen = strlen(s);
	while (j < slen){
		while((s[j]==splitmark || s[j]==' ')  && j<slen ) j++;
		k = 0;
		while((s[j]!= splitmark && s[j] !=' ') && j<slen ) { format[n][k++] = s[j++]; }
		format[n][k]='\0';
		n++;
	}
	return n;    
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


/********************************************************************/

int main(int argc, char *argv[]){	

	int i,k,n, line, ret_eof;
	mut_t tmp,is, tt;
	seq_t	 seq;
	mutseq_t sim[2];
	char ref_name[50], ch0[100],fhead[50];
	char str[1000];
	char *term[20];
	int pos,rand,gt_pos,num0,num8, num9, longin_no, long_indel_baseno;
	int ab_num, bb_num, sdel_num, sin_num, ldel_num, lin_num, subs_num;
	char gt[10];
	int 	sim_len;
	char 	*format[20];
	char 	*strtmp;
	
/*
//	tt=1077461665;
	tt=1950417058;
	i=(tt>>28);
	k=(tt<<4)>>12;
	n=tt&0xff;
	printf("g4=%d   mid=%d:%d:%d:%d  g7-4=%d   g3-0=%d", i, k&0x3,(k>>2)&0x3, (k>>4)&0x3, (k>>6)&0x3, (n>>4)&0xf, n&0xf);
	return 0;
*/
	FILE* ingenfp;
	FILE* invcffp;
	FILE* outsimfp;
	FILE *outlonginfp;

	if(argc<5){printf("Error!Missing parameters...\n\tusage: vcftosim  *.fa  *.vcf  *.sim long_indel_base_no\n");return 0;}
	ingenfp  = FileOpen( argv[1], "r"); 
	invcffp  = FileOpen( argv[2], "r");
	outsimfp = FileOpen( argv[3], "w");
	long_indel_baseno = atoi(argv[4]);
	sprintf(fhead, "%s.idx", argv[3]);
	outlonginfp = FileOpen( fhead, "w");

	INIT_SEQ(seq);
	INIT_MUTSEQ(sim[0]);
	INIT_MUTSEQ(sim[1]);
	strtmp = new char [200000];
	for(i=0; i<20;i++) term[i] = new char [200000];	
	for(i=0; i<20;i++) format[i] = new char [200000];	

	//* read reference fasta */	
	REF_LEN = seq_read_fasta(ingenfp, &seq, ref_name, 0);
	cout<<"seq_read_fasta() return! REF_LEN = "<<REF_LEN<<endl;

	fprintf(outsimfp,"&%s\n", ref_name);
	//initial sim
	sim[0].l = seq.l; sim[1].l = seq.l;
	sim[0].max = seq.max; 
	sim[1].max = seq.max;
	sim[0].s = (mut_t *)calloc(seq.max, sizeof(mut_t));
	sim[1].s = (mut_t *)calloc(seq.max, sizeof(mut_t));

	for (i = 0; i < seq.l; ++i) {
		sim[0].s[i] = sim[1].s[i] = (mut_t)gch_to_num(seq.s[i]);
	}
	//process vcf per line
    line=0;
	while(!feof(invcffp) && (read_line(invcffp, strtmp) >= 0) && strtmp[0] == '#') line++;
 
//	fseek(invcffp,-1,SEEK_CUR);
	ab_num = bb_num = sdel_num = sin_num = ldel_num = lin_num = subs_num = 0;
	longin_no = 0;
	while(!feof(invcffp) ){
		// split one line from *.vcf file
  	 	num0 = splitline(strtmp, term, '\t');
		//find GT
		gt_pos=-1;
		num8=splitline(term[8], format, ':');
		for(i=0;i<num8;i++) { 
			if (!strcmp(format[i],"GT")) {gt_pos=i; break;}
		}
		if (gt_pos<0) printf("gt_pos ERROR????  num=%d,  gt_pos=%d  term[8]=%s\n", num8, gt_pos, term[8]);
		num9=splitline(term[9], format, ':');
		if (gt_pos<num9) {strcpy(gt, format[gt_pos]);}      

		pos=atoi(term[1]);
		k = strlen(term[3])>strlen(term[4])? strlen(term[3]): strlen(term[4]);
		if (pos+k>=REF_LEN) break;

		if(!strcmp(gt,"0/1")){//genotype=AB
			ab_num++;
			if((strlen(term[3])==1) && (strlen(term[4])==1)){//subsitution
				rand = drand48()<0.5?0:1;	
				sim[rand].s[pos]=  ( ( (sim[rand].s[pos]>>4)<<4 ) | gch_to_num(term[4][0]) | SUBSTITUTE | GTYPEAB );
				fprintf(outsimfp, "%10u\t%10u\t%10u\n", pos, sim[0].s[pos], sim[1].s[pos]);
				subs_num++;
			}else if(strlen(term[3])>strlen(term[4])){//deletion
				rand = drand48()<0.5?0:1;
				is=strlen(term[3])-1;
				for(i=0;i<is;i++){
					sim[rand].s[pos+i]=sim[rand].s[pos+i] | DELETE |GTYPEAB;
					fprintf(outsimfp, "%10u\t%10u\t%10u\n", pos+i, sim[0].s[pos+i], sim[1].s[pos+i]);
				}
				if (is<10) sdel_num++; else ldel_num++; 
			}else if((strlen(term[3])<strlen(term[4])) && (strlen(term[4])<=10)){// short insertion
				rand = drand48()<0.5?0:1;
				is=strlen(term[4])-1;
				tmp = 0;
				for(i=0;i<is;i++) tmp = (tmp<<2) | gch_to_num(term[4][i+1]);
				sim[rand].s[pos]= sim[rand].s[pos] | (is<<28) | (tmp<<8) | INSERT | GTYPEAB;
				fprintf(outsimfp, "%10u\t%10u\t%10u\n", pos, sim[0].s[pos], sim[1].s[pos]);
				sin_num++;
			}else if((strlen(term[3])<strlen(term[4])) && (strlen(term[4])>10)){// long insertion
				rand = drand48()<0.5?0:1;
				is=strlen(term[4])-1;
				tmp = long_indel_baseno + longin_no;
				sim[rand].s[pos]= sim[rand].s[pos] | (tmp<<8) | LONGINSERT | INSERT | GTYPEAB;
				fprintf(outsimfp, "%10u\t%10u\t%10u\n", pos, sim[0].s[pos], sim[1].s[pos]);
				fprintf(outlonginfp, "%d %d %d %d %d %d %s\n", tmp, pos, is, 1, 1, 0, &term[4][1]);
				longin_no++; lin_num++;
			}
		}else if(!strcmp(gt,"1/1")){//genotype=BB
			bb_num++;
			if((strlen(term[3])==1) && (strlen(term[4])==1)){//subsitution
					sim[0].s[pos] = sim[1].s[pos] = ((sim[0].s[pos]>>4)<<4)|gch_to_num(term[4][0])|SUBSTITUTE|GTYPEBB;
					fprintf(outsimfp, "%10u\t%10u\t%10u\n", pos, sim[0].s[pos], sim[1].s[pos]);
					subs_num++;
			}else if(strlen(term[3])>strlen(term[4])){//deletion
					is=strlen(term[3])-1;
					for(i=0;i<is;i++){
						sim[0].s[pos+i] = sim[1].s[pos+i] = sim[0].s[pos+i] | DELETE |GTYPEBB;
						fprintf(outsimfp, "%10u\t%10u\t%10u\n", pos+i, sim[0].s[pos+i], sim[1].s[pos+i]);
					}
					if (is<10) sdel_num++; else ldel_num++; 
			}else if((strlen(term[3])<strlen(term[4]))&& (strlen(term[4])<=11)){//short insertion
				is=strlen(term[4])-1;
				tmp = 0;
				for(i=0;i<is;i++) tmp = (tmp<<2) | gch_to_num(term[4][1+i]) ;
				sim[0].s[pos]=sim[1].s[pos]= sim[0].s[pos] | (is<<28) | (tmp<<8) | INSERT | GTYPEBB  ;
				fprintf(outsimfp, "%10u\t%10u\t%10u\n", pos, sim[0].s[pos], sim[1].s[pos]);
				sin_num++;
			}else if((strlen(term[3])<strlen(term[4])) && (strlen(term[4])>10)){// long insertion
				is=strlen(term[4])-1;
				tmp = long_indel_baseno + longin_no;
				sim[0].s[pos]= sim[1].s[pos]= sim[0].s[pos] | (tmp<<8) | LONGINSERT | INSERT | GTYPEAB;
				fprintf(outsimfp, "%10u\t%10u\t%10u\n", pos, sim[0].s[pos], sim[1].s[pos]);
				fprintf(outlonginfp, "%d %d %d %d %d %d %s\n", tmp, pos, is, 2, 1, 0, &term[4][1]);
				longin_no++; lin_num++;
			}				
		} //end of genotype=BB
		if (line%1000==0) printf("line=%d\n",line);
		line++;
		read_line(invcffp, strtmp);
 	}	
	printf("ab_num = %d\tbb_num = %d\tAB+BB(include long INDEL)=%d\t\n",ab_num, bb_num, ab_num+bb_num);
	printf("sdel_num = %d\tsin_num = %d\tsdel_num+sin_num=%d\t\n", sdel_num, sin_num, sdel_num+sin_num);
	printf("ldel_num = %d\tlin_num = %d\tldel_num+lin_num=%d\nsubsitution var_num = %d\n\n",ldel_num, lin_num, ldel_num+lin_num, subs_num);

	CLEAR_SEQ(seq);
	if (sim[0].s) CLEAR_MUTSEQ(sim[0]);
	if (sim[1].s) CLEAR_MUTSEQ(sim[1]);
	delete [] strtmp;
	for(i=0; i<20;i++) delete [] term[i];
	for(i=0; i<20;i++) delete [] format[i];

	fclose(ingenfp); 
	fclose(invcffp);
	fclose(outsimfp);
	fclose(outlonginfp);
	return 0;

}

