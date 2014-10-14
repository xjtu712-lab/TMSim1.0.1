/*readgen--- the program of reads generating for single or paired_end file(*_left.fq, *_right.fq)
 *Based on *.sim files that are generated from normal and subclones. 
 *It form the simulation environment including all kinds of variations and complete reads random sampling. 
 *Last, it can random generate two files of left.fq and right.fq that exit in positive and reverse chain,respectively.*/
 
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

int	MAX_DIS=10000;
int	DALT=10;
int	LREAD_LEN=100;
int	RREAD_LEN=100;
float 	COVER_RATE = 5.0;
float 	ERROR_RATE=0.00;
int 	KEEP_N = 0;
int 	SINGLE_READ = 0;

int 	LONG_INDEL = 0; 
int 	INDEL_NUM = 0;
int	REF_LEN=120000000;
char ChrName[20];

FILE* ingen; 
FILE* insim; 
FILE* outleft;
FILE* outright;
FILE* outresult;

int mainsub();
void initdata();
static void pair_read_generator(int l_beg, int r_end, char *left, char *right);

#define INIT_SEQ(seq) (seq).s = 0; (seq).l = (seq).max = 0
#define CLEAR_SEQ(seq) free( (seq).s )

enum muttype_t {NOCHANGE = 0, INSERT = 0x80, LONGINSERT = 0xf0000000, SUBSTITUTE = 0x40, DELETE = 0xc0, GTYPEBB=0x20, GTYPEAB=0x10};
typedef unsigned int mut_t;
static mut_t mutmsk = (mut_t)0xc0;
static int SEQ_BLOCK_SIZE=512;

typedef struct {
	int l, max;
	char *s;
} seq_t;

typedef struct {
	int l,max;
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
int MAXSS = pow(2, 31)-1;
static double zrand(int max){
	static double r;
	int n;
	n = rand() % max;
	r = 1.0*n/max;
	return r;
}


/**********************************************/
static double zt_stdfun(){
	static double v1,v2,s;
	do {
		v1 = zrand(MAXSS);
		v2 = zrand(MAXSS);
		s = 1.0*cos(2*3.1415926*v1)*sqrt(-2*log(v2));
	} while(s>2 || s<-2);
	return s;
}

/**********************************************/
FILE *FileOpen(char *fileName, const char *mode)
{
	FILE *fp=NULL;
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
	int c, l, max;
	char *p;
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
	l = 0; max = seq->max;
	while (!feof(fp) && (c = fgetc(fp)) ) {
		if (isalpha(c) || c == '-' || c == '.') {
			if (l + 1 >= max) {
				max += SEQ_BLOCK_SIZE;
				seq->s = (char *)realloc(seq->s, sizeof(char) * max);
			}
			seq->s[l++] = (char)c;
		}
	}
	seq->s[l] = 0;
	seq->max = max; 
	seq->l = l;
	fprintf(outresult, "\nseq_read_fasta():\nref_name=[%s] From *.fa file.\n",chrname);
	fprintf(outresult, "ref_strlen=[%d] seq->max = [%d]\n",seq->l, seq->max);
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

/*******************************************************/
void init_sim(seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, char * indelfile, indel_t *indelp)
{
	int fret;
	int no, idno, idbeg, idend, idlen, idtype, idnum, idvar;
	int pos, n, m, i, l;
	char name[50], *idstr, ch=0;
	mut_t s0, s1;
	mutseq_t *ret[2];
	FILE *indelfp;

	ret[0] = hap1; ret[1] = hap2;
	ret[0]->l = seq->l; ret[1]->l = seq->l;
	ret[0]->max = seq->max; ret[1]->max = seq->max;
	ret[0]->s = (mut_t *)calloc(seq->max, sizeof(mut_t));
	ret[1]->s = (mut_t *)calloc(seq->max, sizeof(mut_t));
	m = 0;
	for (i = 0; i < seq->l; ++i) {
		ret[0]->s[i] = ret[1]->s[i] = (mut_t)gch_to_num(seq->s[i]);
		if ( (ret[0]->s[i]&0xf) < 4) m++;
	}
	n=0;
	fret=fscanf(insim, "%s\n", name);
	fret=fscanf(insim, "%d %u %u\n", &pos, &s0, &s1);
	while (fret != EOF){
		n++;
		ret[0]->s[pos] = s0;
		ret[1]->s[pos] = s1;
		fret=fscanf(insim, "%d %u %u\n", &pos, &s0, &s1);
	}
	fprintf(outresult, "\ninit_sim():\nchanged num=[%d] in valid gch [%d] From *.sim file.\n",n, m);

	// handle for long insert area from long indel file
	if (!LONG_INDEL) return ;
	indelfp = FileOpen( indelfile, "rt");
	if (!indelfp) return;
	idstr = new char [100000];
	n=0; 
	while ( !feof(indelfp) && fscanf(indelfp, "%d %d %d %d %d %d %s\n", &idno, &idbeg, &idlen, &idtype, &idnum, &idvar, idstr) != EOF ){
		addinsert_iterm(indelp, idno, idbeg, idlen, idtype, idnum, idvar, idstr);
		n++;
	} //while (reading file .....)
	printf("END:: Long insert num=%d    Long insert table size=%d\n", n, indelp->l);
	delete [] idstr;
	fclose(indelfp);
	return;
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
void press_sim(mutseq_t *hap1, mutseq_t *hap2, seq_t seq[2][2], indel_t *indelp)
{
	int i, j, k, n, m;
	seq_t  *s[2];
	mut_t  msk, c;

	m=0;
	s[0] = &seq[0][0];
	s[1] = &seq[0][1];
	for (i = 0; i < hap1->l; ++i) {
		c = hap1->s[i];
		msk = c & mutmsk;
		if (msk != DELETE && (c & 0xf) <= (KEEP_N?4:3) ){
			if (msk == 0 || msk == SUBSTITUTE ){ 
				put_gch(s[0], s[1], (c&0xf) ); 
				if (msk == SUBSTITUTE) m++; 
			}else if ((c & LONGINSERT) == LONGINSERT){ // long insert 
				int in_no, pos, k, cc, loopmk;
				put_gch(s[0], s[1], (c&0x3) ); m++;
				pos = ((c<<4)>>12);
				loopmk = 0;
				for(j=0; j<indelp->l; j++){
					if (indelp->s[j].indel_no == pos){
						in_no = indelp->s[j].indel_len;
						for (k=0; k<in_no; k++){
							cc = gch_to_num(indelp->s[j].str[k]);
							if (cc < 4){ put_gch(s[0], s[1],  cc ); m++; }
						}
						loopmk = 1;
					}else if (loopmk) break;
				}
			}else if (msk == INSERT) {  // insert
				int k, innum, inch;
				put_gch(s[0], s[1], (c&0x3) ); m++;
				innum = (c>>28);
				inch  = ((c<<4)>>12)<<(12+(10-innum)*2);
				k=0;
				while (k<innum){
					put_gch(s[0], s[1], (inch>>30)&3 ); m++;
					k++; inch = (inch<<2);
				}
			} // end insert
		} else if (msk == DELETE) m++; 
	}
	printf("Enter press_sim():\n\tref0 len=[%d] ref0' len=[%d] for valid gchar.\n", seq[0][0].l,seq[0][1].l);
	printf("\tref0 snp num=[%d] var_rate=[%f] for valid gchar.\n", m, 1.0*m/seq[0][0].l);
	fprintf(outresult, "\npress_sim():\nref0 len=[%d] ref0' len=[%d] for valid gchar.\n",
						seq[0][0].l,seq[0][1].l);
	fprintf(outresult, "ref0 snp num=[%d] var_rate=[%f] for valid gchar.\n", m, 1.0*m/seq[0][0].l);
	m=0;
	s[0] = &seq[1][0];
	s[1] = &seq[1][1];
	for (i = 0; i < hap2->l; ++i) {
		c = hap2->s[i];
		msk = c & mutmsk;
		if (msk != DELETE && (c & 0xf) <= (KEEP_N?4:3) ){
			if (msk == 0 || msk == SUBSTITUTE ) { 
				put_gch(s[0], s[1], (c&0xf) ); 
				if (msk == SUBSTITUTE) m++; 
			}else if ((c & LONGINSERT) == LONGINSERT){ // long insert 
				int in_no, pos, cc, mark;
				put_gch(s[0], s[1], (c&0x3) ); m++;
				pos = ((c<<4)>>12);
				mark = 0;
				for(j=0; j<indelp->l; j++){
					if (indelp->s[j].indel_no == pos){
						in_no = indelp->s[j].indel_len;
						for (k=0; k<in_no; k++){
							cc = gch_to_num(indelp->s[j].str[k]);
							if (cc < 4){ put_gch(s[0], s[1],  cc ); m++; }
						}
						mark = 1;
					}else if (mark) break;
				}
			}else if (msk == INSERT) {  // insert
				int innum, inch;
				put_gch(s[0], s[1], (c&0x3) ); m++;
				innum = (c>>28);
				inch  = ((c<<4)>>12)<<(12+(10-innum)*2);
				k=0;
				while (k<innum){
					put_gch(s[0], s[1], (inch>>30)&3 ); m++;
					k++; inch = (inch<<2);
				}
			} // end insert
		}else if (msk == DELETE) m++;
	}
	printf("\n\tref1 len=[%d] ref1' len=[%d]  for valid gchar.\n",seq[1][0].l,seq[1][1].l);
	printf("\tref1 snp num=[%d] var_rate=[%f] for valid gchar.\n", m, 1.0*m/seq[1][0].l);

	fprintf(outresult, "\nref1 len=[%d] ref1' len=[%d]  for valid gchar.\n",seq[1][0].l,seq[1][1].l);
	fprintf(outresult, "ref1 snp num=[%d] var_rate=[%f] for valid gchar.\n", m, 1.0*m/seq[1][0].l);
	return;
}

/**********************************************/
#define PACKAGE_VERSION "1.0.1"
static int simu_usage()
{	char ss[200];

	printf("\nProgram: ReadGen (the reads generator from *.sim)\n");
	printf("Version: %s \n", PACKAGE_VERSION);
	printf("Contact: Yu Geng <gengyu@stu.xjtu.edu.cn>;Zhongmeng Zhao <zmzhao@mail.xjtu.edu.cn>\n");
	printf("Usage:   readgen [options] <ref.fa> <*.sim> <left.fq> <right.fq>\n\n");
	printf("Options: -d INT		outer distance between the two ends of paired_end reads [%d]\n", MAX_DIS);
	printf("         -s INT     standard deviation of MAX_DIS [%d]\n", DALT);;
	printf("         -l INT     length of left  read [%d]\n", LREAD_LEN);
	printf("         -r INT     length of right read [%d]\n", RREAD_LEN);
	printf("         -c float   cover rate of read [%f]\n", COVER_RATE);
	printf("         -e float   error rate in generate reads [%f]\n", ERROR_RATE);
	printf("         -S       	output single read(not paired_end reads)\n");
	printf("         -k       	keep 'N' character in reads\n");
	printf("         -I <subclonde.sim.idx>   long insert table file from 'TumSim'\n");
	printf("         -O <genread.log>  		  result for runing case\n");
	return 1;
}

/**********************************************/

int gen_pair_reads(seq_t qseq[2][2], char *ref_name);
int gen_single_read(seq_t qseq[2][2], char *ref_name);

/**********************************************/
int main(int argc, char *argv[])
{	char c, cc;
	int  i;
	char ingenf[50];
	char insimf[50];
	char indelf[50];
	char outleftf[50];
	char outrightf[50];
	char outresultf[50];
	seq_t	 seq, gseq[2][2];
	mutseq_t sim[2];
	indel_t  indelp;

	srand(unsigned (time(0)));  // Init Random data 
	srand48(time(0));
	indelf[0]='\0'; outrightf[0]='\0'; outresultf[0]='\0';
	cc=0;
	while ((c = getopt(argc, argv, "d:s:l:r:c:e:I:O:n:Sk")) >= 0) {
		switch (c) {
		case 'd': MAX_DIS = atoi(optarg); if (MAX_DIS<LREAD_LEN+RREAD_LEN) cc=1; break;
		case 's': DALT 	  = atoi(optarg); if (DALT<0||DALT>=MAX_DIS/2) cc=1; break;
		case 'l': LREAD_LEN = atoi(optarg); if (LREAD_LEN<=0 ) cc=1; break;
		case 'r': RREAD_LEN = atoi(optarg); if (RREAD_LEN<=0 ) cc=1; break;
		case 'c': COVER_RATE = atof(optarg); if (COVER_RATE<=0) cc=1; break;
		case 'e': ERROR_RATE = atof(optarg); if (ERROR_RATE<0 || ERROR_RATE > 1) cc=1; break;
		case 'S': SINGLE_READ = 1;  break;
		case 'k': KEEP_N = 1;  break;
		case 'I': LONG_INDEL=1; strcpy(indelf, optarg);break;
		case 'O': strcpy(outresultf, optarg); break;
		default : break;
		}
		if (cc) {optind = argc; break;}
	}
	if (optind != (argc-(SINGLE_READ?3:4)) ) return simu_usage();
	strcpy(ingenf, argv[optind]);
	strcpy(insimf, argv[optind+1]);
	strcpy(outleftf, argv[optind+2]);
	if (!SINGLE_READ) strcpy(outrightf, argv[optind+3]);
	if (strlen(outresultf)==0) sprintf(outresultf, "genread.log"); 
	ingen    = FileOpen( ingenf, "rt"); 
	insim    = FileOpen( insimf, "rt");
	outleft  = FileOpen( outleftf,  "wt");
	if (!SINGLE_READ) outright = FileOpen( outrightf, "wt");
	outresult= FileOpen( outresultf, "wt");
	printf("\nmain():\n\tInput file ingenf=[%s] insimf=[%s]\n", ingenf, insimf);
	printf("\t\tOutput file outleftf=[%s] outrightf=[%s] outresultf=[%s]\n", outleftf,outrightf,outresultf);
	fprintf(outresult, "\nmain():\nInput file ingenf=[%s] insimf=[%s]\n", ingenf, insimf);
	fprintf(outresult, "       \nOutput file outleftf=[%s] outrightf=[%s] outresultf=[%s]\n", 
					outleftf,outrightf,outresultf);
	if (!SINGLE_READ) printf("paired-end reads!\n\n");

	INIT_SEQ(seq);
	INIT_SEQ(sim[0]);
	INIT_SEQ(sim[1]);
	INIT_SEQ(gseq[0][0]);
	INIT_SEQ(gseq[0][1]);
	INIT_SEQ(gseq[1][0]);
	INIT_SEQ(gseq[1][1]);
	INIT_SEQ(indelp);

	//* read reference fasta */	
	REF_LEN = seq_read_fasta(ingen, &seq, ChrName, 0);
	printf("seq_read_fasta() return! REF_LEN = %d\n\n", REF_LEN);

	//* init sim and set from *.sim file */ 
	init_sim(&seq, &sim[0], &sim[1], indelf, &indelp);
	printf("init_sim() return!\n\n");

	//* Press reference to 4 seq */ 
	press_sim(&sim[0], &sim[1], gseq, &indelp);
	printf("press_sim() return!\n\n");

	//generate pair_end reads
	if (SINGLE_READ) gen_single_read(gseq, ChrName);
	else gen_pair_reads(gseq, ChrName);
	printf("gen_pair_reads() or gen_single_read() return!\n");
	printf("\nEnd main()!\n");

	CLEAR_SEQ(seq);
	CLEAR_SEQ(sim[0]);
	CLEAR_SEQ(sim[1]);
	CLEAR_SEQ(gseq[0][0]);
	CLEAR_SEQ(gseq[0][1]);
	CLEAR_SEQ(gseq[1][0]);
	CLEAR_SEQ(gseq[1][1]);
	for (i=0; i<indelp.l; i++) if (indelp.s[i].str) free(indelp.s[i].str);
	CLEAR_SEQ(indelp);

	fclose(ingen); 
	fclose(insim); 
	fclose(outleft);
	if (!SINGLE_READ) fclose(outright);
	fclose(outresult);
	return 0;
}

/**********************************************/
static int out_read(int l_beg, int r_end, char *l_readmk, char *r_readmk, int pr_no, 
				seq_t qseq[2][2], int n1, int n2, char * ref_name)
{
	char left[LREAD_LEN+1], ch;
	char right[RREAD_LEN+1];
	int j, l_zbeg, r_zbeg, err_num, k, n, total_err, lzfmark, rzfmark;
	float f1; 

	err_num=0; total_err=0;
	for (j=0;j<LREAD_LEN;j++){
		if (zrand(MAXSS)<ERROR_RATE){ err_num++; total_err++; }
		if (n2 == 0) left[j] = qseq[n1][0].s[l_beg+j];
		else left[j] = qseq[n1][1].s[r_end-j];
	}
	left[LREAD_LEN] = '\0';
	if (err_num>0){
		for (j=0; j<err_num; j++){
			k = zt_stdfun()*LREAD_LEN;
			if (k<0) k=-k;
			k = k % LREAD_LEN;
			n = LREAD_LEN-1-k;
			k = gch_to_num(left[n]);
			k = (k + (int)(zrand(MAXSS) * 3.0 + 1)) & 3;
			left[n] = num_to_gch(k);
		}
	}

	err_num=0;
	for (j=0;j<RREAD_LEN;j++){
		if (zrand(MAXSS)<ERROR_RATE) { err_num++; total_err++; }
		if (n2 == 0) right[j] = qseq[n1][1].s[r_end-j];
		else right[j] = qseq[n1][0].s[l_beg+j];
	}
	right[RREAD_LEN]= '\0';
	if (err_num>0){
		for (j=0; j<err_num; j++){
			k = zt_stdfun()*RREAD_LEN;
			if (k<0) k=-k;
			k = k % RREAD_LEN;
			n = RREAD_LEN-1-k;
			k = gch_to_num(right[n]);
			k = (k + (int)(zrand(MAXSS) * 3.0 + 1)) & 3;
			right[n] = num_to_gch(k);
		}
	}

	if (n2 == 0) { l_zbeg = l_beg; r_zbeg = r_end - RREAD_LEN + 1; } 
	else { l_zbeg = r_end - RREAD_LEN + 1; r_zbeg = l_beg; }

	fprintf(outleft, "@%s_%d_%d_%d_%d_%d_0:0:0_0:0:0_%x/1\n", ref_name, l_beg+1, r_end+1, l_zbeg+1, n1, n2, pr_no);
	fprintf(outleft, "%s\n", left);
	fprintf(outleft, "+\n%s\n", l_readmk);

	fprintf(outright,"@%s_%d_%d_%d_%d_%d_0:0:0_0:0:0_%x/2\n", ref_name, l_beg+1, r_end+1, r_zbeg+1, n1, n2, pr_no);	
	fprintf(outright, "%s\n", right);
	fprintf(outright, "+\n%s\n", r_readmk);

	return total_err;
}


/**********************************************/
//msk=0--rand, msk==1---beg, msk==2---end
void  beg_end(int msk, int &beg, int &end, int len){ 
	int pair_reads_dis, k;

	k = zt_stdfun()*DALT;
	while (k>2*DALT || k<(-2*DALT)) k = zt_stdfun()*DALT;
	pair_reads_dis = MAX_DIS+k;

	switch(msk){
	case 1:	end = beg + pair_reads_dis - 1;	break;
	case 2: beg = end - pair_reads_dis + 1; break;
	case 0:	//k = (int)(1.0*zrand(MAXSS)*len);
		while (1){
			beg = zrand(MAXSS)*len;
			if(beg + pair_reads_dis - 1 < len) break; 
		}
		end = beg + pair_reads_dis - 1;
		break;
	default: break;
	}
	return;
}


/**********************************************/
int gen_pair_reads(seq_t qseq[2][2], char *ref_n)
{
	char l_readmk[LREAD_LEN+1];//={'I'};
	char r_readmk[RREAD_LEN+1];//={'I'};
	seq_t *s[2];
	int PAIR_READS_NUM;
	int l_beg, r_end, r_beg, l_len,r_len,a_len;
	int i, k,m, n, pr_no, err_num, err_preads; 

	for (i=0;i<LREAD_LEN;i++) l_readmk[i]='I';
	l_readmk[LREAD_LEN]='\0';
	for (i=0;i<RREAD_LEN;i++) r_readmk[i]='I';
	r_readmk[RREAD_LEN]='\0';

	pr_no=0; err_num=0;err_preads=0;
	l_len=qseq[0][0].l;
	r_len=qseq[1][0].l;
	a_len= l_len>r_len? l_len : r_len;

	PAIR_READS_NUM = COVER_RATE * a_len /(LREAD_LEN+RREAD_LEN);
	printf("Enter gen_pair_reads():\n\tPAIR_READS_NUM==[%d]   total cover rate=[%f]\n", PAIR_READS_NUM, COVER_RATE);
	fprintf(outresult, "\ngen_pair_reads():\nPAIR_READS_NUM==[%d]   total cover rate=[%f]\n",
				PAIR_READS_NUM, COVER_RATE);

	//PAIR_READS_NUM = PAIR_READS_NUM - pr_no;
	while (pr_no<PAIR_READS_NUM){ 
		m = zrand(MAXSS)<0.5? 0:1;
		k = zrand(MAXSS)<0.5? 0:1;
		beg_end(0, l_beg, r_end, qseq[m][k].l);
		//generate one pair_ends reads.
		n = out_read(l_beg, r_end, l_readmk, r_readmk, pr_no, qseq, m, k, ref_n);
		err_num += n;
		if (n>0) err_preads++;
		pr_no++; // count pair_end reads order
		if (pr_no % 1000000 == 0) printf("\tpr_no=%d\n", pr_no);
	}  //while (pr_no<PAIR_READS_NUM){
	printf("\tgenerated paired_end read num = %d\n", pr_no);
	m = pr_no*(LREAD_LEN+RREAD_LEN);
	printf("\tpair_end reads total = [%d], Error pair reads =[%d] error rate=[%f] in ref1 + ref2.\n",
						pr_no, err_preads, 1.0*err_preads/pr_no);
	printf("\ttotal gchar = [%d] error gchar num=[%d]  total error rate=[%f].\n", m, err_num, 1.0*err_num/m);
	fprintf(outresult, "pair_end reads total = [%d], Error pair reads =[%d] error rate=[%f] in ref1 + ref2.\n",
						pr_no, err_preads, 1.0*err_preads/pr_no);
	fprintf(outresult, "total gchar = [%d] error gchar num=[%d]  total error rate=[%f].\n", 
						m, err_num, 1.0*err_num/m);

	return 0;

}  

/**********************************************/
static int out_single_read(int beg, char *readmk, int pr_no, seq_t qseq[2][2], int n1, int n2, char * ref_name)
{
	char left[LREAD_LEN+1], ch;
	int j, err_num, k, n, total_err;
	float f1; 

	err_num=0; total_err=0;
	for (j=0;j<LREAD_LEN;j++){
		if (zrand(MAXSS)<ERROR_RATE){ err_num++; total_err++; }
		if (n2==0) left[j] = qseq[n1][0].s[beg+j]; 
		else left[j] = qseq[n1][1].s[beg+LREAD_LEN-1-j];
	}
	left[LREAD_LEN] = '\0';
	if (err_num>0){
		for (j=0; j<err_num; j++){
			k = (int)(zt_stdfun()*LREAD_LEN);
			if (k<0) k=-k;
			k = k % LREAD_LEN;
			n = LREAD_LEN-1-k;
			k = gch_to_num(left[n]);
			k = (k + (int)(zrand(MAXSS) * 3.0 + 1)) & 3;
			left[n] = num_to_gch(k);
		}
	}

	fprintf(outleft, "@%s_%d_%d_%d_0:0:0_0:0:0_%x\n", ref_name, beg+1, n1, n2, pr_no);
	fprintf(outleft, "%s\n", left);
	fprintf(outleft, "+\n%s\n", readmk);

//	fprintf(outleft, "%s\n", left);
	return total_err;
}

/**********************************************/
int gen_single_read(seq_t qseq[2][2], char *ref_n)
{
	char readmk[LREAD_LEN+1];//={'I'};
	seq_t *s[2];
	int PAIR_READS_NUM;
	int beg, l_len, r_len, a_len;
	int i, k,m, n, pr_no, err_num, err_preads; 

	for (i=0;i<LREAD_LEN;i++) readmk[i]='I';
	readmk[LREAD_LEN]='\0';
	srand(unsigned (time(0)));  // Init Random data 

	l_len=qseq[0][0].l;
	r_len=qseq[1][0].l;
	a_len= l_len>r_len? l_len : r_len;

	PAIR_READS_NUM = COVER_RATE * a_len /LREAD_LEN;
	fprintf(outresult, "\ngen_single_read():\nSingle read num==[%d]   total cover rate=[%f]\n",
				PAIR_READS_NUM, COVER_RATE);
	pr_no=0; err_num=0;err_preads=0;

	//PAIR_READS_NUM = PAIR_READS_NUM - pr_no;
	while (pr_no<PAIR_READS_NUM){
		m=zrand(MAXSS)<0.5? 0:1;
		k=zrand(MAXSS)<0.5? 0:1;
		do{
			beg = zrand(MAXSS) * qseq[m][k].l;
		}while(beg+LREAD_LEN >= qseq[m][k].l);
		//generate one pair_ends reads.
		n = out_single_read(beg, readmk, pr_no, qseq, m, k, ref_n);
		err_num += n;
		if (n>0) err_preads++;
		pr_no++; // count pair_end reads order
	}  //for (i=0;i<READS_NUM){
	m = pr_no*LREAD_LEN;
	fprintf(outresult, "single reads total = [%d], Error reads =[%d] error rate=[%f] in ref1 + ref2.\n",
						pr_no, err_preads, 1.0*err_preads/pr_no);
	fprintf(outresult, "total gchar = [%d] error gchar num=[%d]  total error rate=[%f].\n", 
						m, err_num, 1.0*err_num/m);
	return 0;

}  

