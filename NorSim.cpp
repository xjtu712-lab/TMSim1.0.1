
/*norsim---the program of normal simulator for data generator file(*.sim, *_AB_chged.idx)
 *At present ,most of studies adopt the joint analysis of tumor-normal pair data in cancer.    
 *According to this characteristic ,this program realizes to set variation of sites 
 *that exits in human individual genomic sequences of normal cells .*/


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
float	GV_MUT_RATE=0.001;
float	ERR_RATE=0.00;
float 	INDEL_FRAC=0.0;
float	INDEL_EXTEND = 0.3;
float 	DEL_RATE= 0.5;
float	BB_RATE = 0.33333;

int 	LONG_INDEL = 0;
int 	OUT_AB = 0;

FILE* ingenfp=0;
FILE* indelfp=0;
FILE* outindelfp=0;
FILE* outsimfp=0;
FILE* outabfp=0;
FILE* outresult=0;

#define INIT_SEQ(seq) (seq).s = 0; (seq).l = (seq).max = 0
#define CLEAR_SEQ(seq) free( (seq).s )
#define INIT_MUTSEQ(mseq) (mseq).s = 0; (mseq).l = (mseq).max = 0; (mseq).gnum = (mseq).vnum = 0
#define CLEAR_MUTSEQ(mseq) free( (mseq).s )

enum muttype_t {NOCHANGE = 0, INSERT = 0x80, SUBSTITUTE = 0x40, DELETE = 0xc0, GTYPEBB=0x20, GTYPEAB=0x10, LONGINMK = 0xf0000000 };

typedef unsigned int mut_t;

static mut_t mutmsk = (mut_t)0xc0;
static int SEQ_BLOCK_SIZE = 512;

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
	//cout<<endl<<"ref_name="<<chrname<<endl;
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
	//cout<<endl<<"seq->l= "<<seq->l<<"  max= "<<max<<"  ACGT.num= "<<n<<endl;
	return l;
}

/**********************************************/
int set_pos_sim(mutseq_t *hap1, mutseq_t *hap2, int k){ // return chged num
	int deleting = 0;
	mutseq_t *ret[2];
	int tmp, num, m, max, chgnum;
	double r;
	int c0, c1, c;

	srand48(time(0)); 
	ret[0] = hap1; ret[1] = hap2; m=0;
	max = hap1->l;
	c0 = ret[0]->s[k];
	c1 = ret[1]->s[k];
	chgnum = 1;
	if (drand48() >= INDEL_FRAC) { // substitution
		r = drand48();
		c = (c0 + (int)(r * 3.0 + 1)) & 3;
		if (drand48() < BB_RATE) { //hom
			ret[0]->s[k] = SUBSTITUTE|GTYPEBB|c;
			ret[1]->s[k] = SUBSTITUTE|GTYPEBB|c;
		} else { // het
			tmp = drand48()<0.5?0:1;
			ret[tmp]->s[k] = SUBSTITUTE|GTYPEAB|c;
		}
	} else { // indel
		if (drand48() < DEL_RATE) { // deletion
			if (drand48() < BB_RATE ) { // hom-del
				ret[0]->s[k] |= DELETE|GTYPEBB;
				ret[1]->s[k] |= DELETE|GTYPEBB;
	            deleting = 3;
			} else { // het-del
				deleting = drand48()<0.5? 1:2;
				ret[deleting-1]->s[k] |= DELETE|GTYPEAB;
			}
			k++;
			c0 = ret[0]->s[k];
			c1 = ret[1]->s[k];
			while(drand48()<INDEL_EXTEND  && k < max && c0 < 4 && c1 < 4 ) {
				if (deleting<3) ret[deleting-1]->s[k] |= DELETE|GTYPEAB; 
        		else ret[1]->s[k] = ret[0]->s[k] |= DELETE|GTYPEBB; 
				chgnum++;
				k++; 
				c0 = ret[0]->s[k];
				c1 = ret[1]->s[k];
			}
		} else { // insertion
			int num_ins = 0, ins = 0;
			do {
				num_ins++;
				ins = (ins << 2) | (int)(drand48()*4.0);
			} while (num_ins < 10 && (drand48() < INDEL_EXTEND) );
			chgnum = num_ins; //insert number ;
			if (drand48() < BB_RATE ) { // hom-ins
				ret[0]->s[k] = ret[1]->s[k] = (num_ins << 28) | (ins << 8) | INSERT | GTYPEBB | c0;
			} else { // het-ins
				tmp = drand48()<0.5?0:1;
				ret[tmp]->s[k] = (num_ins << 28) | (ins << 8) | INSERT | GTYPEAB | c0;
			}
		}
	} //end indel // mutation
	return chgnum;
}

/*******************************************************/
void set_mut_sim(mutseq_t *hap1, mutseq_t *hap2)
{
	mutseq_t *ret[2];
	int n, max, num, k;
	double r;
	mut_t c0, c1, c;

	srand48(time(0)); 
	ret[0] = hap1; ret[1] = hap2;
	max = ret[0]->l;
	num = ret[0]->gnum * GV_MUT_RATE - ret[0]->vnum;
	n = 0;
	while(n < num){ //set mutation and not normal AB mode
		k = drand48()*max;
		c0 = ret[0]->s[k];
		c1 = ret[1]->s[k];
		while ( c0 >= 4 || c1 >= 4){
			k = drand48()*max;
			c0 = ret[0]->s[k];
			c1 = ret[1]->s[k];
		}
		n += set_pos_sim( hap1, hap2, k);
	}  //end of while(n < num)
	if (num > 0) { ret[0]->vnum += num; ret[1]->vnum += num;}
	printf("\tset_mut_sim(): max_len=%d REAL_len=%d mut_num = %d total mut=%d\n",
				 max, ret[0]->gnum, num, ret[0]->vnum);
	fprintf(outresult,"set_mut_sim(): max_len=%d REAL_len=%d mut_num = %d total mut=%d\n",
				 max, ret[0]->gnum, num, ret[0]->vnum);
	return;
}

/*******************************************************/
int addinsert_iterm(indel_t *indelp, int idno, int idbeg, int idlen, int gtype, int idnum, int idvar, char *idstr, int *n2){
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
	*n2 = *n2 + idlen;
	return 0;
}


/*******************************************************/
void init_indel_sim( seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, indel_t *indelp)
{
	int fret, no, idno, idbeg, idlen, gtype, idnum, idvar;
	char ch=0, name[50], *idstr, *idmark; //[100000];
	int i,j,k, n1, n2, h, n, del_no, in_no;
	mutseq_t *ret[2];

	ret[0] = hap1; ret[1] = hap2;
	ret[0]->l = seq->l; ret[1]->l = seq->l;
	ret[0]->max = seq->max; 
	ret[1]->max = seq->max;
	ret[0]->s = (mut_t *)calloc(seq->max, sizeof(mut_t));
	ret[1]->s = (mut_t *)calloc(seq->max, sizeof(mut_t));

	n=0;
	for (i = 0; i < seq->l; ++i) {
		ret[0]->s[i] = ret[1]->s[i] = (mut_t)gch_to_num(seq->s[i]);
		if ( (ret[0]->s[i]&0xf) < 4) n++;
	}
	ret[0]->gnum = n;  ret[1]->gnum = n;
	ret[0]->vnum = 0;  ret[1]->vnum = 0;
	
	// handle for long delete area from long indel file
	if (!indelfp || !LONG_INDEL) return;
	idstr = new char [100000];
	idmark = new char [100000];
	n1=0; n2 = 0; del_no = 0; in_no = 1;
	while ( fscanf(indelfp, "%d %d %d %d %d %d %s\n", &idno, &idbeg, &idlen, &gtype, &idnum, &idvar, idstr) != EOF ){
		if (idlen<0) {  // long delete from long indel file!!!
			idlen=-idlen;
			if (idbeg<0 || idbeg+idlen>=REF_LEN) gtype=-1;
			switch(gtype){
			case 0: break;
			case 1: if (zrand(MAXSS)<0.5) for (int i=0; i<idlen; i++) ret[0]->s[idbeg+i] |= DELETE|GTYPEAB;
				else for (int i=0; i<idlen; i++) ret[1]->s[idbeg+i] |= DELETE|GTYPEAB;
				n1 = n1 + idlen;
				del_no ++;
				break;
			case 2:
				for (int i=0; i<idlen; i++){ ret[0]->s[idbeg+i] = ret[1]->s[idbeg+i] |= DELETE|GTYPEBB;}
				n1 = n1 + idlen;
				del_no ++;
				break;
			default: break;
			}
			if (gtype == -1) printf("\tInvalid DEL  No= %d   beg: %d    len: %d\n", del_no, idbeg, idlen);
		}else {  //long Insert from long indel file
			if (idbeg<0 || idbeg+idlen>=REF_LEN || idnum<1) gtype=-1;
			switch(gtype){
			case 0: break;
			case 1: h = (int)(zrand(MAXSS)*2);
				ret[h]->s[idbeg] = (ret[h]->s[idbeg] & 0x03 );
				ret[h]->s[idbeg] |= (LONGINMK|INSERT|GTYPEAB|(in_no<<8));
				break;
			case 2: 
				ret[0]->s[idbeg] = (ret[0]->s[idbeg] & 0x03 );
				ret[1]->s[idbeg] = (ret[1]->s[idbeg] & 0x03 );
				ret[0]->s[idbeg] |= (LONGINMK|INSERT|GTYPEBB|(in_no<<8));
				ret[1]->s[idbeg] |= (LONGINMK|INSERT|GTYPEBB|(in_no<<8));
				break;
			default: break;
			}
			if (gtype > 0) {
				int res_num, curr_num;
				if (idvar == 0){  //repeat without var genos
					for(i=0; i<idnum; i++) addinsert_iterm(indelp, idno, idbeg, idlen, gtype, idnum, idvar, idstr, &n2);
				}else { //repeat with var genos
					for(i=0; i<strlen(idstr); i++) idmark[i] = 0;
					res_num = idnum;
					while(res_num>4){
						k = zrand(MAXSS)*(0.5*res_num-1)+1;
						for(i=0; i<k; i++) addinsert_iterm(indelp, idno, idbeg, idlen, gtype, idnum, idvar, idstr, &n2);
						res_num = res_num - k;
						j = zrand(MAXSS)*strlen(idstr);
						while (idmark[j]) j = zrand(MAXSS)*strlen(idstr);
						idmark[j] = 1;
						i = (gch_to_num(idstr[j]) + (int)(zrand(MAXSS)* 3.0 + 1)) & 3;
						idstr[j] = num_to_gch(i);
					}
					if (res_num>1){
						k = zrand(MAXSS)*(res_num-1)+1;
						for(i=0; i<k; i++) addinsert_iterm(indelp, idno, idbeg, idlen, gtype, idnum, idvar, idstr, &n2);
						res_num = res_num - k;
						j = zrand(MAXSS)*strlen(idstr);
						while (idmark[j]) j = zrand(MAXSS)*strlen(idstr);
						idmark[j] = 1;
						i = (gch_to_num(idstr[j]) + (int)(zrand(MAXSS)* 3.0 + 1)) & 3;
						idstr[j] = num_to_gch(i);
						for(i=0; i<res_num; i++) addinsert_iterm(indelp, idno, idbeg, idlen, gtype, idnum, idvar, idstr, &n2);
					}else if (res_num == 1) addinsert_iterm(indelp, idno, idbeg, idlen, gtype, idnum, idvar, idstr, &n2);
				}//if(ifvar ==0 )...
				in_no++;
			}else printf("\tInvalid INSERT  No= %d   beg: %d     len: %d\n", in_no, idbeg, idlen);
		}//if...else...
	} //while (reading file .....)
	printf("\tEND:: Long del_no=%d  del_len=%d;  Long insert num=%d  in_lens=%d;   Long insert table size=%d\n", del_no, n1, in_no-1, n2, indelp->l);
	delete [] idstr;
	delete [] idmark;

	return;
}

/**********************************************/
void print_simfile(mutseq_t *hap1, mutseq_t *hap2, char *refname, indel_t *indelp)
{
	int i, n, n0,n1,cnt0,cnt1, ab_no;

	n=n0=n1=cnt0=cnt1=0; ab_no=0;
	fprintf(outsimfp,"&%s\n", refname);
	for (i = 0; i < hap1->l; ++i) {
		if ((hap1->s[i]&mutmsk)!=DELETE && (hap1->s[i] & 0xf)<4 ){
			n0++;
			if ((hap1->s[i] & mutmsk)==INSERT ) n0=n0+(hap1->s[i]>>28);
		}
		if ((hap2->s[i] & mutmsk)!=DELETE && (hap2->s[i] & 0xf)<4 ){
			n1++;
			if ((hap2->s[i] & mutmsk)==INSERT ) n1=n1+(hap2->s[i]>>28);
		}

		if ((hap1->s[i] & mutmsk) || (hap2->s[i] & mutmsk) ) {
			n++;
			if( OUT_AB && ( (hap1->s[i] & 0xf0)==0x50 || (hap2->s[i] & 0xf0)== 0x50 ) ){
				fprintf(outabfp, "%10u\t%10u\t%10u\n", i, hap1->s[i], hap2->s[i]);
				ab_no++;
			}
			if (hap1->s[i] & mutmsk) cnt0++;
			if (hap2->s[i] & mutmsk) cnt1++;
			fprintf(outsimfp, "%10u\t%10u\t%10u\n", i, hap1->s[i], hap2->s[i]);
		}
	}
	// out long insert TABLE from long indel file
	if (indelfp && LONG_INDEL && indelp->l){
		for (i = 0; i < indelp->l; ++i) {
			fprintf(outindelfp, "%d %d %d %d %d %d \n%s\n", indelp->s[i].indel_no, indelp->s[i].indel_beg, 
					indelp->s[i].indel_len, indelp->s[i].indel_type, indelp->s[i].indel_num, indelp->s[i].indel_var, indelp->s[i].str);
		}
	}
	
	printf("\t ref1_len=%d  ref1 var num = %d Var rate = %.10f\n\n",n0,cnt0,1.0*cnt0/n0);
	printf("\t ref2_len=%d  ref2 var num = %d Var rate = %.10f\n\n",n1,cnt1,1.0*cnt1/n1);
	printf("\t Ref_len =%d  all  var num = %d Var rate = %.10f AB_mode num= %d \n\n", REF_LEN,n,1.0*n/REF_LEN,ab_no);
	fprintf(outresult,"\n ref1_len=%d  ref1 var num = %d Var rate = %.10f\n\n",n0,cnt0,1.0*cnt0/n0);
	fprintf(outresult,"\n ref2_len=%d  ref2 var num = %d Var rate = %.10f\n\n",n1,cnt1,1.0*cnt1/n1);
	fprintf(outresult,"\n Ref_len =%d  all  var num = %d Var rate = %.10f AB_mode num= %d \n\n", 
			REF_LEN,n,1.0*n/REF_LEN,ab_no);
	return;
}

/**********************************************/
#define PACKAGE_VERSION "1.0.1"
static int simu_usage()
{
	printf("\n\nProgram: NorSim (the Normal Simulator for the data generator)\n");
	printf("Version: %s \n", PACKAGE_VERSION);
	printf("Contact: Yu Geng <gengyu@stu.xjtu.edu.cn>;Zhongmeng Zhao <zmzhao@mail.xjtu.edu.cn>\n");
	printf("Usage:   norsim [options] <ref.fa> <nor.sim>\n");
	printf("Options: -r FLOAT	mutation rate of GV [%12.10f]\n", GV_MUT_RATE); //[0.001]
	printf("         -R FLOAT      	fraction of indels [%10.6f]\n",INDEL_FRAC); //0.1
	printf("         -X FLOAT      	probability an indel is extended [%10.6f]\n",INDEL_EXTEND);//0.3
	printf("         -D FLOAT      	delete rate in indel [%10.6f]\n", DEL_RATE);//0.5
	printf("         -B FLOAT      	BB rate in mutation [%10.6f]\n", BB_RATE);//0.33333
	printf("         -I <delpos.txt>  	input long indel set file\n");
	printf("         -A <nor_AB.idx>  	output the positions of AB mutation type in normal\n");
	printf("         -o <nor_simout.txt>  	result for runing case\n");
	return 1;
}

/**********************************************/
/**********************************************/
int main(int argc, char *argv[])
{	char c, cc, fhead[50];
	int i;
	char ingenf[50];
	char indelf[50];
	char outindelf[50];
	char outsimf[50];
	char outabf[50];
	char outresultf[50];
	seq_t	 seq;
	mutseq_t sim[2];
	char ref_name[50];
	indel_t  indelp;

	srand48(time(0)); 
	outresultf[0]='\0';
	indelf[0]='\0';
	cc=0;
	while ((c = getopt(argc, argv, "r:R:X:D:B:I:A:o:")) >= 0) {
		switch (c) {
		case 'r': GV_MUT_RATE     = atof(optarg); if (GV_MUT_RATE<0 || GV_MUT_RATE>1) cc=1; break;
		case 'R': INDEL_FRAC   = atof(optarg); if (INDEL_FRAC<0 || INDEL_FRAC>1) cc=1; break;
		case 'X': INDEL_EXTEND = atof(optarg); if (INDEL_EXTEND<0 || INDEL_EXTEND>1) cc =1; break;
		case 'D': DEL_RATE = atof(optarg);     if (DEL_RATE<0 || DEL_RATE>1) cc =1; break;
		case 'B': BB_RATE = atof(optarg);      if (BB_RATE<0 || BB_RATE>1) cc =1; break;
		case 'I': LONG_INDEL = 1; strcpy(indelf, optarg);break;
		case 'A': OUT_AB = 1; 	strcpy(outabf, optarg);break;
		case 'o': strcpy(outresultf, optarg);break;
		default : break; 
		}
		if (cc) {optind = argc; break;}
	}
	if (optind != (argc-2) ) return simu_usage();
	strcpy(ingenf, argv[optind]);
	strcpy(outsimf, argv[optind+1]);
	for (i=0; i<strlen(outsimf); i++) if (outsimf[i] == '.') break; 
	if (i<strlen(outsimf) ) { strncpy(fhead, &outsimf[0], i); fhead[i]='\0'; }
	else strcpy(fhead, outsimf); 
	if (strlen(outresultf)==0) sprintf(outresultf, "%s_simout.txt", fhead); 
	if (LONG_INDEL){
		indelfp = FileOpen( indelf,  "r");
		sprintf(outindelf, "%s.idx", outsimf);
		outindelfp = FileOpen( outindelf,  "w"); 
	} else indelfp = NULL;
	if (OUT_AB )	outabfp = FileOpen( outabf,  "w"); else outabfp = NULL;
	ingenfp  = FileOpen( ingenf,  "r"); 
	outsimfp = FileOpen( outsimf,   "w");
	outresult = FileOpen( outresultf, "w");

	INIT_SEQ(seq);
	INIT_MUTSEQ(sim[0]);
	INIT_MUTSEQ(sim[1]);
	INIT_SEQ(indelp);

	//* read reference fasta */	
	REF_LEN = seq_read_fasta(ingenfp, &seq, ref_name, 0);
	printf("seq_read_fasta() return! REF_LEN = %d\n\n", REF_LEN);

	//* init sim and set long DEL */ 
	init_indel_sim(&seq, &sim[0], &sim[1], &indelp);
	printf("init_del_sim return!\n\n");

	//* set mutation for simsulator */
	set_mut_sim( &sim[0], &sim[1] );
	printf("set_mut_sim() return!\n\n");

	//out put simulator file
	print_simfile( &sim[0], &sim[1], ref_name, &indelp);
	printf("print_simfile() return!\n\n");

	CLEAR_SEQ(seq);
	CLEAR_MUTSEQ(sim[0]);
	CLEAR_MUTSEQ(sim[1]);
	for (i=0; i<indelp.l; i++) if (indelp.s[i].str) free(indelp.s[i].str);
	CLEAR_SEQ(indelp);

	if(indelfp) fclose(indelfp); 
	if(outindelfp) fclose(outindelfp); 
	if(outabfp) fclose(outabfp);
	fclose(ingenfp); 
 	fclose(outsimfp);
	fclose(outresult);
	printf("NorSim meet end!\n\n");

	return 0;
}

