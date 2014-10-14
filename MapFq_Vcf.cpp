#include <stdio.h>
#include <stdlib.h>
#include <iostream> 
#include <string> 
#include <fstream>
#include <string.h>
using namespace std;

const int REF_MAX = 250000000;
const int TERM_LEN = 200;  
int REF_LEN=REF_MAX;

char ref[REF_MAX];
unsigned char *same;
unsigned char *diff;
unsigned char *cover;
unsigned char *genotype_site;
char  *diffchar;
char  *normal;

char ch0[TERM_LEN];
char ch1[TERM_LEN];
char ch2[TERM_LEN];
char ch3[TERM_LEN];
char ch4[TERM_LEN];
char ch5[TERM_LEN];
char ch6[TERM_LEN];
char ch7[TERM_LEN];
char ch8[TERM_LEN];
char ch9[TERM_LEN];
char *pos;


int reads_num =0;
int reads_match_num=0;
int ret_eof;
char fhead[50];	


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


/**************************************************************/
static int match_reads(int lbeg, int zfmk, char *reads){
    int j,len;
	char ch[200];

	len=strlen(reads);
	if (zfmk){
		for (j=0; j<len; j++) ch[j] = gch_rev( reads[len-1-j] );
		reads = ch;
	}	

	for(j=0; j < len; j++){
		cover[lbeg+j]++;
		if (_toupper(ref[lbeg+j]) != _toupper(reads[j]) ){
			diff[lbeg+j]++;
			diffchar[lbeg+j] = reads[j];
		}else same[lbeg+j]++;
	}
	return 0;
} // end of match_reads(int beg, char *lrmark, char *reads)


/*******************************************************/
 void read_sam(char *filename, int lrmk){
	    int lbeg, rbeg, gaplen, rend, j, k, zbeg, zfmk, mark;
		char outhead[50];

		FILE * fp=fopen(filename,"r");
		ret_eof = fscanf(fp, "%s", ch0);
		reads_num=0;
	 	reads_match_num =0;
		sprintf(outhead,"@%s_",fhead);
		
		while((ret_eof = fscanf(fp, "%s", ch0)) != EOF){

			while ( !(pos=strstr(ch0, outhead)) && (ret_eof != EOF) ) ret_eof = fscanf(fp, "%s", ch0);
			//while ( (strcmp(ch0, outhead)!=0)&& (ret_eof != EOF) ) ret_eof = fscanf(fp, "%s", ch0);
			if (ret_eof == EOF) break;

			j=6; k=0;
			while (ch0[j] != '_') { ch1[k++] = ch0[j++]; };
			ch1[k] = '\0';
			lbeg=atoi(ch1);
	
			k=0; j++;
			while (ch0[j] != '_') { ch2[k++] = ch0[j++]; }
		    ch2[k] = '\0';
		    rend=atoi(ch2);

			k=0; j++;
			while (ch0[j] != '_') { ch3[k++] = ch0[j++]; }
		    ch3[k] = '\0';
		    zbeg=atoi(ch3);

			k=0; j++;
			while (ch0[j] != '_') { ch4[k++] = ch0[j++]; }
		    ch4[k] = '\0';
		    zfmk = atoi(ch4);


			ret_eof=fscanf(fp, "%s %s %s", ch1, ch2, ch3);
		 	if (ret_eof==EOF) break;
		 	reads_num++;

			if (lrmk == 0) { if (zfmk == 0)  mark = 0; else mark = 1; }
			else { if (zfmk == 0)  mark = 1; else mark = 0; }
			if (zbeg>=0 && zbeg <REF_LEN){
				match_reads(zbeg-1, mark, ch1);
				reads_match_num++;
			}

			if (reads_num%10000000==0) {cout<<"readno="<<reads_num<<"  read_match_num="<<reads_match_num<<endl;}

		} 

	fclose(fp);
	
}
 /**********************************************************/

int main(int argc,char *argv[])
{ 
	
	int i,j,k,m;
	char outfn[100];
	char refname[50];
	
	int ac=0;
	char gt[10];
   	if(argc<5) {printf("Error:Missing parameters!\nUsage: mapfq  ref.fa  left.fq  right.fq  fq.vcf\n");  return 0; }
	
	FILE *inref = fopen(argv[1], "r");

	FILE* outvcf=fopen(argv[4],"w");

	strcpy(refname,argv[1]);
	for(i=0;i<strlen(refname);i++) if(refname[i]=='.')break;
	if(i<strlen(refname)){strncpy(fhead,&refname[0],i);fhead[i]='\0';}
	else strcpy(fhead,refname);
	//cout<<"fhead1="<<fhead<<endl;
    

	fscanf(inref,"%s\n",ch0);
	i=0;
	while(fscanf(inref,"%s\n",&ref[i])!=EOF){i=i+50;}
	fclose(inref);
	REF_LEN=strlen(ref);
	cout<<"REF_LEN="<<REF_LEN<<endl;

	same  = new unsigned char[REF_LEN];
	diff  = new unsigned char[REF_LEN];
	cover = new unsigned char[REF_LEN];
	diffchar = new char[REF_LEN];
	normal = new char[REF_LEN];
	genotype_site=new unsigned char[REF_LEN];

	for(k=0; k<REF_LEN; k++){same[k]=0; diff[k]=0; cover[k]=0;}

	read_sam(argv[2], 0); // 0---left.fq   1---right.fq

	read_sam(argv[3], 1);

	for(k=0; k < REF_LEN; k++){
		 if (diff[k] == 0) {normal[k] = ref[k]; genotype_site[k]=0; } // AA TYPE
		 else if (diff[k] == cover[k]) {
					normal[k] = diffchar[k];
					genotype_site[k]=2; // BB TYPE
			   }else if (diff[k] > same[k]) {normal[k] = diffchar[k];genotype_site[k]=1;} // AB TYPE
						else {normal[k] = ref[k]; genotype_site[k]=1; } // AB TYPE
 
	 }  

	cout<<"reads_num="<<reads_num<<"    reads_match_num="<<reads_match_num<<"    NUM-MATCH="<<reads_num-reads_match_num<<endl;
	
	fprintf(outvcf,"##fileformat=VCFv4.1\n");
	fprintf(outvcf,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n");
	fprintf(outvcf,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n");
	fprintf(outvcf,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n");
	fprintf(outvcf,"##INFO=<ID=RU,Number=1,Type=String,Description=\"Tandem repeat unit (bases)\">\n");
	fprintf(outvcf,"##INFO=<ID=STR,Number=1,Type=Flag,Description=\"Variant is a short tandem repeat\">\n");
	fprintf(outvcf,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n");
	fprintf(outvcf,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(outvcf,"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n");
	fprintf(outvcf,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n");
	fprintf(outvcf,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n");

//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO							FORMAT	NA12878
//  9	    ***	.	*	*	.	  	.		AC=*;AF=*;AN=2;DP=1;RU=*;STR	GT:AD:DP	1/1:0,1:19	***	.	*	*	.	.	AC=*;AF=*;AN=2;DP=1;RU=*;STR	GT:AD:DP	1/1:0,1:1
	

	m=0;
	for(k=0; k < REF_LEN; k++) {
        if(genotype_site[k] ) {   
					
			//if(genotype_site[k] ) {  
				if(genotype_site[k]==1){ac=1;strcpy(gt,"0/1");}else {ac=2;strcpy(gt,"1/1");}	
				fprintf(outvcf,"%s\t%d\t.\t%c\t%c\t.\t.\tAC=%d;AF=%5.2f;AN=2;DP=%d\tGT:AD:DP\t%s:%d,%d:%d\n",fhead,k,ref[k],diffchar[k],
						ac,1.0*diff[k]/cover[k],cover[k],gt,same[k],diff[k],cover[k]);
	  		 
        }
		m++;
    }
   
   
	delete [] same;
	delete [] diff;
	delete [] cover;
	delete [] diffchar;
	delete [] genotype_site;
	delete [] normal;

	fclose(outvcf);

	return 0 ;
}












