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

char Ref_Name[50];
char ref_str[50];
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
FILE* outzbeg;

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

		FILE * fp=fopen(filename,"r");
		ret_eof = fscanf(fp, "%s", ch0);
		reads_num=0;
	 	reads_match_num =0;

		while((ret_eof = fscanf(fp, "%s", ch0)) != EOF){
			while ( !(pos=strstr(ch0, ref_str)) && (ret_eof != EOF) ) ret_eof = fscanf(fp, "%s", ch0);
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
			}else fprintf(outzbeg, "%s %s %s %s\n", ch0, ch1, ch2, ch3);

			if (reads_num%10000000==0) {cout<<"readno="<<reads_num<<"  read_match_num="<<reads_match_num<<endl;}

		} 

	fclose(fp);
	
}
 /**********************************************************/

int main(int argc,char *argv[])
{ 
	
	int i,j,k,m;
	char outfn[100];

   	if(argc<5) {printf("Error:Missing parameters!\nUsage: mapfq *.fa left.fq right.fq fqout.txt\n");  return 0; }
	
	FILE *inref = fopen(argv[1], "r");
	FILE* outMatch = fopen(argv[4], "w");
	sprintf(outfn, "zbeg_%s", argv[4] );
	outzbeg  = fopen(outfn, "w");

	fscanf(inref,"%s\n",ch0);
	strcpy(Ref_Name, &ch0[1]);
	sprintf(ref_str, "%s_", Ref_Name);
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
	
	m=0;
	for(k=0; k < REF_LEN; k++) {
        //if(genotype_site[k] ) {   
			fprintf(outMatch,"%9d  %3d %3d %3d %c %c %1d\n",k, diff[k],same[k],cover[k],normal[k],ref[k],genotype_site[k]);
	  		 m++;
       // }
    }
   
   
	delete [] same;
	delete [] diff;
	delete [] cover;
	delete [] diffchar;
	delete [] genotype_site;
	delete [] normal;

	fclose(outMatch);
	fclose(outzbeg);

	return 0 ;
}












