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

int reads_num,reads_match_num,ret_eof;


/**************************************************************/
static int match_reads(int lbeg, char *reads){
    int j,len;
	len=strlen(reads);

	for(j=0; j < len; j++){
			cover[lbeg+j]++;
			if (_toupper(ref[lbeg+j]) != _toupper(reads[j]) ){diff[lbeg+j]++;}
			else same[lbeg+j]++;
	}
	return 0;
} // end of match_reads(int beg, char *lrmark, char *reads)


/*******************************************************/
void read_sam(FILE *fp){
int lbeg, rbeg, gaplen;
		ret_eof = fscanf(fp, "%s", ch0);
		reads_num=0;
	 	reads_match_num =0;
		while(ret_eof != EOF){
			ret_eof = fscanf(fp, "%s", ch0);
			while ( !(pos=strstr(ch0, ref_str)) && (ret_eof != EOF) ) ret_eof = fscanf(fp, "%s", ch0);
			if (ret_eof == EOF) break;
			ret_eof=fscanf(fp, "%s %s %s %s %s %s %s %s %s ", ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9);
		 	
			if (ret_eof==EOF) break;
			lbeg = atoi(ch3)-1;
		 	gaplen = atoi(ch8);
		 	reads_num++;
			if (gaplen != 0 && strcmp(ch5,"100M")==0){ 
				match_reads(lbeg, ch9);
			 	reads_match_num++;
			}
			if (reads_num%1000000==0) {cout<<"readno="<<reads_num<<"  ch3= "<<ch3<<endl;}

			ret_eof = fscanf(fp, "%s ", ch0);
			while ( !(pos=strstr(ch0, ref_str)) && (ret_eof != EOF) ) ret_eof = fscanf(fp, "%s", ch0);
			if (ret_eof == EOF) break;

			ret_eof=fscanf(fp, "%s %s %s %s %s %s %s %s %s ", ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9);
			if (ret_eof==EOF) break;
			lbeg = atoi(ch3)-1;
			gaplen = atoi(ch8);
			reads_num++;

			if ( gaplen !=0 && strcmp(ch5,"100M")==0){ 
				match_reads(lbeg, ch9);
			 	reads_match_num++;
			}
			if (reads_num%1000000==0) {cout<<"readno="<<reads_num<<"  ch3= "<<ch3<<endl;}

		} 

	fclose(fp);

}
 /**********************************************************/

int main(int argc,char *argv[])
{ 
	int * n;
	int len, gtype;
	int i,j,k,c,m;
	int beg,end,dis;
	int mb,me;

	if (argc < 6)  {printf("Error:Missing parameters!\nusage: cg  refrence.fa  normal.sam  tumor.sam  normal_geno.txt  tumor_geno.txt \n");return 0;}

	FILE *inref = fopen(argv[1], "r");
	FILE *innorsam = fopen(argv[2],"r");
	FILE *intursam = fopen(argv[3],"r");
	FILE *outNgentype = fopen(argv[4],"w");
	FILE *outTgentype = fopen(argv[5],"w");
	
	fscanf(inref,"%s\n",ch0);
	strcpy(Ref_Name, &ch0[1]);
	sprintf(ref_str, "%s_", Ref_Name);
	i=0;
	while(fscanf(inref,"%s\n",&ref[i])!=EOF){i=i+50;}
	fclose(inref);
	REF_LEN=strlen(ref);
	cout<<REF_LEN<<endl;

	same  = new unsigned char[REF_LEN];
	diff  = new unsigned char[REF_LEN];
	cover = new unsigned char[REF_LEN];

	for(k=0; k<REF_LEN; k++){same[k]=0; diff[k]=0; cover[k]=0;}
	read_sam(innorsam);
	for(k=0; k < REF_LEN; k++){
		if (diff[k] == 0) {gtype=0; } // AA TYPE
			 else if (diff[k] == cover[k]) gtype=2; // BB TYPE
				  else  gtype=1; // AB TYPE
	 	if (gtype) fprintf(outNgentype,"%10d %2d\n", k, gtype);
	 }  
	fclose(outNgentype);

	for(k=0; k<REF_LEN; k++){same[k]=0; diff[k]=0; cover[k]=0;}
	read_sam(intursam);
    for(k=0; k < REF_LEN; k++){
		if (diff[k] == 0) {gtype=0; } // AA TYPE
			 else if (diff[k] == cover[k]) gtype=2; // BB TYPE
				  else  gtype=1; // AB TYPE
	 	if (gtype) fprintf(outTgentype,"%10d %2d\n", k, gtype);
		//if (gtype) cout<<"g["<<k<<"]="<<gtype<<endl;
	 }  
	fclose(outTgentype);

	cout<<"genotype complete!"<<endl;
   
	delete [] same;
	delete [] diff;
	delete [] cover;

	return 0 ;
}












