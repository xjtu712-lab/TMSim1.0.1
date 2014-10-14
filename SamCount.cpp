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
char ch9_1[TERM_LEN];
char ch9_2[TERM_LEN];
char *pos;

int reads_num,reads_match_num,ret_eof;

/*******************************************************/
FILE *FileOpen(char *fileName,char *mode){
	FILE *fp=NULL;
	fp=fopen(fileName,mode);
	if(fp==NULL){
		fprintf(stdout,"Error:Cannot Open the file %s for%s !\n",fileName,mode);
		fflush(stdout);
		exit(0);
	}
	return fp;
}


/**************************************************************/
static int match_diff(int lbeg, char *reads){
    int j,len, n;
	len=strlen(reads);
	n =0;
	for(j=0; j < len; j++){
		if (_toupper(ref[lbeg+j]) != _toupper(reads[j]) ){
			n++;
		}
	}
	return n;
} // end of match_reads(int beg, char *lrmark, char *reads)


/**************************************************************/
static int match_reads(int lbeg, char *reads){
    int j,len;
	len=strlen(reads);

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
 void read_sam(char *filename){
		FILE * fp=fopen(filename,"r");
		char errfile[100];
	    int lbeg1,lbeg2, rbeg, gaplen, n, obeg, oend, j, k;
		char ch;
		int m1,m2;
		
		ret_eof = fscanf(fp, "%s", ch0);
		sprintf(errfile, "err_%s", filename);
		FILE * errfp=fopen(errfile,"w");

		reads_num=0;
	 	reads_match_num =0;

		while(ret_eof != EOF){
			ret_eof = fscanf(fp, "%s", ch0);
			while ( !(pos=strstr(ch0, ref_str)) && (ret_eof != EOF) ) ret_eof = fscanf(fp, "%s", ch0);
			if (ret_eof == EOF) break;

			j=5; k=0;
			while (ch0[j] != '_') { ch1[k++] = ch0[j++]; };
			ch1[k] = '\0';
			obeg=atoi(ch1);
	
			k=0; j++;
			while (ch0[j] != '_') { ch2[k++] = ch0[j++]; }
		    ch2[k] = '\0';
		    oend=atoi(ch2);

			ret_eof=fscanf(fp, "%s %s %s %s %s %s %s %s %s ", ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9_1);
			if (ret_eof==EOF) break;

			lbeg1= atoi(ch3)-1;
		 	gaplen = atoi(ch8);
            if(gaplen<0)gaplen=-gaplen;
		 	reads_num++;
			
			m1=0;
			if (gaplen != 0 && strcmp(ch5,"100M")==0 ){ 
                m1=1;
             
			}else {
				n = match_diff(lbeg1, ch9_1);
				if (lbeg1+1 != obeg && lbeg1+1 != oend-99 ){
					ch = ref[lbeg1+100];
					ref[lbeg1+100] = '\0';
					fprintf(errfp, "%s %s %s %s %s %s %s %s %s\n\t\t%s n=%d\n\t\t%s\n", ch0, ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9_1, n, &ref[lbeg1]);
					ref[lbeg1+100] = ch;
				}
			}

			ret_eof = fscanf(fp, "%s ", ch0);
			while ( !(pos=strstr(ch0, ref_str)) && (ret_eof != EOF) ) ret_eof = fscanf(fp, "%s", ch0);
			if (ret_eof == EOF) break;

			j=5; k=0;
			while (ch0[j] != '_') { ch1[k++] = ch0[j++]; };
			ch1[k] = '\0';
			obeg=atoi(ch1);
	
			k=0; j++;
			while (ch0[j] != '_') { ch2[k++] = ch0[j++]; }
		    ch2[k] = '\0';
		    oend=atoi(ch2);

			ret_eof=fscanf(fp, "%s %s %s %s %s %s %s %s %s ", ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9_2);
			if (ret_eof==EOF) break;
			lbeg2 = atoi(ch3)-1;
			gaplen = atoi(ch8);
			reads_num++;

			m2=0;
			if (gaplen != 0 && strcmp(ch5,"100M")==0 ){ 
				m2=1;					

			}else{
				n = match_diff(lbeg2, ch9_2);
				if ( lbeg2+1 != obeg && lbeg2+1 != oend-99 ){
					ch = ref[lbeg2+100];
					ref[lbeg2+100] = '\0';
					fprintf(errfp, "%s %s %s %s %s %s %s %s %s\n\t\t%s n=%d\n\t\t%s\n", ch0, ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9_2, n, &ref[lbeg2]);
					ref[lbeg2+100] = ch;
				}
			}
			if (reads_num%10000000==0) {cout<<"readno="<<reads_num<<"  reads_match_num= "<<reads_match_num<<endl;}
			if(m1 && m2) {
				match_reads(lbeg1, ch9_1);
				match_reads(lbeg2, ch9_2);
				reads_match_num += 2;
			}

		} 

	fclose(fp);
	fclose(errfp);
	
}
 /**********************************************************/

int main(int argc,char *argv[])
{ 
	int i,j,k,m;
   	if(argc<4) {printf("Error:Missing parameters!\nUsage: samcount  ref.fa  *.sam  sammatchout.txt\n\n"); return 0; }
	
	FILE *inref = fopen(argv[1], "r");
	FILE *insam = fopen(argv[2], "r");
	FILE* outMatch=fopen(argv[3], "w");

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

	for(k=0; k<REF_LEN; k++){ref[k] = _toupper(ref[k]); same[k]=0; diff[k]=0; cover[k]=0;}

	read_sam(argv[2]);

	for(k=0; k < REF_LEN; k++){
		 if (diff[k] == 0) {normal[k] = ref[k]; genotype_site[k]=0; } // AA TYPE
		 else if (diff[k] == cover[k]) {
					normal[k] = diffchar[k];
					genotype_site[k]=2; // BB TYPE
			   }else if (diff[k] > same[k]) {normal[k] = diffchar[k];genotype_site[k]=1;} // AB TYPE
						else {normal[k] = ref[k]; genotype_site[k]=1; } // AB TYPE
 
	 }  

	cout<<"reads_num="<<reads_num<<"reads_match_num="<<reads_match_num<<"NUM-MATCH="<<reads_num-reads_match_num<<endl;

	m=0;
	for(k=0; k < REF_LEN; k++) {

		if(genotype_site[k]>0 && cover[k]>=3 && diff[k]>1){
			fprintf(outMatch,"%9d  %3d %3d %3d %c %c %1d\n",k, diff[k],same[k],cover[k],normal[k],ref[k],genotype_site[k]);
        	m++;
        }
    }
   
   
	delete [] same;
	delete [] diff;
	delete [] cover;
	delete [] diffchar;
	delete [] genotype_site;
	delete [] normal;

	fclose(outMatch);

	return 0 ;
}












