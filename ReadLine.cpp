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

const int TERM_LEN = 200;  
const int READ_LEN = 100;  

char temp[2048];


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


/**********************************************/
int MAXSS = pow(2, 31)-1;
static double zrand(int max){
	static double r;
	int n;
	n = rand() % max;
	r = 1.0*n/max;
	return r;
}

int readline(FILE *fp, char *s){
int j;
	j=0;
	while(!feof(fp) && (s[j] = fgetc(fp)) != '\n') j++;
	s[j]='\0';
	if (j==0 && feof(fp)) return -1;
	else return j;    
}

int splitline(char *s, char term[][1024]){
int j, k, n, slen;
	j=0; n = 0;
	slen = strlen(s);
	while (j < slen){
		while((s[j]==' ' || s[j]=='\t' || s[j]=='\r' || s[j]=='\n') && j<slen ) j++;
		k = 0;
		while((s[j]!=' ' && s[j]!='\t' && s[j]!='\r' && s[j]!='\n') && j<slen ) { term[n][k++] = s[j++]; }
		term[n][k]='\0';
		n++;
	}
	return n;    
}


 /**********************************************************/

int main(int argc,char *argv[])
{ 
	int i, j, begl, endl;
   	if(argc<5) {
		printf("Error: Missing parameters!\n");
		printf("\tUsage: source-file begin-line end-line out-file \n");
		printf("\tNote: line no begin from 0\n\n"); 
		return 0; 
	}
	FILE *inf  = fopen(argv[1], "r");
	FILE *outf =fopen(argv[4], "w");
	begl = atoi(argv[2]);
	endl = atoi(argv[3]);

	i=0;
	while(!feof(inf) && i< begl && readline(inf, temp) ) i++;
	j = 0;
	while(!feof(inf) && i<= endl && readline(inf, temp) ){
		fprintf(outf,"%s\n", temp);
		i++; j++;
	}
	printf("read line = %d   read to = %d\n", j, i);
	fclose(inf);
	fclose(outf);
	return 0 ;
}












