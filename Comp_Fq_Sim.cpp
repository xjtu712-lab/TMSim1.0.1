#include <stdio.h>
#include <stdlib.h>
#include <iostream> 
#include <string> 
#include <fstream>
#include <string.h>
using namespace std;

const int REF_MAX = 250000000;

typedef struct {
	int no;
	int diff, same, cover, gtype;
	char nor, ref;
} fq_t;

typedef struct {
	unsigned int u1,u2;
	int index;
} sim_t;


 /**********************************************************/

int main(int argc,char *argv[])
{ 

	int i,j,k,m, index, sim_len, fq_len;
	fq_t *fqdata; //[REF_MAX]; 
	sim_t *simdata;  //[REF_MAX];
	char refname[100];
	int ret_eof;

   	if(argc<4) { printf("Error:Missing parameters!\nusage: compfqsim  *.sim  fqout.txt fq_sim_matchout.txt\n"); return 0; }
	
	fqdata  = new fq_t[REF_MAX];
	simdata  = new sim_t[REF_MAX];

	FILE *insim = fopen(argv[1], "r");
	FILE *infq = fopen(argv[2], "r");
	FILE* outMatch = fopen(argv[3], "w");

	ret_eof = fscanf(insim, "%s\n", refname); 
	k = 0;
	while( ret_eof != EOF ) {
		ret_eof=fscanf(insim, "%d %u %u\n", &simdata[k].index, &simdata[k].u1, &simdata[k].u2);
		k++;
	}
	sim_len = k;
	fclose(insim);
	cout<<"sim_len = "<<sim_len<<endl;

	k = 0;
	while( (ret_eof=fscanf(infq, "%d %d %d %d %c %c %d\n", &fqdata[k].no, &fqdata[k].diff, &fqdata[k].same, &fqdata[k].cover, &fqdata[k].nor, &fqdata[k].ref, &fqdata[k].gtype) ) != EOF ) k++;
	fq_len = k;
	fclose(infq);
	cout<<"fq_len = "<<fq_len<<endl;

	for(k=0; k<sim_len; k++){
		m = 0;
		for(j = 0; j<fq_len; j++){
			if (simdata[k].index == fqdata[j].no) { 
				m = 1; 
				fqdata[j].gtype = -fqdata[j].gtype;
				break;
			}
		}
		if (m == 0) simdata[k].index = -simdata[k].index;
	}


	for(k=0; k<sim_len; k++){
		 if (simdata[k].index<0) fprintf(outMatch, "%d %u %u\n", simdata[k].index, simdata[k].u1, simdata[k].u2);
	}
	fprintf(outMatch, "==========================================================\n");
	for(k=0; k<fq_len; k++){
		if (fqdata[k].gtype > 0)	fprintf(outMatch, "%10d %3d %3d %3d %c %c %1d\n", 
				fqdata[k].no, fqdata[k].diff, fqdata[k].same, fqdata[k].cover, fqdata[k].nor, fqdata[k].ref, fqdata[k].gtype);
	}


	fclose(outMatch);
	delete [] fqdata;
	delete [] simdata;
	return 0 ;
}












