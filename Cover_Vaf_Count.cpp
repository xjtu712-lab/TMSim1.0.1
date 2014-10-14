#include <stdio.h>
#include <stdlib.h>
#include <iostream> 
#include <string> 
#include <fstream>
#include <string.h>
using namespace std;

int covernum[2000];
int main(int argc,char * argv[]){
	int i,j,max,max_index,reof, n, m;
	const int LEN=200;
	char tnor, tref;
	int tindex, tdiff, tsame, tcover, tgeno;
    int freq[200];
	int k=0;
	float maf;

	if(argc<4) printf("Error:Missing parameters.\n");

	FILE * inmatch = fopen(argv[1],"r");
	FILE * outcover = fopen(argv[2],"w");
	FILE * outmaf = fopen(argv[3],"w");
 	max=0;
	for(i=0; i<2000; i++) covernum[i] = 0;
	for(i=0; i<100; i++) freq[i] = 0;
	i=0; n =0; m = 0;
	tcover = 0;
	while( fscanf(inmatch,"%d %d %d %d %c %c %d", &tindex, &tdiff, &tsame, &tcover, &tnor, &tref, &tgeno) != EOF ){
			covernum[tcover] ++;
			if(tcover > 0){
				maf = 1.0*tdiff/tcover;
	            j = (int)(maf*40+0.5);
				if (j<=100)	freq[j] += 1;
				m++;
			}
			tcover = 0;
	}
	fclose(inmatch);
	
	max=0;
	max_index =0;
    for(n=0; n<2000; n++){  
        if(max < covernum[n]){
        	max = covernum[n]; 
			max_index = n;
		}
		fprintf(outcover,"%5d\t %10d\n", n, covernum[n]);
     }
    
	for(i=0;i<100;i++) fprintf(outmaf,"%5.2f \t %10d\n",1.0*i/40,freq[i]);
		
	fclose(outcover);
	fclose(outmaf);
	return 0;
}
