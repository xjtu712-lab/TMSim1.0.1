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

int main(int argc,char *argv[]){

    int i,j,k,min,n, m, mark;
	char tnor, tref;
	int tindex, tdiff, tsame, tcover, tgeno;
	float maf;
    int * index;
	int * cover;

	if(argc<3) {printf("Error:Missing parameters.\nusage:cs matchout.txt cover.txt\n"); return 0;}

	FILE * inmatch = fopen(argv[1],"r");
	FILE * outcover = fopen(argv[2],"w");
    index = new int[100000000];
    cover = new int[100000000];
	
 	n=0;m=0;
	while( fscanf(inmatch,"%d %d %d %d %c %c %d", &tindex, &tdiff, &tsame, &tcover, &tnor, &tref, &tgeno) != EOF ){
			
			if(n%1000== 0)  {
				cover[m]=tcover;
                m++;
			}
			n++;
	}

	for(i=0;i<m/100;i++){
		for(j=i*100; j<i*100+100-1; j++){
			for(k=j+1; k<i*100+100; k++){
				if(cover[j] > cover[k]){  min=cover[k]; cover[k]=cover[j]; cover[j]=min;}
			}
			        
		}
	   fprintf(outcover,"%5d\n", cover[i*100+50]);
	}
	fclose(inmatch);
	fclose(outcover);
    delete []index;
	delete []cover;
	return 0;

}
