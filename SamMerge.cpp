/*merge many *.sam to one file.
*parameters:s1.sam;s2.sam;.....sn.sam;merge.sam*/

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

/*****************************************************************/

int main(int argc, char *argv[])
{ 
int i, m;
char c;
char ch[120];
char s[100], ref_str[50];
	FILE *insam = fopen(argv[2],"r");
	FILE *outmergsam = fopen(argv[argc-1],"w");//density_area sam file

	if(argc<4) {printf("Error : Missing parameters.\nusage: sm  ref_name *1.sam ... *n.sam all.sam\n"); return 0;}
	while(!feof(insam) && ((c=fgetc(insam))=='@') ){
		i=0;
		while( (c = fgetc(insam)) != '\n') { s[i++] = c; }	
		s[i]='\0';
		fprintf(outmergsam,"@%s \n",s);
	}
    fclose(insam);

	for(i=0;i<100;i++) ch[i]='I';
	ch[100]='\0';

	sprintf(ref_str, "%s_", argv[1]);
	reads_num=0;
	reads_match_num=0;
	ch0[0]=c;
	for(i=2;i<(argc-1);i++){
		insam = fopen(argv[i], "r");
		ret_eof = fscanf(insam, "%s", &ch0[1]);
		while (ret_eof != EOF){
			while ( !(pos=strstr(ch0, ref_str)) && (ret_eof != EOF) ) ret_eof = fscanf(insam, "%s", ch0);
			if (ret_eof == EOF) break;
			ret_eof=fscanf(insam, "%s %s %s %s %s %s %s %s %s ", ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9);
			m=atoi(ch8);
			if (m<0) m=-m;
			reads_num++;
			if (reads_num%1000000==0) {
				cout<<"readno="<<reads_num<<"  ch3= "<<ch3<<"  read_match_num="<<reads_match_num<<endl;
			}
			if( (m < 10200)&& strcmp(ch5,"100M")==0 ){
				reads_match_num++;
				fprintf(outmergsam,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",ch0,ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9);
				fprintf(outmergsam,"%s\n",ch);
			}
			if (reads_match_num > 100) break;
			ret_eof = fscanf(insam, "%s", ch0);
		}
		cout<<endl<<"reads_match_num="<<reads_match_num<<endl;
	    fclose(insam);
	}
	fclose(outmergsam);
	return 0 ;
}

