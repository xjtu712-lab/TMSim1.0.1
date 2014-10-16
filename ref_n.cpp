
//把参考序列中表识为N的碱基删掉
#include <stdio.h>
#include <stdlib.h>
#include <iostream> 
#include <string> 
#include <fstream>
#include <string.h>
using namespace std;

const int LEN = 141213431;


int main(int argc,char *argv[]){
 char *ref;
 ref=new char[LEN];
 int i,j,m,ref_len;
  char ref_name[50];
  char new_name[50];
	if(argc<2){printf("Error:Missing ref_name\n usage:n_ref *.fa\n");return 0;}
	
	for(int n=0;n<strlen(argv[1])-1;n++){
       if(argv[1][n]!='\0' && argv[1][n]!='.') ref_name[n]=argv[1][n];
		else {ref_name[n++]='\0';break;}
	}
	printf("refname:%s\n",ref_name);
	sprintf(new_name,"%s_n.fa",ref_name);
	printf("rnewname:%s\n",new_name);

 FILE * inref=fopen(argv[1],"r");
 FILE * outref=fopen(new_name,"w");
  
 fscanf(inref,"%s\n",ref);

 i=0;
 while(fscanf(inref, "%s\n", &ref[i]) != EOF ) {i=i+50;}
 ref_len=strlen(ref);

 m=0;
 for(j=0;j<ref_len;j++)
	if(_toupper(ref[j]) !='N') ref[m++]=ref[j];

 ref_len=m;

 fprintf(outref,">%s",ref_name);
 for(i=0;i<ref_len; i++){
   if ((i%50) == 0) fprintf(outref,"\n%c",ref[i]);
		else fprintf(outref,"%c",ref[i]);
 }

//delete ref;
fclose(inref);
fclose(outref);
delete []ref;
return 0;

}
