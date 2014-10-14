/*In tumor.sam,identify density area
*four parameters:tumor.sam;density.txt;begin position of density;end position of density */


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
int READ_SIZE = 100;

char Ref_Name[50];
char ref_str[50];
char ref[REF_MAX];
unsigned char * g_n_site;
unsigned char * g_t_site;

char ch0[TERM_LEN];
char *pos;

int reads_num,reads_match_num,ret_eof;



int readline(FILE *fp, char *s){
int j;
	j=0;
	while(!feof(fp) && (s[j] = fgetc(fp)) != '\n') j++;
	s[j]='\0';
	if (j==0 && feof(fp)) return -1;
	else return j;    
}

int splitline(char *s, char **term){
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
static  void match_density_sam(char *tumorfile, char *densityfile, int beg, int end){
	char c, num;	
	int  i,j, tno1, tno2, insize1, insize2, beg1, beg2;
	char s[TERM_LEN];
	char *tmp1, *tmp2, *term1[50], *term2[50];

	FILE *intursam = fopen(tumorfile,"r");
	FILE *outdensam = fopen(densityfile,"w");//density_area sam file
	tmp1 = new char [2048];
	tmp2 = new char [2048];
	for(i=0; i<50;i++) term1[i] = new char [1024];
	for(i=0; i<50;i++) term2[i] = new char [1024];
	i=0;
	while(!feof(intursam) && readline(intursam, tmp1)>=0 ){
		if (tmp1[0] == '@') fprintf(outdensam,"@%s\n",tmp1);
		else break;
	}
	reads_num=0;
	reads_match_num=0;
	while ( !feof(intursam) ){
		tno1 = splitline(tmp1, term1);
		if (readline(intursam, tmp2)<0) break;
		tno2 = splitline(tmp2, term2);
		insize1 = atoi(term1[8]);
		if (insize1<0) insize1 = -insize1;
		beg1 = atoi(term1[3])-1;

		insize2 = atoi(term2[8]);
		if (insize2<0) insize2 = -insize2;
		beg2 = atoi(term2[3])-1;
		reads_num++;

		if( ( (beg1 >= beg) && (beg1 <= end) || (beg1+READ_SIZE-1 >= beg) && (beg1+READ_SIZE-1 <= end)) && (insize1 < 10200) ){
			if( ( (beg2 >= beg) && (beg2 <= end) || (beg2+READ_SIZE-1 >= beg) && (beg2+READ_SIZE-1 <= end)) && (insize2 < 10200) ){
				reads_match_num++;
				fprintf(outdensam,"%s\n",tmp1);
				fprintf(outdensam,"%s\n",tmp2);
			}
		}
		if (reads_match_num > 100) break;

		if (readline(intursam, tmp1)<0) break;
	}
	cout<<endl<<"reads_match_num="<<reads_match_num<<endl;
    fclose(outdensam);
	fclose(intursam);
	delete [] tmp1;
	delete [] tmp2;
	for(i=0; i<50;i++) if (term1[i])delete term1[i];
	for(i=0; i<50;i++) if (term2[i])delete term2[i];

}


int main(int argc, char *argv[])
{ 
	int *n;
	int len, gtype;
	int i,j,k,c,m;
	int beg,end,dis;
	int mb,me;
	
	/*if (argc < 4)  {printf("Error:Missing parameters!\nusage: ds  refrence.fa  tumor.sam  density.sam \n");return 0;}

	FILE *inref = fopen( argv[1], "r"); 
	FILE *inNgentype = fopen("nor_genotype.txt","r");
	FILE *inTgentype = fopen("tum_genotype.txt","r");
	FILE *outdensity = fopen("density.txt","w");
	
	fscanf(inref,"%s\n",ch0);
	strcpy(Ref_Name, &ch0[1]);
	sprintf(ref_str, "%s_", Ref_Name);
	i=0;
	while(fscanf(inref,"%s\n",&ref[i])!=EOF){i=i+50;}
	fclose(inref);
	REF_LEN=strlen(ref);
	cout<<REF_LEN<<endl;

	g_n_site = new unsigned char[REF_LEN];
	g_t_site = new unsigned char[REF_LEN];
	n = new int[(int)(REF_LEN/1000+1)];
	for(k=0; k<REF_LEN; k++){g_n_site[k]=0; g_t_site[k]=0;}

	while(fscanf(inNgentype,"%10d %2d\n", &k, &gtype)!=EOF){ g_n_site[k] = gtype;}
	while(fscanf(inTgentype,"%10d %2d\n", &k, &gtype)!=EOF){ g_t_site[k] = gtype;}
	fclose(inNgentype);
	fclose(inTgentype);
		
	for(i=0; i<REF_LEN; i+=1000) {
		k=0;   		
		for(j=i;j<i+1000;j++){
		    if(j>=REF_LEN) break;
			if(g_n_site[j]==0 && g_t_site[j]) k++; 
		}
		len=(i+1)/1000;
		n[len]=k;
		//if (k>0) cout<<"n["<<len<<"]="<<k<<endl;
	}
	
	len=(REF_LEN)/1000+1;
 	m=0; mb = 0; me = 0;
	for(i=0;i<len;i++){
 		 if(n[i]>=2){
			beg=i;
			c=0;
			for(j=beg;j<len;j++,i++){
          		if(n[j]<2) { 
					if (!c) end=j; 
					c++;
					if (c>3) break;
				}else {c=0;}
			}
			dis=end-beg;
			if (dis>m) {m=dis;mb=beg; me=end;}
			fprintf(outdensity,"beg=[%10d] end=[%10d] dis=[%10d]\n",beg,end,dis);
			cout<<"beg="<<beg<<"   end="<<end<<"   dis="<<dis<<endl;

  		}
	}
	cout<<"mb="<<mb<<" me="<<me<<" m="<<m<<endl;
    fclose(outdensity);

	delete [] g_n_site;
	delete [] g_t_site;
	delete [] n;
*/
	//if (argc >= 4) match_density_sam(argv[2], argv[3], mb*1000, me*1000);

	match_density_sam(argv[2], argv[3], 10000*1000, 110000*1000);
	//else {printf("parameter ERRor. too samll!\n usage: ds  refrence.fa  tumor.sam density.sam ");return 0;}
   
	return 0 ;
}

