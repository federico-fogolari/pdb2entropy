/***** includes *********************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/***** local includes *********************/

#define MAX_SEQ_LEN 2048
char res_char[27] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
char res_name[26][4] =
 {
 "ALA",
 "BBB",
 "CYS",
 "ASP",
 "GLU",
 "PHE",
 "GLY",
 "HIS",
 "ILE",
 "JJJ",
 "LYS",
 "LEU",
 "MET",
 "ASN",
 "OOO",
 "PRO",
 "GLN",
 "ARG",
 "SER",
 "THR",
 "SEC",
 "VAL",
 "TRP",
 "XXX",
 "TYR",
 "ZZZ"
 };
double ent[26][6]=
 /*"ALA"*/ {{-2.4 , 0.5 , -0.6 , 0.0 , -1.0 , 0.0 },
 /*"BBB"*/ { 0.0 , 0.0 , 0.0 , 0.0 , 0.0, 0.0},
 /*"CYS"*/ {-3.6 , 4.7 , -2.1 , 1.5 , -2.3 , 1.3},
 /*"ASP"*/ {-4.5 , 5.5 , -2.7 , 2.6 , -2.9 , 2.4},
 /*"GLU"*/ {-5.2 , 11.3 , -3.5 , 7.4 , -4.0 , 6.7},
 /*"PHE"*/ {-4.9 , 7.0 , -3.2 , 3.6 , -3.7 , 3.1 },
 /*"GLY"*/ {-1.9 , 0.5 , -0.4 , 0.0 , -0.5 , 0.0 },
 /*"HIS"*/ {-4.4 , 8.5 , -2.7 , 4.6 , -3.1 , 4.0 },
 /*"ILE"*/ {-6.6 , 4.8 , -4.9 , 2.3 , -4.9 , 2.2},
 /*"JJJ"*/ {  0.0 , 0.0 , 0.0 , 0.0 , 0.0, 0.0},
 /*"LYS"*/ {-7.5 , 16.8 , -6.0 , 12.5 , -6.3 , 11.6},
 /*"LEU"*/ {-6.3 , 4.5 , -4.6 , 2.1 , -5.2 , 1.8 },
 /*"MET"*/ {-6.1 , 14.7 , -4.5 , 10.1 , -4.9 , 9.1 },
 /*"ASN"*/ { -4.7 , 6.6 , -3.0 , 3.3 , -3.1 , 3.1 },
 /*"OOO"*/ {   0.0 , 0.0 , 0.0 , 0.0 , 0.0, 0.0},
 /*"PRO"*/ {-0.8 , 0.0 , -0.8 , 0.0 , 0.0 , 0.0 },
 /*"GLN"*/ { -5.5 , 13.7 , -3.9 , 9.3 , -4.4 , 8.3 },
 /*"ARG"*/ {-6.9 , 18.2 , -5.3 , 13.6 , -5.7 , 12.7 },
 /*"SER"*/ { -4.6 , 8.0 , -3.0 , 4.2 , -3.2 , 3.9 },
 /*"THR"*/ {-5.1 , 7.2 , -3.4 , 3.7 , -3.6 , 3.5 },
 /*"SEC"*/ { -3.6 , 4.7 , -2.1 , 1.5 , -2.3 , 1.3},
 /*"VAL"*/ {-4.6 , 2.2 , -3.0 , 0.5 , -3.0 , 0.5 },
 /*"TRP"*/ {-4.8 , 9.2 , -3.1 , 5.1 , -4.1 , 23.5 },
 /*"XXX"*/ {   0.0 , 0.0 , 0.0 , 0.0 , 0.0, 0.0},
 /*"TYR"*/ {-5.4 , 17.2 , -3.9 , 12.0 , -4.2 , 10.9 },
 /*"ZZZ"*/ {    0.0 , 0.0 , 0.0 , 0.0 , 0.0, 0.0}};


struct Flag_par {  
	       char file_in_seq[120];
               char file_out_ent[120];
               int verbose;
	     } Flag_par;

void init_flag_par(struct Flag_par *flag_par);
void check_cmd_line(int argc, char *argv[], struct Flag_par *flag_par);
void print_info_flag_par(struct Flag_par flag_par);

/**** MAIN *****/

int main(int argc, char *argv[]) {
      FILE *fp;
      char buf[1024], buf1[10],buf2[10];
      struct Flag_par flag_par;
      int i, j, k, l, m, n;
      char seq[MAX_SEQ_LEN];
      double totent = 0;

        /* Initialize options and parameters for this program */
        init_flag_par(&flag_par);      
        /* Read options from comman line */
        check_cmd_line(argc, argv, &flag_par);      
        /* Print information on options and parameters */
        print_info_flag_par(flag_par);      
        
        /* system 1....*/
        read_seq(flag_par.file_in_seq,seq);
        n = strlen(seq);
        printf("sequence read:\n");
        for(i=0; i<strlen(seq); i++)
        {
        printf("%c",seq[i]);
        if((i+1)%50 == 0 || i == (n-1))
        printf("\n");
        }
        
        fp = fopen(flag_par.file_out_ent,"w");

        printf("__________________________________________________________\n\n");
        printf("Entropy of the unfolded protein (in R units) computed from\n");
        printf("table 1 in Fogolari et al., PLOS One, 10, e0132356, 2015\n");
        printf("__________________________________________________________\n\n");
        fprintf(fp,"__________________________________________________________\n\n");
        fprintf(fp,"Entropy of the unfolded protein (in R units) computed from\n");
        fprintf(fp,"table 1 in Fogolari et al., PLOS One, 10, e0132356, 2015\n");
        fprintf(fp,"__________________________________________________________\n\n");



        printf("Resn aa entropy (R units)\n");
        fprintf(fp,"Resn aa entropy (R units)\n");
        for(i=0; i<strlen(seq); i++)
        {
         for(j=0;j<26;j++)
          if(seq[i] == res_char[j]) break;
        if(j == 26) 
        {
        printf("%5i %c not found\n",(i+1), seq[i]);
        fprintf(fp,"%5i %c not found\n",(i+1), seq[i]);
        }
        else if(i == 0)
        {
        printf("%5i %c %8.4lf\n",(i+1), seq[i], ent[j][2]);
        fprintf(fp,"%5i %c %8.4lf\n",(i+1), seq[i], ent[j][2]);
        totent = totent + ent[j][2];
        }
        else if(i == (n-1))
        {
        printf("%5i %c %8.4lf\n",(i+1), seq[i], ent[j][4]);
        fprintf(fp,"%5i %c %8.4lf\n",(i+1), seq[i], ent[j][4]);
        totent = totent + ent[j][4];
        }
        else 
        {
        printf("%5i %c %8.4lf\n",(i+1), seq[i], ent[j][0]);
        fprintf(fp,"%5i %c %8.4lf\n",(i+1), seq[i], ent[j][0]);
        totent = totent + ent[j][0];
        }
        }
        printf("Total unfolded entropy (R units): %11.4lf\n",totent);
        fprintf(fp,"Total unfolded entropy (R units): %11.4lf\n",totent);

}
void check_cmd_line(int argc, char *argv[], struct Flag_par *flag_par)
{
	int i;

	if(argc < 3) 
	{
	printf("Usage:\n"); 
	printf("./seq2unfent file_seq_in file_ent_out  [Options]\n"); 
	printf("Options:\n"); 
	printf("-v (verbose mode)"); 
	printf("\n"); 
	exit(1);
	}
        
        

        strcpy((*flag_par).file_in_seq, argv[1]);
        strcpy((*flag_par).file_out_ent, argv[2]);


	for (i = 3; i < argc; i++) {
		if (!strncmp(argv[i],"-v",2)) (*flag_par).verbose = 1;
		else 
                {
                 printf("I don't know option %s\n", argv[i]);
                 exit(2);
                }
        }
	for (i = 3; i < argc; i++) {
		if (!strncmp(argv[i],"-v",2)) (*flag_par).verbose = 1;
		else 
                {
                 printf("I don't know option %s\n", argv[i]);
                 exit(2);
                }
        }

}

void init_flag_par(struct Flag_par *flag_par)
{
(*flag_par).verbose=0;
}

void print_info_flag_par(struct Flag_par flag_par)
{
	printf("########################################################\n"); 
	printf("# ACCORDING TO YOUR PARAMETER CHOICES:                 #\n"); 
	printf("########################################################\n\n"); 

        printf("seq_file_in: %s\n", flag_par.file_in_seq);
        printf("ent_file_out: %s\n", flag_par.file_out_ent);
        printf("\n");
	if(flag_par.verbose) 
		printf("I will be verbose\n");
	else
		printf("I will NOT be verbose\n");
	printf("########################################################\n\n"); 
}

read_seq(char *filename, char *seq)
{
FILE *fp;
char buf[1024];
int i, l=0;

fp = fopen (filename, "r");
        while(fgets(buf,120,fp) != NULL)
        {
        if (!(buf[0] == '>'))
        for(i = 0; i< strlen(buf); i++)
           if(isalpha(buf[i]))
           {
           seq[l++] = toupper(buf[i]);
           }
        }
seq[l]='\0';
fclose(fp);
}

