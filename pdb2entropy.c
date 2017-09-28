#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_ANG_PER_RESIDUE 10
#define HUGE_INT 1000000
#define MAX_N_RES 1000
#define MAX_N_TORS 50
#define MAX_N_NEXT 3
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_set_num_threads(num_threads) 0
#endif



struct Entropy {
int n_single;
int n_pair;
int n_nn;
double **h1;
int *i1;
int *i2;
double **h2;
double **mi;
double *total;
} Entropy;

struct System {
int n_atoms;
struct Atom *atoms;
int n_residues;
struct Residue *residues;
int n_chains;
struct Chain *chains;
int n_segments;
struct Segment *segments;
int n_models;
struct Model *models;
} System;

struct Next_def
{
char res[5];
char at1[5], at2[5];
double d;
} Next_def;

struct Tors_def
{
char res[5];
char tors_name[8];
char at1[5], at2[5], at3[5], at4[5];
double per;
} Tors_def;

struct Defs {
int n_res;
char res[MAX_N_RES][5];
int bn[MAX_N_RES],en[MAX_N_RES];
int bt[MAX_N_RES],et[MAX_N_RES];
int n_next[MAX_N_RES];
struct Next_def next[MAX_N_RES][MAX_N_NEXT];
int n_tors[MAX_N_RES];
struct Tors_def tors[MAX_N_RES][MAX_N_TORS];
} Defs;

struct Atom {
char at_name[5];
char alt_loc;
char res_name[4];
char chain;
char element[3];
int model;
int at_n;
int res_n;
char res_ins;
double coor[3];
double occ,temp;
char segid[5];
char pdb_chrg[3];
} Atom;

struct Residue {
char res_name[4];
char chain;
char res_ins;
int res_n;
int model;
char seg_name[5];
int beg, end ;
int n_alt_loc;
int prev, next;
} Residue;

struct Chain { char ch_name;
int beg, end;
int model;
char seg_name[5];
} Chain;


struct Segment { char seg_name[5];
int beg, end;
int model;
} Segment;

struct Model { int model_n;
int beg, end;
} Model;

int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}
# 132 "../../pdb2entropy.c"
struct Flag_par {
        char file_in_pdb[120];
        char file_in_def[120];
               char file_out[120];
               char pdb_outfile[120];
               double cutoff;
               double minres;
               int n, flag_n, list;
               int res,mi,kmi;
               int nt;
               int rt;
               int wp;
               int verbose;
      } Flag_par;

struct Tors {
        double phi;
 char tors_name[8];
 char res_name[4];
 char chain;
 int model;
 int res_n;
 double per;
 char res_ins;
 char segid[5];
} Tors;


struct Tors_res4nn {
        int n_models;
        int n_ang;
        double **phi;
        double **v;
 char **tors_name;
 char res_name[4];
 double *per;
 int res_n;
        char res_ins;
        char chain;
 char segid[5];
        double ref;
} Tors_res4nn;

void read_PDB_atoms(FILE *fp1, int *n_atoms, struct Atom *atoms);
void read_atom_pdb(char *buf, struct Atom *atom);
void pdb2tors(struct System system, struct Defs defs, struct Tors_res4nn *tors_res);

double distv(double *r1, double *r2);
double dist_ang(double *a1,double *a2, int n, double *per);
void init_flag_par(struct Flag_par *flag_par);
void check_cmd_line(int argc, char *argv[], struct Flag_par *flag_par);
void print_info_flag_par(struct Flag_par flag_par);
int comp (const void * elem1, const void * elem2);
void copy_atom(struct Atom *atom_to, struct Atom *atom_from);
void make_system(struct System *system);
void eliminate_alt_loc(struct System *system);
int is_next(struct System system, int i, int j, struct Next_def next_def);
double torsionv(double *x1, double *x2, double *x3, double *x4);
void kruskal(struct Entropy *entropy, struct Flag_par flag_par);
int compedge(const void *elem1, const void *elem2);
void suppos(double **ref, double **mob, int len, double *t, double **R, double *rmsd);
void quat_to_rot(double *q, double **R);
void CalcQuarticCoeffs(double **coords1, double **coords2, int len, double *coeff);
void CalcQuarticCoeffs_M(double **coords1, double **coords2, int len, double *coeff, double **M);
double CoordsInnerProd(double **coords, int len);
double QCProot(double *coeff, double guess, double delta);
double eval_horn_NR_corrxn(double *c, double x);
double eval_horn_quart(double *c, double x);
double eval_horn_quart_deriv(double *c, double x);

FILE *file_open(char *fname,char *acc);


int main(int argc, char *argv[]) {
      FILE *fp_in_1, *fp_in_2, *fp_out_1;
      char buf[1024], tmpres[5];
      struct Flag_par flag_par;
      struct System system;
      struct Defs defs;
      struct Entropy entropy;
      int i, j,ii,jj,k, kk, l, m, n_tors,n_atoms_per_model,n_res_per_model,n_res_per_model_mi,n_tors_per_model,K, num_threads;
      double **phit,c,L;
      double *ent_k, *d, *d_mean, *ld_mean, *ent_k_tot, **h1, **h2,logdk;
      double **R, *t, **coords1,**coords2, rmsd;
      struct Tors *tors;
      struct Tors_res4nn *tors_res, *tors_mi, tors_mi2;


        init_flag_par(&flag_par);

        check_cmd_line(argc, argv, &flag_par);
        if(flag_par.nt < 1)
        num_threads = omp_thread_count();
        else num_threads = flag_par.nt;
        omp_set_num_threads(num_threads);

        print_info_flag_par(flag_par);
        srand48(1);

 fp_in_1 = file_open(flag_par.file_in_pdb,"r");
        system.n_atoms = 0;
        while(fgets(buf,120,fp_in_1) != NULL )
            if(!strncmp(buf,"ATOM  ",strlen("ATOM  "))) system.n_atoms++;
        printf("%i atoms found in file %s\n",system.n_atoms,flag_par.file_in_pdb);
        system.atoms=calloc(system.n_atoms, sizeof(struct Atom));
        if(system.atoms == NULL)
        {printf("I could not allocate memory for %i atoms... exiting...\n", system.n_atoms); exit(0);}
        rewind(fp_in_1);
        read_PDB_atoms(fp_in_1, &(system.n_atoms), system.atoms);
        fclose(fp_in_1);
        for(i = 0; i< system.n_atoms; i++)
        for(k = 0; k< 3; k++)
        system.atoms[i].coor[k] += (drand48() - 0.5) * 0.001;
        make_system(&system);


  fp_in_2 = file_open(flag_par.file_in_def,"r");
          for(i=0;i<MAX_N_RES;i++)
          defs.n_next[i] = defs.n_tors[i] = 0;
          defs.n_res = 0;

          while(fgets(buf,120,fp_in_2) != NULL )
            {
            sscanf(buf,"%*s %s", tmpres);
            for(i = 0; i<defs.n_res; i++)
               if(!strcmp(tmpres,defs.res[i])) break;
            if(i == defs.n_res)
              {
              strcpy(defs.res[defs.n_res],tmpres);
              defs.n_res++;
              }
            if(!strncmp(buf,"NEXT",strlen("NEXT"))) {
              sscanf(buf,"%*s %s %s %s %lf",
              defs.next[i][defs.n_next[i]].res,
              defs.next[i][defs.n_next[i]].at1,
              defs.next[i][defs.n_next[i]].at2,
              &(defs.next[i][defs.n_next[i]].d)
              );
              defs.n_next[i]++;
              }
            else if(!strncmp(buf,"TORS",strlen("TORS"))) {
              sscanf(buf,"%*s %s %s %s %s %s %s %lf",
              defs.tors[i][defs.n_tors[i]].res,
              defs.tors[i][defs.n_tors[i]].tors_name,
              defs.tors[i][defs.n_tors[i]].at1,
              defs.tors[i][defs.n_tors[i]].at2,
              defs.tors[i][defs.n_tors[i]].at3,
              defs.tors[i][defs.n_tors[i]].at4,
              &(defs.tors[i][defs.n_tors[i]].per)
              );

              defs.tors[i][defs.n_tors[i]].per = 360.0/defs.tors[i][defs.n_tors[i]].per;
              defs.n_tors[i]++;
              }
            }
          fclose(fp_in_2);
printf("%i residues read in %s\n", defs.n_res, flag_par.file_in_def);





         if(system.n_residues % system.n_models == 0 && system.n_atoms % system.n_models == 0)
               {
               n_res_per_model = system.n_residues / system.n_models;
               n_atoms_per_model = system.n_atoms / system.n_models;
               tors_res=calloc(n_res_per_model, sizeof(Tors_res4nn));
               printf("OK: n. atoms and residues modulo n_models == 0\n");
               printf("n. residues per model = %i\n",n_res_per_model);
               }
          else {printf("n. atoms (%i) and residues (%i) modulo n_models (%i) != 0... exiting\n",system.n_atoms, system.n_residues,system.n_models); exit(0);}

     for(i=0; i< n_res_per_model; i++)
      {
       for(k=0; k<defs.n_res; k++)
       if(!strcmp(system.residues[i].res_name,defs.res[k])) break;
       if(k < defs.n_res)
        {
        tors_res[i].n_models = system.n_models;
 strcpy(tors_res[i].res_name, system.residues[i].res_name);
 strcpy(tors_res[i].segid, system.residues[i].seg_name);
        tors_res[i].res_ins = system.residues[i].res_ins;
        tors_res[i].chain = system.residues[i].chain;
        tors_res[i].res_n = system.residues[i].res_n;
 tors_res[i].phi = calloc(defs.n_tors[k],sizeof(double *));
          for(j=0; j<defs.n_tors[k]; j++)
    tors_res[i].phi[j] = calloc(system.n_models,sizeof(double));
 tors_res[i].v = calloc(defs.n_tors[k],sizeof(double *));
          for(j=0; j<defs.n_tors[k]; j++)
    tors_res[i].v[j] = calloc(3,sizeof(double));
 tors_res[i].per = calloc(defs.n_tors[k],sizeof(double));
 tors_res[i].tors_name = calloc(defs.n_tors[k],sizeof(char *));
          for(j=0; j<defs.n_tors[k]; j++)
             tors_res[i].tors_name[j] = calloc(8, sizeof(char));
        }
       }




    R = malloc(3 * sizeof(double *));
    for(i=0; i <3; i++)
        R[i] = malloc(3 * sizeof(double));
        t = malloc(3 * sizeof(double));
    coords1 = malloc(n_atoms_per_model * sizeof(double *));
    coords2 = malloc(n_atoms_per_model * sizeof(double *));
    for(i=system.models[0].beg; i<= system.models[0].end; i++)
      coords1[i] = system.atoms[i].coor;
    if(flag_par.rt)
    for(j=1; j<system.n_models; j++)
      {
      for(i=system.models[j].beg; i<= system.models[j].end; i++)
      coords2[i-system.models[j].beg] = system.atoms[i].coor;
      suppos(coords1, coords2, n_atoms_per_model, t, R, &rmsd);
      }
     if(flag_par.wp)
     {
     fp_out_1=fopen(flag_par.pdb_outfile,"w");
     for(j=0; j<system.n_models; j++)
      {
      fprintf(fp_out_1,"MODEL  %5i\n",j);
      for(i=system.models[j].beg; i<= system.models[j].end; i++)
 fprintf(fp_out_1,"ATOM  %5i %4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2lf%6.2lf      %4s%2s%2s\n",
             system.atoms[i].at_n, system.atoms[i].at_name, system.atoms[i].alt_loc, system.atoms[i].res_name, system.atoms[i].chain,
             system.atoms[i].res_n, system.atoms[i].res_ins,system.atoms[i].coor[0], system.atoms[i].coor[1], system.atoms[i].coor[2],
             system.atoms[i].occ, system.atoms[i].temp, system.atoms[i].segid, system.atoms[i].element, system.atoms[i].pdb_chrg);
      fprintf(fp_out_1,"ENDMDL\n");
      }
     fclose(fp_out_1);
     }
     pdb2tors(system, defs, tors_res);

    free(system.atoms);
    free(system.residues);
    free(system.chains);
    free(system.models);
    free(system.segments);
fp_out_1 = file_open(flag_par.file_out,"w");

     if(flag_par.list)
{
for(m = 0; m< n_res_per_model; m++)
for(j = 0; j< tors_res[m].n_ang; j++)
{
fprintf(fp_out_1,"%4s %5i %c %c %4s %8s --- ",
tors_res[m].res_name,
tors_res[m].res_n,
tors_res[m].res_ins,
tors_res[m].chain,
tors_res[m].segid,
tors_res[m].tors_name[j]);
for(i = 0; i < tors_res[m].n_models; i++)
fprintf(fp_out_1,"%9.4lf ", tors_res[m].phi[j][i]);
fprintf(fp_out_1,"\n");
}
}



K = flag_par.n + 1;
if(K > system.n_models)
{
K = system.n_models;
flag_par.n = K - 1;
}
ent_k = calloc(K-1,sizeof(double));
ent_k_tot = calloc(K-1,sizeof(double));
d_mean = calloc(K,sizeof(double));
ld_mean = calloc(K,sizeof(double));
entropy.n_nn = flag_par.n;


if(flag_par.res && !flag_par.mi)
{
fprintf(fp_out_1,"          k   entropy     average n.n.  Residue number ins.\n");
fprintf(fp_out_1,"              (R units)   distance      chain segid\n");
fprintf(fp_out_1,"                          (radians)\n");

for(m = 0; m < n_res_per_model; m++)
if(tors_res[m].n_ang > 0)
{
phit = calloc(system.n_models, sizeof(double *));
for(i=0; i<system.n_models; i++)
{
phit[i] = calloc(tors_res[m].n_ang, sizeof(double));
for(j = 0; j< tors_res[m].n_ang; j++)
phit[i][j] = tors_res[m].phi[j][i];
}

for(i=0; i<K; i++)
d_mean[i] = ld_mean[i] = 0;
for(i=0; i< K-1; i++)
ent_k[i] = 0;
#pragma omp parallel for num_threads(num_threads) private(j,k,d,logdk) shared(ent_k, d_mean, ld_mean)
for(i=0; i<system.n_models; i++)
{
if(flag_par.verbose)
{
if(i == 0)
printf("STARTING with model %i for residue %i\n", i, m);
else if(i%1000 == 0 && i != 0)
printf("DOING model %i for residue %i\n", i,m);
}
d = calloc(system.n_models,sizeof(double));
for(j=0; j<system.n_models; j++)
{
d[j] = dist_ang(phit[i],phit[j],tors_res[m].n_ang, tors_res[m].per);
d[j] = d[j] * M_PI/180.0;
}
qsort(d,system.n_models,sizeof(double), comp);
for(k = 1; k<K; k++)
{
if(d[k] < flag_par.minres)
{
logdk = log(flag_par.minres);

}
else
logdk = log(d[k]);
#pragma omp critical
{
ent_k[k-1] = ent_k[k-1] + logdk;
d_mean[k] = d_mean[k] + d[k];
ld_mean[k] = ld_mean[k] + logdk;
}
}
free(d);
}

for(k = 1; k<K; k++)
ent_k[k-1] = ent_k[k-1] * ((double) tors_res[m].n_ang / (double) system.n_models);

for(k=0,c=0.0; k<tors_res[m].n_ang; k++)
c = c - log( tors_res[m].per[k] * M_PI /180.0);

c = c +((double) tors_res[m].n_ang) * log(M_PI)/2.0 - lgamma(1.0 + ((double) tors_res[m].n_ang)/2.0) + 0.5722 + log( (double) system.n_models);

for(k = 1,L=0; k <= K - 1; k++)
{
ent_k[k-1] = ent_k[k-1] + c - L;
L = L + 1.0/(double) k;
d_mean[k] = d_mean[k]/(double) system.n_models;
ld_mean[k] = ld_mean[k]/(double) system.n_models;
ent_k_tot[k-1] = ent_k_tot[k-1] + ent_k[k-1];
}
for(k = 0; k < K - 1; k++)
{
fprintf(fp_out_1,"ent_k %5i %9.2lf     %9.2lf --- %4s %5i %c %c %4s\n", k+1, ent_k[k], d_mean[k+1],
tors_res[m].res_name,
tors_res[m].res_n,
tors_res[m].res_ins,
tors_res[m].chain,
tors_res[m].segid);
}
}
fprintf(fp_out_1,"                k   Total entropy\n");
fprintf(fp_out_1,"                    (R units)\n");
for(k = 0; k<K-1; k++)
fprintf(fp_out_1,"Total entropy %5i %10.2lf\n", k+1, ent_k_tot[k]);
fclose(fp_out_1);
}
else if(flag_par.mi)
{


l = 0;
for(i = 0 ; i< n_res_per_model; i++)
{
l = l + (int) floor((double) tors_res[i].n_ang/ (double) flag_par.kmi);
if( (tors_res[i].n_ang % flag_par.kmi) != 0) l++;
}
n_res_per_model_mi = l;
entropy.n_single = n_res_per_model_mi;

entropy.h1 = calloc(n_res_per_model_mi, sizeof(double *));
for(l=0; l<n_res_per_model_mi; l++)
entropy.h1[l] = calloc(K, sizeof(double *));
tors_mi=calloc(n_res_per_model_mi, sizeof(Tors_res4nn));
for(i = 0 ; i< n_res_per_model_mi; i++)
{
tors_mi[i].n_models = system.n_models;
tors_mi[i].phi = calloc(flag_par.kmi,sizeof(double *));
tors_mi[i].v = calloc(flag_par.kmi,sizeof(double *));
tors_mi[i].per = calloc(flag_par.kmi,sizeof(double));
tors_mi[i].tors_name = calloc(flag_par.kmi,sizeof(char *));
for(j=0; j<flag_par.kmi; j++)
tors_mi[i].tors_name[j] = calloc(8, sizeof(char));
}

for(i = 0, j = 0; i< n_res_per_model; i++)
{
for(k = 0; k< tors_res[i].n_ang; k++)
{
if(flag_par.verbose)
printf("res_mi %i - res %i - k %i\n", j, i, k+1);
strcpy(tors_mi[j].res_name, tors_res[i].res_name);
strcpy(tors_mi[j].segid, tors_res[i].segid);
strcpy(tors_mi[j].tors_name[k%flag_par.kmi], tors_res[i].tors_name[k]);
tors_mi[j].res_ins = tors_res[i].res_ins;
tors_mi[j].chain = tors_res[i].chain;
tors_mi[j].res_n = tors_res[i].res_n;
tors_mi[j].per[k%flag_par.kmi] = tors_res[i].per[k];
tors_mi[j].phi[k%flag_par.kmi] = tors_res[i].phi[k];
tors_mi[j].v[k%flag_par.kmi] = tors_res[i].v[k];
tors_mi[j].n_ang = k%flag_par.kmi + 1;
if((k%flag_par.kmi) == (flag_par.kmi - 1) || k == (tors_res[i].n_ang - 1)) j++;
}
}
for(m = 0; m < n_res_per_model_mi; m++)
{
fprintf(fp_out_1,"group %4i: %4s %5i %c %c %4s ", m,
tors_mi[m].res_name,
tors_mi[m].res_n,
tors_mi[m].res_ins,
tors_mi[m].chain,
tors_mi[m].segid);
{
for(j = 0; j< tors_mi[m].n_ang; j++)
fprintf(fp_out_1,"%8s ",tors_mi[m].tors_name[j]);
fprintf(fp_out_1,"\n");
}
}
fprintf(fp_out_1,"          k  entropy   average n.n.  group Residue number ins.\n");
fprintf(fp_out_1,"             (R units) distance        n.  chain segid\n");
fprintf(fp_out_1,"                       (radians)\n");

phit = calloc(system.n_models, sizeof(double *));
for(m = 0; m < n_res_per_model_mi; m++)
{
for(i=0; i<system.n_models; i++)
{
phit[i] = calloc(tors_mi[m].n_ang, sizeof(double));
for(j = 0; j< tors_mi[m].n_ang; j++)
phit[i][j] = tors_mi[m].phi[j][i];
}
for(i=0; i<K-1; i++)
ent_k[i] = 0;
for(i=0; i<K; i++)
d_mean[i] = ld_mean[i] = 0;

#pragma omp parallel for num_threads(num_threads) private(j,k,d,logdk) shared(ent_k, d_mean, ld_mean)
for(i=0; i<system.n_models; i++)
{
if(flag_par.verbose)
{
if(i%1000 == 0 && i == 0)
printf("STARTING with model %i for residue %i\n", i,m);
else if(i%1000 == 0 && i != 0)
printf("DOING model %i for residue %i\n", i,m);
}
d = calloc(system.n_models,sizeof(double));
for(j=0; j<system.n_models; j++)
{
d[j] = dist_ang(phit[i],phit[j],tors_mi[m].n_ang, tors_mi[m].per);
d[j] = d[j] * M_PI/180.0;
}
qsort(d,system.n_models,sizeof(double), comp);


for(k = 1; k<K; k++)
{
if(d[k] < flag_par.minres)
logdk = log(flag_par.minres);
else
logdk = log(d[k]);
#pragma omp critical
{
ent_k[k-1] = ent_k[k-1] + logdk ;
d_mean[k] = d_mean[k] + d[k];
ld_mean[k] = ld_mean[k] + logdk;
}
}
free(d);
}

for(k = 1; k<K; k++)
ent_k[k-1] = ent_k[k-1] * ((double) tors_mi[m].n_ang / (double) system.n_models);

for(k=0,c=0.0; k<tors_mi[m].n_ang; k++)
c = c - log( tors_mi[m].per[k] * M_PI /180.0);

c = c +((double) tors_mi[m].n_ang) * log(M_PI)/2.0 - lgamma(1.0 + ((double) tors_mi[m].n_ang)/2.0) + 0.5722 + log( (double) system.n_models);

for(k = 1,L=0; k < K; k++)
{
ent_k[k-1] = ent_k[k-1] + c - L;
L = L + 1.0/(double) k;
d_mean[k] = d_mean[k]/(double) system.n_models;
ld_mean[k] = ld_mean[k]/(double) system.n_models;
ent_k_tot[k-1] = ent_k_tot[k-1] + ent_k[k-1];
}
for(k = 0; k < K - 1; k++)
{
entropy.h1[m][k] = ent_k[k];
fprintf(fp_out_1,"ent_k %5i %9.2lf %9.2lf --- %5i %4s %5i %c %c %4s\n", k+1, ent_k[k], d_mean[k+1], m,
tors_mi[m].res_name,
tors_mi[m].res_n,
tors_mi[m].res_ins,
tors_mi[m].chain,
tors_mi[m].segid);
}
}


tors_mi2.n_models = system.n_models;
tors_mi2.phi = calloc(2*flag_par.kmi,sizeof(double *));
tors_mi2.per = calloc(2*flag_par.kmi,sizeof(double));
phit = calloc(system.n_models, sizeof(double *));
for(i = 0; i<system.n_models; i++)
phit[i] = calloc(2*flag_par.kmi,sizeof(double));


for(i = 0; i< n_res_per_model_mi; i++)
for(j = i+1; j< n_res_per_model_mi; j++)
{
for(l = 0; l<tors_mi[i].n_ang; l++)
for(m = 0; m<tors_mi[j].n_ang; m++)
if(distv(tors_mi[i].v[l], tors_mi[j].v[m]) <= flag_par.cutoff)
{
l = tors_mi[i].n_ang + 1;
m = tors_mi[j].n_ang + 1;
entropy.n_pair++;
}
}
printf("entropy.n_pair %i \n", entropy.n_pair);
entropy.i1 = calloc(entropy.n_pair,sizeof(int));
entropy.i2 = calloc(entropy.n_pair,sizeof(int));
entropy.h2 = calloc(entropy.n_pair,sizeof(double *));
entropy.total = calloc(entropy.n_pair,sizeof(double));
for(l = 0; l<entropy.n_pair; l++)
entropy.h2[l] = calloc(K, sizeof(double));
entropy.mi = calloc(entropy.n_pair,sizeof(double *));
for(l = 0; l<entropy.n_pair; l++)
entropy.mi[l] = calloc(K, sizeof(double));
kk = 0;
fprintf(fp_out_1,"          k entropy    average n.n. group_1 group_2\n");
fprintf(fp_out_1,"            (R units)  distance\n");
fprintf(fp_out_1,"                       (radians)\n");

for(ii = 0; ii< n_res_per_model_mi; ii++)
for(jj = ii+1; jj< n_res_per_model_mi; jj++)
for(l = 0; l<tors_mi[ii].n_ang; l++)
for(m = 0; m<tors_mi[jj].n_ang; m++)
if(distv(tors_mi[ii].v[l], tors_mi[jj].v[m]) <= flag_par.cutoff)
{
entropy.i1[kk] = ii;
entropy.i2[kk] = jj;

l = tors_mi[ii].n_ang + 1;
m = tors_mi[jj].n_ang + 1;
tors_mi2.n_ang = tors_mi[ii].n_ang + tors_mi[jj].n_ang;
for(k=0,i = 0; k < tors_mi[ii].n_ang; k++)
{
tors_mi2.phi[i] = tors_mi[ii].phi[k];
tors_mi2.per[i++] = tors_mi[ii].per[k];
}
for(k=0; k < tors_mi[jj].n_ang; k++)
{
tors_mi2.phi[i] = tors_mi[jj].phi[k];
tors_mi2.per[i++] = tors_mi[jj].per[k];
}

for(i = 0; i<system.n_models; i++)
for(j = 0; j< tors_mi2.n_ang; j++)
phit[i][j] = tors_mi2.phi[j][i];

for(i=0; i<K; i++)
ent_k[i] = d_mean[i] = ld_mean[i] = 0;

#pragma omp parallel for num_threads(num_threads) private(j,k,d,logdk) shared(ent_k, d_mean, ld_mean)
for(i=0; i<system.n_models; i++)
{
if(flag_par.verbose)
{
if(i%1000 == 0 && i == 0)
printf("STARTING with model %i for residue pair %i / %i\n", i,kk,entropy.n_pair);
else if(i%1000 == 0 && i != 0)
printf("DOING model %i for residue %i\n", i,m);
}
d = calloc(system.n_models,sizeof(double));
for(j=0; j<system.n_models; j++)
{
d[j] = dist_ang(phit[i],phit[j],tors_mi2.n_ang, tors_mi2.per);
d[j] = d[j] * M_PI/180.0;
}
qsort(d,system.n_models,sizeof(double), comp);

for(k = 1; k<K; k++)
{
if(d[k] < flag_par.minres)
logdk = log(flag_par.minres);
else
logdk = log(d[k]);
#pragma omp critical
{
ent_k[k-1] = ent_k[k-1] + logdk;
d_mean[k] = d_mean[k] + d[k];
ld_mean[k] = ld_mean[k] + logdk;
}
}
free(d);
}






for(k = 1; k<K; k++)
ent_k[k-1] = ent_k[k-1] * ((double) tors_mi2.n_ang / (double) system.n_models);

for(k=0,c=0.0; k<tors_mi2.n_ang; k++)
c = c - log( tors_mi2.per[k] * M_PI /180.0);

c = c +((double) tors_mi2.n_ang) * log(M_PI)/2.0 - lgamma(1.0 + ((double) tors_mi2.n_ang)/2.0) + 0.5722 + log( (double) system.n_models);

for(k = 1,L=0; k <= K - 1; k++)
{
ent_k[k-1] = ent_k[k-1] + c - L;
L = L + 1.0/(double) k;
d_mean[k] = d_mean[k]/(double) system.n_models;
ld_mean[k] = ld_mean[k]/(double) system.n_models;
ent_k_tot[k-1] = ent_k_tot[k-1] + ent_k[k-1];
}

for(k = 0; k < K - 1; k++)
{
entropy.h2[kk][k] = ent_k[k];
entropy.mi[kk][k] = entropy.h2[kk][k] - entropy.h1[entropy.i1[kk]][k] - entropy.h1[entropy.i2[kk]][k];

fprintf(fp_out_1,"MI    %5i %9.2lf %9.2lf --- %5i %5i\n", k+1, -entropy.mi[kk][k], d_mean[k+1], entropy.i1[kk], entropy.i2[kk]);
}
kk++;
}

kruskal(&entropy, flag_par);
fprintf(fp_out_1,"                 k   Total entropy\n");
fprintf(fp_out_1,"                     (R units)\n");
for(k = 0; k<K-1; k++)
fprintf(fp_out_1,"Total entropy %5i %10.2lf \n", k+1, entropy.total[k]);
}

}
void check_cmd_line(int argc, char *argv[], struct Flag_par *flag_par)
{
 int i;
        char tmp[100];
        char extension[100];

 if(argc < 4)
 {
 printf("Usage:\n");
 printf("./pdb2entropy pdb_infile def_infile outfile [Options]\n");
 printf("Options:\n");
 printf("-n (max k neighbours for listing entropies (20 default))\n");
 printf("-mi (compute entropy from mist)\n");
 printf("-kmi k (compute entropy from mist approx. at 2k (default k 1))\n");
 printf("-c X (cutoff distance (Angstrom) for MI)\n");
 printf("-mr X (minimum resolution (radians) assumed to avoid log(0), 5e-4 default)\n");
 printf("-nt X (number of threads to be used, if less than 1 the program finds the number of threads available)\n");
 printf("-rt (superpose all structures to the first one)\n");
 printf("-wp pdb_file (write superimposed strcutures in pdb_file)\n");
 printf("-l (list computed angles)\n");
 printf("-v (verbose mode)\n");
 printf("\n");
 exit(1);
 }



        strcpy((*flag_par).file_in_pdb, argv[1]);
        strcpy((*flag_par).file_in_def, argv[2]);
        strcpy((*flag_par).file_out, argv[3]);


 for (i = 4; i < argc; i++) {
  if (!strncmp(argv[i],"-v",3)) (*flag_par).verbose = 1;
  else if (!strncmp(argv[i],"-c",3)) (*flag_par).cutoff = atoi(argv[++i]);
  else if (!strncmp(argv[i],"-mr",4)) (*flag_par).minres = atof(argv[++i]);
  else if (!strncmp(argv[i],"-mi",4)) (*flag_par).mi = 1;
  else if (!strncmp(argv[i],"-nt",4)) (*flag_par).nt = atoi(argv[++i]);
  else if (!strncmp(argv[i],"-kmi",5)) (*flag_par).kmi = atoi(argv[++i]);
  else if (!strncmp(argv[i],"-rt",4)) (*flag_par).rt = 1;
  else if (!strncmp(argv[i],"-wp",4))
                 {
                 flag_par->wp = 1;
                 strcpy(flag_par->pdb_outfile,argv[++i]);
                 }
  else if (!strncmp(argv[i],"-n",3))
                 {
                 (*flag_par).flag_n = 1;
                 (*flag_par).n = atoi(argv[++i]);
                 }
  else if (!strncmp(argv[i],"-l",3)) (*flag_par).list = 1;
  else
                {
                 printf("I don't know option %s\n", argv[i]);
                 exit(2);
                }
}
printf("\n########################################################\n\n");

}

void init_flag_par(struct Flag_par *flag_par)
{
(*flag_par).list=0;
(*flag_par).flag_n=1;
(*flag_par).nt=1;
(*flag_par).rt=0;
(*flag_par).wp=0;
(*flag_par).n=20;
(*flag_par).verbose=0;
(*flag_par).cutoff=6.0;
(*flag_par).minres=1e-10;
(*flag_par).mi=0;
(*flag_par).kmi=1;
(*flag_par).res=1;
}

void print_info_flag_par(struct Flag_par flag_par)
{
 printf("########################################################\n");
 printf("# ACCORDING TO YOUR PARAMETER CHOICES:                 #\n");
 printf("########################################################\n\n");

        printf("pdb file: %s\n", flag_par.file_in_pdb);
        printf("out file: %s\n", flag_par.file_out);
        if(flag_par.nt > 0)
        printf("I will use %i threads\n", flag_par.nt);
        else
        printf("I will use all threads available\n");
        printf("I will print entropy only for the first %i neighbours\n", flag_par.n);
        if(flag_par.list)
        printf("I will print torsion angles\n");
        if(flag_par.rt)
        printf("I will superimpose all structures on the first one\n");
        if(flag_par.rt && flag_par.wp)
        printf("I will write superimposed structures in file %s\n", flag_par.pdb_outfile);
        if(flag_par.mi)
        {
        printf("I will compute entropy based on single entropies and mutual information\n");
        printf("Torsions are grouped %i by %i\n", flag_par.kmi, flag_par.kmi);
        printf("Mutual information computed only for those torsions whose central bonds are closer than %lf\n", flag_par.cutoff);
        }
        else
        printf("I will compute entropy based on single residues\n");
        printf("To avoid zeros in the distances I assume %e radian minimum distance\n", flag_par.minres);
 if(flag_par.verbose)
  printf("I will be verbose\n");
 else
  printf("I will NOT be verbose\n");
 printf("########################################################\n\n");
}

void pdb2tors(struct System system, struct Defs defs, struct Tors_res4nn *tors_res)
{
int i, j, k, l, n, id, prev, next, n_res_per_model, n_atoms_per_model;
int i1, i2, i3, i4;
double *v1, *v2, *v3, *v4, phi, t[3];
char *at;

for(i=0; i<system.n_residues; i++)
{
       for(k=0; k<defs.n_res; k++)
       if(!strcmp(system.residues[i].res_name,defs.res[k])) break;
       id = k;
prev = next = -1;
if(i < system.n_residues - 1)
for(j=0; j< defs.n_next[k];j++)
if(is_next(system,i,i+1,defs.next[k][j]))
{

system.residues[i].next = i+1;
system.residues[i+1].prev = i;
}
}
n_res_per_model = system.n_residues/system.n_models;
n_atoms_per_model = system.n_atoms/system.n_models;
for(i=0; i<n_res_per_model; i++)
{
     tors_res[i].n_ang = 0;

       for(k=0; k<defs.n_res; k++)
       if(!strcmp(system.residues[i].res_name,defs.res[k])) break;
       id = k;
       for(j = 0; j < defs.n_tors[k]; j++)
       {
        v1=v2=v3=v4=NULL;

        if((defs.tors[k][j].at1[0] == '-'))
         {
         if(system.residues[i].prev == (i-1))
           for(l = system.residues[i-1].beg; l <= system.residues[i-1].end; l++)
              {
              at = defs.tors[k][j].at1 + 1;
              if(!strcmp(system.atoms[l].at_name,at))
              {
              v1 = system.atoms[l].coor; i1 = l;
              break;
              }
              }
         }
          else
           for(l = system.residues[i].beg; l <= system.residues[i].end; l++)
              {
              at = defs.tors[k][j].at1;
              if(!strcmp(system.atoms[l].at_name,at))
              {
              v1 = system.atoms[l].coor; i1 = l;
              break;
              }
              }
        if((defs.tors[k][j].at2[0] == '-'))
         {
         if(system.residues[i].prev == (i-1))
           for(l = system.residues[i-1].beg; l <= system.residues[i-1].end; l++)
              {
              at = defs.tors[k][j].at2 + 1;
              if(!strcmp(system.atoms[l].at_name,at))
              {
              v2 = system.atoms[l].coor; i2 = l;
              break;
              }
              }
          }
          else
           for(l = system.residues[i].beg; l <= system.residues[i].end; l++)
              {
              at = defs.tors[k][j].at2;
              if(!strcmp(system.atoms[l].at_name,at))
              {
              v2 = system.atoms[l].coor; i2 = l;
              break;
              }
              }
        if((defs.tors[k][j].at3[0] == '+'))
         {
         if(system.residues[i].next == (i+1))
           for(l = system.residues[i+1].beg; l <= system.residues[i+1].end; l++)
              {
              at = defs.tors[k][j].at3 + 1;
              if(!strcmp(system.atoms[l].at_name,at))
              {
              v3 = system.atoms[l].coor; i3 = l;
              break;
              }
              }
          }
          else
           for(l = system.residues[i].beg; l <= system.residues[i].end; l++)
              {
              at = defs.tors[k][j].at3;
              if(!strcmp(system.atoms[l].at_name,at))
              {
              v3 = system.atoms[l].coor; i3 = l;
              break;
              }

              }
        if((defs.tors[k][j].at4[0] == '+'))
         {
         if(system.residues[i].next == (i+1))
           for(l = system.residues[i+1].beg; l <= system.residues[i+1].end; l++)
              {
              at = defs.tors[k][j].at4 + 1;
              if(!strcmp(system.atoms[l].at_name,at))
              {
              v4 = system.atoms[l].coor; i4 = l;
              break;
              }
              }
          }
          else
           for(l = system.residues[i].beg; l <= system.residues[i].end; l++)
              {
              at = defs.tors[k][j].at4;
              if(!strcmp(system.atoms[l].at_name,at))
              {
              v4 = system.atoms[l].coor; i4 = l;
              break;
              }
              }
         if ( (v1 != NULL) && (v2 != NULL) && (v3 != NULL) && (v4 != NULL))
         {
         for(l = 0; l< system.n_models; l++)
           {
           v1 = system.atoms[i1 + l * n_atoms_per_model].coor;
           v2 = system.atoms[i2 + l * n_atoms_per_model].coor;
           v3 = system.atoms[i3 + l * n_atoms_per_model].coor;
           v4 = system.atoms[i4 + l * n_atoms_per_model].coor;
           tors_res[i].phi[tors_res[i].n_ang][l] = torsionv(v1,v2,v3,v4);
           for(n=0; n<3; n++)
           tors_res[i].v[tors_res[i].n_ang][n] += ((v2[n] + v3[n])*0.5);
           }
           for(n=0; n<3; n++)
           tors_res[i].v[tors_res[i].n_ang][n] = tors_res[i].v[tors_res[i].n_ang][n]/(double) system.n_models;
         strcpy(tors_res[i].tors_name[tors_res[i].n_ang],defs.tors[k][j].tors_name);
         tors_res[i].per[tors_res[i].n_ang] = defs.tors[k][j].per;
         tors_res[i].chain = system.residues[i].chain;
         tors_res[i].res_n = system.residues[i].res_n;
         tors_res[i].res_ins = system.residues[i].res_ins;
         strcpy(tors_res[i].segid, system.residues[i].seg_name);
         tors_res[i].n_ang++;
         }
         }

}
}

double dist_ang(double *a1,double *a2, int n, double *per)
{

double t = 0.0, d;
int i;






for(i = 0,t = 0; i<n; i++)
{
while(a1[i] > per[i]/2.0) a1[i] = a1[i] - per[i];
while(a1[i] < -per[i]/2.0) a1[i] = a1[i] + per[i];
while(a2[i] > per[i]/2.0) a2[i] = a2[i] - per[i];
while(a2[i] < -per[i]/2.0) a2[i] = a2[i] + per[i];
d = fabs(a1[i] - a2[i]);
if(d > per[i]/2.0) d = per[i] - d;




d = d*d;
t = t + d;
}
return sqrt(t);
}

int comp (const void * elem1, const void * elem2) {
    double f1 = *((double *)elem1);
    double f2 = *((double *)elem2);
    if (f1 > f2) return 1;
    if (f1 < f2) return -1;
    return 0;
}

FILE *file_open(char *fname,char *acc) {
    FILE *fp;
    fp =fopen(fname,acc);
    if (fp==NULL)
 {
        fprintf(stderr,"unable to open file \"%s\"... \n",fname);
        exit(1);
    }
    return(fp);
}

void read_PDB_atoms(FILE *fp1, int *n_atoms, struct Atom *atoms)
{
 char buf[120];
 int i=0;
        int mod_id=1;

 while(fgets(buf,120,fp1) != NULL )
     if(!strncmp("ATOM",buf,4))
            {
      atoms[i].model=mod_id;
   read_atom_pdb(buf, &atoms[i]);
            i++;
     }
     else
       if(!strncmp("ENDMDL",buf,6))
              mod_id++;
 *n_atoms = i;
}

void read_atom_pdb(char *buf, struct Atom *atom)
{

    char at_rec[5];
    char tok[10];

    strncpy(tok,buf,4);
    tok[4] = '\0';
    sscanf(tok,"%s", at_rec);
    if(strncmp("ATOM",at_rec,4))
    {
     printf("What is supposed to be an ATOM line does not start with string ATOM... exiting...\n");
        exit(1);
    }

    strncpy(tok,buf + 6,5);
    tok[5] = '\0';
    sscanf(tok,"%i",&(atom->at_n));

    strncpy(tok,buf + 12,4);
    tok[4] = '\0';
    sscanf(tok,"%s", atom->at_name);

    strncpy(tok,buf + 16,1);
    tok[1] = '\0';
    if(sscanf(tok,"%c", &(atom->alt_loc)) == -1) atom->alt_loc=' ';



 strncpy(tok,buf + 17,3);
    tok[3] = '\0';
    sscanf(tok,"%s", atom->res_name);

 strncpy(tok,buf + 21,1);
    tok[1] = '\0';
    if(sscanf(tok,"%c", &(atom->chain)) == EOF) atom->chain = ' ';

    strncpy(tok,buf + 22,4);
    tok[4] = '\0';
    sscanf(tok,"%i", &(atom->res_n));

 strncpy(tok,buf + 26,1);
    tok[1] = '\0';
    if (sscanf(tok,"%c", &(atom->res_ins)) == EOF) atom->res_ins=' ';

    strncpy(tok,buf + 30,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[0]));

 strncpy(tok,buf + 38,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[1]));

    strncpy(tok,buf + 46,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[2]));

    strncpy(tok,buf + 54,6);
    tok[6] = '\0';
    sscanf(tok,"%lf", &(atom->occ));

 strncpy(tok,buf + 60,6);
    tok[6] = '\0';
    sscanf(tok,"%lf", &(atom->temp));

    if(strlen(buf) > 76)
    {
       strncpy(tok,buf + 72,4);
        tok[4] = '\0';
  if (sscanf(tok,"%s", (atom->segid)) == EOF)
   strcpy(atom->segid,"    ");
    }
    else strcpy(atom->segid,"    ");


    if(strlen(buf) > 78)
 {
        strncpy(tok,buf + 76,2);
        tok[2] = '\0';
        if (sscanf(tok,"%s", (atom->element)) == EOF)
   strcpy(atom->element,"UN");
 }

    if(strlen(buf) > 80)
 {
        strncpy(tok,buf + 78,2);
        tok[2] = '\0';
        if (sscanf(tok,"%s", (atom->pdb_chrg)) == EOF)
   strcpy(atom->pdb_chrg,"  ");
 }

}



int cmp_atoms(const void *p1, const void *p2)
{
 struct Atom A_atom, B_atom;
 int check = 0 ;

 A_atom = *((struct Atom *)p1);
 B_atom = *((struct Atom *)p2);




 if( A_atom.model < B_atom.model) check = -1;
 else if (A_atom.model == B_atom.model)
    {
        check = 0;
     if( strcmp(A_atom.segid,B_atom.segid) < 0) check = -1;
     else if (!strcmp(A_atom.segid,B_atom.segid))
        {
         check = 0;
            if( (int) A_atom.chain < (int) B_atom.chain ) check = -1;
            else if( A_atom.chain == B_atom.chain )
            {
             check = 0;
                if( A_atom.res_n < B_atom.res_n ) check = -1;
                else if(A_atom.res_n == B_atom.res_n )
                {
                 check = 0;
                    if((int) A_atom.res_ins < (int) B_atom.res_ins ) check = -1;
                    else if( A_atom.res_ins == B_atom.res_ins )
                    {
                     check = 0;
                        if( strcmp(A_atom.at_name,B_atom.at_name) < 0) check = -1;
                        else if (!strcmp(A_atom.at_name,B_atom.at_name))
                        {
                        check = 0;
                        if( (int) A_atom.alt_loc < (int) B_atom.alt_loc ) check = -1;
                        else if( A_atom.alt_loc == B_atom.alt_loc ) check = 0;
                        else check = 1;
   }
                        else check = 1;
                    }
                        else check = 1;
                }
                        else check = 1;
            }
                        else check = 1;
        }
                        else check = 1;
 }
        else check = 1;
 return check;

}


void hpsort(struct Atom *ra, int n);


void make_system(struct System *system)
{
 int n_atoms = system->n_atoms;
 int *p_n_residues = &(system->n_residues);
 int *p_n_chains = &(system->n_chains);
 int *p_n_segments = &(system->n_segments);
 int *p_n_models = &(system->n_models);

 int imodel=-1, isegment=-1, ichain=-1, iresidue=-1;
 char res_ins = '*';
 char chain = '*';
 char segid[5] = "****";
 int model = -1;
 int res_n = -HUGE_INT;

 int i;

 hpsort(system->atoms, system->n_atoms);

        eliminate_alt_loc(system);

        printf("after eliminate_alt_loc %i atoms left\n", system->n_atoms);

 n_atoms = system->n_atoms;



             iresidue = ichain = isegment = imodel = 0;
 for(i=0; i<n_atoms; i++)
             {
  if(system->atoms[i].model != model)
  {
   imodel++;
   isegment++;
   ichain++;
   iresidue++;
                }
                else if(strcmp(system->atoms[i].segid,segid))
                        {
                                isegment++;
                                ichain++;
                                iresidue++;
                        }
                 else if(system->atoms[i].chain != chain)
                        {
                                ichain++;
                                iresidue++;
                        }
          else if(system->atoms[i].res_n != res_n)
                                iresidue++;
                 else if(system->atoms[i].res_ins != res_ins)
           iresidue++;

         model = system->atoms[i].model;
         chain = system->atoms[i].chain;
         res_n = system->atoms[i].res_n;
         res_ins = system->atoms[i].res_ins;
             strcpy(segid, system->atoms[i].segid);
                }

                 system->residues = calloc(iresidue,sizeof(struct Residue));
        if(system->residues == NULL)
        {printf("I could not allocate memory for %i residues... exiting...\n", iresidue); exit(0);}
                 system->chains = calloc(ichain,sizeof(struct Chain));
                 system->segments = calloc(isegment,sizeof(struct Segment));
                 system->models = calloc(imodel,sizeof(struct Model));

             iresidue = ichain = isegment = imodel = -1;
 res_ins = '*';
 chain = '*';
 strcpy(segid,"****");
 model = -1;
 res_n = -HUGE_INT;
 for (i=0; i<n_atoms; i++)
 {
  if(system->atoms[i].model != model)
  {
   imodel++;
   isegment++;
   ichain++;
   iresidue++;

   strcpy(system->residues[iresidue].res_name, system->atoms[i].res_name);
   system->residues[iresidue].res_n = system->atoms[i].res_n;
   system->residues[iresidue].res_ins = system->atoms[i].res_ins;
   system->residues[iresidue].chain = system->atoms[i].chain;
   system->residues[iresidue].model = system->atoms[i].model;
   strcpy(system->residues[iresidue].seg_name, system->atoms[i].segid);
   system->residues[iresidue].beg = i;
   if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
         chain = system->atoms[i].chain;
         res_n = system->atoms[i].res_n;
         res_ins = system->atoms[i].res_ins;
             strcpy(segid, system->atoms[i].segid);

   model=system->atoms[i].model;

   system->chains[ichain].ch_name = chain;
   system->chains[ichain].beg = i;
                  system->chains[ichain].model=system->atoms[i].model;
   strcpy(system->chains[ichain].seg_name, system->atoms[i].segid);
   if (ichain != 0) system->chains[ichain - 1].end = i-1;

   strcpy(system->segments[isegment].seg_name, segid);
         system->segments[isegment].model=system->atoms[i].model;
         system->segments[isegment].beg = i;
    if (isegment != 0) system->segments[isegment - 1].end = i-1;

         system->models[imodel].beg=i;
         system->models[imodel].model_n=model;
         if (imodel != 0) system->models[imodel - 1].end =i-1;



  }
  else if(strcmp(system->atoms[i].segid,segid))
   {
    isegment++;
      ichain++;
      iresidue++;

    strcpy(system->residues[iresidue].res_name, system->atoms[i].res_name);
    system->residues[iresidue].res_n = system->atoms[i].res_n;
    system->residues[iresidue].res_ins = system->atoms[i].res_ins;
    system->residues[iresidue].chain = system->atoms[i].chain;
    system->residues[iresidue].model = system->atoms[i].model;
    strcpy(system->residues[iresidue].seg_name, system->atoms[i].segid);
    system->residues[iresidue].beg = i;

    if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
       chain = system->atoms[i].chain;
          res_n = system->atoms[i].res_n;
           res_ins = system->atoms[i].res_ins;



              strcpy(segid, system->atoms[i].segid);

    system->chains[ichain].ch_name = chain;
    system->chains[ichain].beg = i;
          system->chains[ichain].model=system->atoms[i].model;
    strcpy(system->chains[ichain].seg_name, system->atoms[i].segid);
    if (ichain != 0) system->chains[ichain - 1].end = i-1;

         strcpy(system->segments[isegment].seg_name, segid);
          system->segments[isegment].model=system->atoms[i].model;
           system->segments[isegment].beg = i;
     if (isegment != 0) system->segments[isegment - 1].end = i-1;



   }
   else if(system->atoms[i].chain != chain)
    {
     ichain++;
     iresidue++;

     strcpy(system->residues[iresidue].res_name, system->atoms[i].res_name);
     system->residues[iresidue].res_n = system->atoms[i].res_n;
     system->residues[iresidue].res_ins = system->atoms[i].res_ins;
     system->residues[iresidue].chain = system->atoms[i].chain;
     system->residues[iresidue].model = system->atoms[i].model;
     strcpy(system->residues[iresidue].seg_name, system->atoms[i].segid);
     system->residues[iresidue].beg = i;
     if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
           chain = system->atoms[i].chain;
            res_n = system->atoms[i].res_n;
           res_ins = system->atoms[i].res_ins;
     system->chains[ichain].ch_name = chain;
     system->chains[ichain].beg = i;
           system->chains[ichain].model=system->atoms[i].model;
     strcpy(system->chains[ichain].seg_name, system->atoms[i].segid);
     if (ichain != 0) system->chains[ichain - 1].end = i-1;


    }
    else if(system->atoms[i].res_n != res_n)
     {
     iresidue++;

     strcpy(system->residues[iresidue].res_name, system->atoms[i].res_name);
     res_n = system->atoms[i].res_n;
     res_ins = system->atoms[i].res_ins;
     system->residues[iresidue].res_n = res_n;
     system->residues[iresidue].res_ins = res_ins;
     system->residues[iresidue].chain = chain;
     system->residues[iresidue].model = system->atoms[i].model;
     strcpy(system->residues[iresidue].seg_name, system->atoms[i].segid);
     system->residues[iresidue].beg = i;

     if (iresidue != 0) system->residues[iresidue - 1].end = i-1;


     }
     else if(system->atoms[i].res_ins != res_ins)
     {
     iresidue++;
     strcpy(system->residues[iresidue].res_name, system->atoms[i].res_name);
     res_n = system->atoms[i].res_n;
     res_ins = system->atoms[i].res_ins;
     system->residues[iresidue].res_n = res_n;
     system->residues[iresidue].res_ins = res_ins;
     system->residues[iresidue].chain = chain;
     system->residues[iresidue].model = system->atoms[i].model;
     strcpy(system->residues[iresidue].seg_name, system->atoms[i].segid);
     system->residues[iresidue].beg = i;
     if (iresidue != 0) system->residues[iresidue - 1].end = i-1;


     }
 }
if(n_atoms != 0)
    {
    system->segments[isegment].end = i-1;
 system->chains[ichain].end = i-1;
 system->residues[iresidue].end = i-1;
    system->models[imodel].end = i-1;
    }
 *p_n_chains = ichain+1;
 *p_n_segments = isegment+1;
 *p_n_residues = iresidue+1;
 *p_n_models = imodel+1;

 printf("########################################################\n");
 printf("# SYSTEM:                                              #\n");
 printf("########################################################\n\n");
 printf("atoms = %8i, residues = %8i, chains = %8i, segments = %8i, models = %8i\n\n", n_atoms, *p_n_residues, *p_n_chains,*p_n_segments,*p_n_models );
 printf("########################################################\n\n");
}


void hpsort(struct Atom *ra, int n)
{
    int N, i, parent, child;
    struct Atom rra;
    N = n;
    i = n/2;
    for (;;) {
        if (i > 0) {
            i--;
            copy_atom(&(rra),&(ra[i]));

        } else {
            n--;
            if (n == 0) return;
            copy_atom(&(rra),&(ra[n]));

            copy_atom(&(ra[n]),&(ra[0]));

        }

        parent = i;
        child = i*2 + 1;

        while (child < n) {

            if (child + 1 < n && (cmp_atoms(&(ra[child + 1]),&(ra[child])) > 0)) {
                child++;
            }
            if (cmp_atoms(&(ra[child]),&(rra)) > 0) {
                copy_atom(&(ra[parent]),&(ra[child]));

                parent = child;

                child = parent*2+1;
            } else {
                break;
            }
        }
       copy_atom(&(ra[parent]),&(rra));

    }
}


void copy_atom(struct Atom *atom_to, struct Atom *atom_from)
{
int i;
strcpy((*atom_to).at_name,(*atom_from).at_name);
(*atom_to).alt_loc = (*atom_from).alt_loc;
strcpy((*atom_to).res_name,(*atom_from).res_name);
(*atom_to).chain = (*atom_from).chain;
(*atom_to).model = (*atom_from).model;
(*atom_to).at_n = (*atom_from).at_n;
(*atom_to).res_n = (*atom_from).res_n;
(*atom_to).res_ins = (*atom_from).res_ins;
for(i=0;i<3;i++)
(*atom_to).coor[i] = (*atom_from).coor[i];
(*atom_to).occ = (*atom_from).occ;
(*atom_to).temp = (*atom_from).temp;
strcpy((*atom_to).segid,(*atom_from).segid);
strcpy((*atom_to).pdb_chrg,(*atom_from).pdb_chrg);
}

void eliminate_alt_loc(struct System *system)
{
int i,j;
char buf[120];
j = 1;
for (i=1; i< (*system).n_atoms; i++)
{
if(
strcmp((*system).atoms[i-1].at_name,(*system).atoms[i].at_name) ||
strcmp((*system).atoms[i-1].res_name,(*system).atoms[i].res_name) ||
strcmp((*system).atoms[i-1].segid,(*system).atoms[i].segid) ||
((*system).atoms[i-1].chain != (*system).atoms[i].chain) ||
((*system).atoms[i-1].model != (*system).atoms[i].model) ||
((*system).atoms[i-1].res_n != (*system).atoms[i].res_n) ||
((*system).atoms[i-1].res_ins != (*system).atoms[i].res_ins)
)
{
copy_atom(&((*system).atoms[j]), &((*system).atoms[i]));

(*system).atoms[j].alt_loc = ' ';
j++;
}
}
if ((*system).n_atoms > 0)
(*system).n_atoms = j;
else
{
printf("n. atoms less than 1.... exiting \n");
exit(0);
}
}


int is_next(struct System system, int i, int j, struct Next_def next_def)
{
int k, p = 0;
double *v1, *v2;
v1 = v2 = NULL;
if(
(system.residues[i].chain == system.residues[j].chain ) &&
(system.residues[i].model == system.residues[j].model ) &&
!strcmp(system.residues[i].seg_name,system.residues[j].seg_name )
)
{
for(k = system.residues[i].beg; k <= system.residues[i].end; k++)
if(!strcmp(system.atoms[k].at_name,next_def.at1)) v1 = system.atoms[k].coor;
for(k = system.residues[j].beg; k <= system.residues[j].end; k++)
if(!strcmp(system.atoms[k].at_name,next_def.at2)) v2 = system.atoms[k].coor;
}
if(v1 == NULL || v2 == NULL) p = 0;
else if(sqrt(
(v1[0] -v2[0])*(v1[0] -v2[0]) +
(v1[1] -v2[1])*(v1[1] -v2[1]) +
(v1[2] -v2[2])*(v1[2] -v2[2])) <= next_def.d) p = 1;
return p;
}

void scale(double *v, double a)
{
 v[0] = v[0] * a;
 v[1] = v[1] * a;
 v[2] = v[2] * a;
}

void diffv(double *v, double *r2, double *r1)
{
 v[0] = r2[0] - r1[0];
 v[1] = r2[1] - r1[1];
 v[2] = r2[2] - r1[2];
}


void sumv(double *v, double *r2, double *r1)
{
 v[0] = r2[0] + r1[0];
 v[1] = r2[1] + r1[1];
 v[2] = r2[2] + r1[2];
}


double dotv(double *x1, double *x2)
{
    double d;
    d = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
    return d;
}


double modv(double *x1)
{
    double d;
    d = sqrt(x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2]);
    return d;
}


void vecv(double *x, double *x1, double *x2)
{
    x[0] = x1[1]*x2[2]-x1[2]*x2[1];
    x[1] = x1[2]*x2[0]-x1[0]*x2[2];
    x[2] = x1[0]*x2[1]-x1[1]*x2[0];
}


double distv(double *r1, double *r2)
{
 double d,v[3];
 diffv(v,r2,r1);
 d = dotv(v, v );
 d = sqrt(d);
 return d;
}


double torsionv(double *x1, double *x2, double *x3, double *x4)
{
        double d21,d32,d43,x21[3],x32[3],x43[3],bxc[3];
        double ac,ab,bc,abxc, t;
        int i;


        diffv(x21, x2, x1);
        diffv(x32, x3, x2);
        diffv(x43, x4, x3);
        d21 = sqrt(dotv(x21,x21));
        d32 = sqrt(dotv(x32,x32));
        d43 = sqrt(dotv(x43,x43));

        for(i = 0; i < 3; i++)
        {
         x21[i] = x21[i]/d21;
         x32[i] = x32[i]/d32;
         x43[i] = x43[i]/d43;
        }

        ab = dotv(x21,x32);
        bc = dotv(x32,x43);
        ac = dotv(x21,x43);

        vecv(bxc, x32, x43);
        abxc=dotv(x21, bxc);

        t = 180.0 * atan2(abxc, -ac + ab*bc )/M_PI;
        return t;

}

struct node {
      int val;
 double weight;
};

struct edge {
  double weight;
     int u, v;
};

void kruskal(struct Entropy *entropy, struct Flag_par flag_par)
{
int i, j, k, kk, v1, v2, l,counts,r;
double weight;
  struct edge *edges =((struct edge *) malloc(entropy->n_pair*sizeof(struct edge)));
  struct edge *MST=((struct edge *) malloc(entropy->n_single*sizeof(struct edge)));
int set[entropy->n_single];
for(kk=0;kk< entropy->n_nn; kk++)
{
counts = 0;
entropy->total[kk] = 0.0;
  for (i = 0; i < entropy->n_pair; i++)
     {
     edges[i].u = (*entropy).i1[i];
     edges[i].v = (*entropy).i2[i];
     edges[i].weight= (*entropy).mi[i][kk];
     }


for(i = 0; i < entropy->n_single ; i++)
      (*entropy).total[kk] += (*entropy).h1[i][kk];
if(flag_par.verbose)
printf("total entropy single [%i] = %lf\n", kk, (*entropy).total[kk]);
  for(i = 0; i < entropy->n_single; i++){
     set[i]=i;
  }

  qsort(edges,entropy->n_pair,sizeof(struct edge), compedge);

  for( i = 0; i < entropy->n_pair; i++){

     j=edges[i].u;
     k=edges[i].v;

  if(set[j] != set[k]){

  MST[counts].u=edges[i].u;
  MST[counts].v=edges[i].v;
  MST[counts].weight=edges[i].weight;
  counts++;



r=set[k];
for ( l = 0; l < entropy->n_single; l++)
if ( set[l] == r) set[l] = set[j];
}
}



if(entropy->n_pair > 0)
  for(i = 0; i < entropy->n_single - 1; i++)
   {
   if(flag_par.verbose)
   printf("MST[%i]: %4i %4i %10.3lf\n", kk + 1, MST[i].u, MST[i].v, MST[i].weight);
      entropy->total[kk] += MST[i].weight;
   }
 }
}

int compedge (const void *elem1, const void *elem2) {
    struct edge f1 = *((struct edge *)elem1);
    struct edge f2 = *((struct edge *)elem2);
    if (f1.weight > f2.weight) return 1;
    if (f1.weight < f2.weight) return -1;
    return 0;
}


void suppos(double **ref, double **mob, int len, double *t, double **R, double *rmsd)
{
    double innerprod;
    double lambdamax;
    double coeff[3];
    double **x1, **x2, c1[3], c2[3];
    double **M;
    double M_copy[4][4];
    double ev[4], v[4], mod;
    double max;
    int i, j, k, found;
    int ind_max[4];


    M = malloc(4 * sizeof(double *));
    for(i=0; i <4; i++)
        M[i] = malloc(4 * sizeof(double));

    x1 = malloc(len * sizeof(double *));
    x2 = malloc(len * sizeof(double *));
    for(i=0; i <len; i++)
        {
        x1[i] = malloc(3 * sizeof(double));
        x2[i] = malloc(3 * sizeof(double));
        }

    for(j=0; j <3; j++)
    {
    c1[j] = c2[j] = 0.0;
    for(i=0; i <len; i++)
      {
      c1[j] = c1[j] + ref[i][j];
      c2[j] = c2[j] + mob[i][j];
      }
    c1[j] = c1[j] / (double) len;
    c2[j] = c2[j] / (double) len;
    for(i=0; i <len; i++)
      {
      x1[i][j] = ref[i][j] - c1[j];
      x2[i][j] = mob[i][j] - c2[j];
      }
    }
    for(j=0; j <3; j++)
    t[j] = c1[j] - c2[j];
    innerprod = CoordsInnerProd(x1, len) + CoordsInnerProd(x2, len);
    CalcQuarticCoeffs_M(x1, x2, len, coeff, M);

    lambdamax = QCProot(coeff, 0.5 * innerprod, 1e-3);
    mod = (innerprod - (2.0 * lambdamax))/(double) len;
    if(mod >= 0.0) *rmsd = sqrt(mod);
    else if(fabs(mod) < 1e-4) *rmsd = 0.0;




    for(i=0; i<4; i++)
    for(j=0; j<4; j++)
        M_copy[i][j] = M[i][j];

    for(j=0; j<4; j++)
    M[j][j] = M[j][j] - lambdamax;
# 1854 "../../pdb2entropy.c"
    for(j=1; j<4; j++) ind_max[j] = -1;
    for(j=1; j<4; j++)
      {
      max = 0;
      for(i=0; i<4; i++)
      if((fabs(M[i][j])) > max)
      {
      found = 0;
      for(k=1; k<=j; k++)
      if( ind_max[k] == i) found = 1;
      if(!found)
      {
      max = fabs(M[i][j]);
      ind_max[j] = i;

      }
      }

      for( k = 0; k < 4; k++)
      if(k != ind_max[j])
      {
      max = M[k][j];
      for(i = 0; i<4; i++)
      M[k][i] = M[k][i] - M[ind_max[j]][i] * max/M[ind_max[j]][j];
      }
      }

      for(j=3; j>0; j--)
        {
        ev[j] = -M[ind_max[j]][0]/M[ind_max[j]][j];
        for(k=1; k<j; k++)
           M[ind_max[k]][0] = M[ind_max[k]][0] - M[ind_max[k]][j] *ev[j];
        }
      ev[0] = 1.0;
      mod = 0.0;
      for(j=0; j<4; j++)
      mod += ev[j]*ev[j];
      mod = sqrt(mod);

      for(j=0; j<4; j++)
        ev[j] = ev[j]/mod;
# 1906 "../../pdb2entropy.c"
    quat_to_rot(ev,R);





    for(i=0; i <len; i++)
    {
    for(j=0; j < 3; j++)
     c2[j] = mob[i][j];
    for(j=0; j < 3; j++)
     mob[i][j] = 0.0;
    for(j=0; j < 3; j++)
    for(k=0; k < 3; k++)
    mob[i][j] = mob[i][j] + R[j][k]*c2[k];
    for(j=0; j < 3; j++)
    mob[i][j] = mob[i][j] + t[j];
    }
    for(i=0; i <len; i++)
        {
        free(x1[i]);
        free(x2[i]);
        }
    free(x1);
    free(x2);

    for(i=0; i <4; i++)
        free(M[i]);
    free(M);
}
# 2041 "../../pdb2entropy.c"
double CoordsInnerProd(double **coords, int len)
{
    int i;
    double sum, tmpx, tmpy, tmpz;

    sum = 0.0;
    for (i = 0; i < len; ++i)
    {
        tmpx = coords[i][0];
        tmpy = coords[i][1];
        tmpz = coords[i][2];
        sum += (tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
    }

    return(sum);
}

void quat_to_rot(double *q, double **R)
{
R[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
R[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
R[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
R[0][1] = 2.0*(q[1]*q[2] + q[0]*q[3]);
R[1][0] = 2.0*(q[1]*q[2] - q[0]*q[3]);
R[0][2] = 2.0*(q[1]*q[3] - q[0]*q[2]);
R[2][0] = 2.0*(q[1]*q[3] + q[0]*q[2]);
R[1][2] = 2.0*(q[2]*q[3] + q[0]*q[1]);
R[2][1] = 2.0*(q[2]*q[3] - q[0]*q[1]);
}

void CalcQuarticCoeffs(double **coords1, double **coords2, int len, double *coeff)
{
    double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
                    SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
                    SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
                    SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    double x1, x2, y1, y2, z1, z2;
    int i;

    Sxx = Sxy = Sxz = Syx = Syy = Syz = Szx = Szy = Szz = 0.0;
    for (i = 0; i < len; ++i)
    {
        x1 = coords1[i][0];
        y1 = coords1[i][1];
        z1 = coords1[i][2];
        x2 = coords2[i][0];
        y2 = coords2[i][1];
        z2 = coords2[i][2];

        Sxx += (x1 * x2);
        Sxy += (x1 * y2);
        Sxz += (x1 * z2);

        Syx += (y1 * x2);
        Syy += (y1 * y2);
        Syz += (y1 * z2);

        Szx += (z1 * x2);
        Szy += (z1 * y2);
        Szz += (z1 * z2);
    }

    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;

    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;

    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;



    coeff[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    coeff[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

    SxzpSzx = Sxz+Szx;
    SyzpSzy = Syz+Szy;
    SxypSyx = Sxy+Syx;
    SyzmSzy = Syz-Szy;
    SxzmSzx = Sxz-Szx;
    SxymSyx = Sxy-Syx;
    SxxpSyy = Sxx+Syy;
    SxxmSyy = Sxx-Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

    coeff[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
             + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
             + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
             + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));
}

void CalcQuarticCoeffs_M(double **coords1, double **coords2, int len, double *coeff, double **M)
{
    double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
                    SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
                    SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
                    SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    double x1, x2, y1, y2, z1, z2;
    int i;
    Sxx = Sxy = Sxz = Syx = Syy = Syz = Szx = Szy = Szz = 0.0;
    for (i = 0; i < len; ++i)
    {
        x1 = coords1[i][0];
        y1 = coords1[i][1];
        z1 = coords1[i][2];
        x2 = coords2[i][0];
        y2 = coords2[i][1];
        z2 = coords2[i][2];

        Sxx += (x1 * x2);
        Sxy += (x1 * y2);
        Sxz += (x1 * z2);

        Syx += (y1 * x2);
        Syy += (y1 * y2);
        Syz += (y1 * z2);

        Szx += (z1 * x2);
        Szy += (z1 * y2);
        Szz += (z1 * z2);
    }
M[0][0] = Sxx + Syy + Szz;
M[1][1] = Sxx - Syy - Szz;
M[2][2] = -Sxx + Syy - Szz;
M[3][3] = -Sxx - Syy + Szz;
M[0][1] = M[1][0] = Syz - Szy;
M[0][2] = M[2][0] = Szx - Sxz;
M[0][3] = M[3][0] = Sxy - Syx;
M[1][2] = M[2][1] = Sxy + Syx;
M[1][3] = M[3][1] = Szx + Sxz;
M[2][3] = M[3][2] = Syz + Szy;
    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;

    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;

    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;



    coeff[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    coeff[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

    SxzpSzx = Sxz+Szx;
    SyzpSzy = Syz+Szy;
    SxypSyx = Sxy+Syx;
    SyzmSzy = Syz-Szy;
    SxzmSzx = Sxz-Szx;
    SxymSyx = Sxy-Syx;
    SxxpSyy = Sxx+Syy;
    SxxmSyy = Sxx-Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

    coeff[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
             + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
             + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
             + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));
}


double QCProot(double *coeff, double guess, double delta)
{
    int i;
    double oldg;
    for (i = 0; i < 50; ++i)
    {
        oldg = guess;

        guess -= eval_horn_NR_corrxn(coeff, guess);

        if ( (fabs(guess - oldg) < fabs(delta*guess)) || fabs(guess+oldg) < 10*delta)
            return(guess);

    }

    fprintf(stderr,
            "\n\n ERROR21: Newton-Raphson root-finding in \'QCProot()\' did not converge \n");
    fprintf(stderr,
            "last rmsd in iteration: %f %f %f\n", oldg ,guess, delta);
    exit(EXIT_FAILURE);
}




double eval_horn_NR_corrxn(double *c, double x)
{
    double x2 = x*x;
    double b = (x2 + c[2])*x;
    double a = b + c[1];

    return((a*x + c[0])/(2.0*x2*x + b + a));
}



double eval_horn_quart(double *c, double x)
{
    return(((x*x + c[2])*x + c[1])*x + c[0]);
}




double eval_horn_quart_deriv(double *c, double x)
{
    return(2.0*(2.0*x*x + c[2])*x + c[1]);
}
