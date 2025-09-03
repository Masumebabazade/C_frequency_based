#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <errno.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define INITIAL_CAP 1024
#define KEY_MULT 1000000000ULL

// structure for sample label
typedef struct {
    char *iid;
    char *city;
    int city_idx;
} Sample;

// genotype data for samples
typedef struct {
    int n_samples;
    int n_positions;
    uint64_t *positions;       // encoded positions
    unsigned char **genotypes; // [n_samples][n_positions]
} GenotypeData;

// informative file structure
typedef struct {
    int idx_i; // city i index
    int idx_j; // city j index
    int n;     // number of rows
    uint64_t *positions; // encoded positions
    unsigned char *values; // genotype value 0/1/2
    double *p_i;
    double *p_j;
} InfoFile;

// mapping city names to indices
typedef struct CityMap {
    char *name;
    int idx;
    struct CityMap *next;
} CityMap;

static CityMap *city_map = NULL;
static int num_cities = 0;
static char **cities = NULL;

static void add_city(const char *name) {
    CityMap *curr = city_map;
    while (curr) {
        if (strcmp(curr->name, name) == 0) return;
        curr = curr->next;
    }
    CityMap *node = (CityMap*)malloc(sizeof(CityMap));
    node->name = strdup(name);
    node->idx = num_cities++;
    node->next = city_map;
    city_map = node;
    cities = (char**)realloc(cities, num_cities * sizeof(char*));
    cities[num_cities-1] = node->name;
}

static int city_index(const char *name) {
    CityMap *curr = city_map;
    while (curr) {
        if (strcmp(curr->name, name) == 0) return curr->idx;
        curr = curr->next;
    }
    return -1;
}

// parse labels file
static Sample *load_labels(const char *label_path, int *out_n) {
    FILE *fp = fopen(label_path, "r");
    if (!fp) { perror("open labels"); return NULL; }
    char line[1024];
    int cap = 128;
    Sample *samples = (Sample*)malloc(cap * sizeof(Sample));
    int n = 0;
    while (fgets(line, sizeof(line), fp)) {
        char *iid = strtok(line, "\t \n\r");
        char *city = strtok(NULL, "\t \n\r");
        if (!iid || !city) continue;
        if (n == cap) { cap *= 2; samples = (Sample*)realloc(samples, cap * sizeof(Sample)); }
        samples[n].iid = strdup(iid);
        samples[n].city = strdup(city);
        add_city(city);
        n++;
    }
    fclose(fp);
    for (int i=0;i<n;i++) {
        samples[i].city_idx = city_index(samples[i].city);
    }
    *out_n = n;
    return samples;
}

// encode chrom,pos to key
static uint64_t encode_pos(const char *chrom, const char *pos_str) {
    unsigned long chr=0;
    if (strcmp(chrom,"X")==0) chr=23;
    else if (strcmp(chrom,"Y")==0) chr=24;
    else if (strcmp(chrom,"MT")==0 || strcmp(chrom,"M")==0) chr=25;
    else chr = strtoul(chrom, NULL, 10);
    unsigned long pos = strtoul(pos_str, NULL, 10);
    return chr*KEY_MULT + pos;
}

// binary search in positions
static int find_position(uint64_t *arr, int n, uint64_t key) {
    int l=0,r=n-1;
    while (l<=r) {
        int m=(l+r)/2;
        if (arr[m]==key) return m;
        else if (arr[m]<key) l=m+1; else r=m-1;
    }
    return -1;
}

// load VCF genotypes for samples
static GenotypeData *load_vcf(const char *vcf_path, Sample *samples, int n_samples) {
    FILE *fp = fopen(vcf_path, "r");
    if (!fp) { perror("open vcf"); return NULL; }
    char *line = NULL; size_t len=0;
    ssize_t read;
    int *sample_map = (int*)malloc(n_samples*sizeof(int));
    for (int i=0;i<n_samples;i++) sample_map[i]=-1;
    GenotypeData *gd = (GenotypeData*)calloc(1,sizeof(GenotypeData));
    gd->n_samples = n_samples;
    gd->positions = NULL;
    gd->n_positions = 0;
    gd->genotypes = (unsigned char**)malloc(n_samples*sizeof(unsigned char*));
    for (int i=0;i<n_samples;i++) gd->genotypes[i]=NULL;
    int capacity = 0;
    while ((read = getline(&line,&len,fp)) != -1) {
        if (line[0]=='#') {
            if (strncmp(line, "#CHROM",6)==0) {
                // parse header to map samples
                char *token; int col=0; int idx=0;
                for (token=strtok(line, "\t\n"); token; token=strtok(NULL,"\t\n")) {
                    if (col>=9) {
                        for (int s=0;s<n_samples;s++) {
                            if (strcmp(samples[s].iid, token)==0) {
                                sample_map[s]=col-9;
                                break;
                            }
                        }
                    }
                    col++;
                }
            }
            continue;
        }
        // data line
        char *fields[1024];
        int nf=0;
        char *tok=strtok(line,"\t\n");
        while (tok && nf<1024) { fields[nf++]=tok; tok=strtok(NULL,"\t\n"); }
        if (nf<9) continue;
        uint64_t key=encode_pos(fields[0],fields[1]);
        if (gd->n_positions==capacity) {
            capacity = capacity?capacity*2:INITIAL_CAP;
            gd->positions = (uint64_t*)realloc(gd->positions, capacity*sizeof(uint64_t));
            for (int s=0;s<n_samples;s++) gd->genotypes[s]=(unsigned char*)realloc(gd->genotypes[s], capacity*sizeof(unsigned char));
        }
        gd->positions[gd->n_positions]=key;
        for (int s=0;s<n_samples;s++) {
            int col = sample_map[s];
            unsigned char val=3; // missing default
            if (col>=0 && col<nf) {
                char *gt = fields[col];
                // genotype like 0|1 or 0/0
                int a=gt[0]-'0';
                int b=gt[2]-'0';
                if (a>=0 && a<=2 && b>=0 && b<=2) val=a+b; else val=3;
            }
            gd->genotypes[s][gd->n_positions]=val;
        }
        gd->n_positions++;
    }
    free(line); fclose(fp); free(sample_map);
    return gd;
}

// load FST and compute Q matrix
static double **load_fst_Q(const char *fst_path) {
    int n=num_cities;
    double **hudson = (double**)malloc(n*sizeof(double*));
    double **Q = (double**)malloc(n*sizeof(double*));
    for (int i=0;i<n;i++) {
        hudson[i]=(double*)calloc(n,sizeof(double));
        Q[i]=(double*)calloc(n,sizeof(double));
    }
    FILE *fp=fopen(fst_path,"r");
    if(!fp){perror("open fst"); return NULL;}
    char line[1024];
    /* read and discard header */
    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        for (int i = 0; i < n; ++i) {
            free(hudson[i]);
            free(Q[i]);
        }
        free(hudson);
        free(Q);
        return NULL;
    }
    while(fgets(line,sizeof(line),fp)){
        char pop1[64], pop2[64];
        double val;
        if(sscanf(line,"%63s%63s%lf",pop1,pop2,&val)!=3) continue;
        int i=city_index(pop1); int j=city_index(pop2);
        if(i<0||j<0) continue;
        hudson[i][j]=hudson[j][i]=val;
    }
    fclose(fp);
    // compute Q matrix column-wise
    for(int j=0;j<n;j++){
        double mean=0, std=0;
        for(int i=0;i<n;i++) mean+=hudson[i][j];
        mean/=n;
        for(int i=0;i<n;i++) {
            double diff=hudson[i][j]-mean; std+=diff*diff; }
        std=sqrt(std/n);
        if(std==0) std=1;
        for(int i=0;i<n;i++){
            double z=(hudson[i][j]-mean)/std;
            double q=0.5*(1.0+erf(z/sqrt(2.0)));
            Q[i][j]=q;
        }
    }
    for(int i=0;i<n;i++) free(hudson[i]);
    free(hudson);
    return Q;
}

// load informative files
static InfoFile ***load_info_files(const char *dir_path) {
    int n=num_cities;
    InfoFile ***table=(InfoFile***)malloc(n*sizeof(InfoFile**));
    for(int i=0;i<n;i++){ table[i]=(InfoFile**)malloc(n*sizeof(InfoFile*)); for(int j=0;j<n;j++) table[i][j]=NULL; }
    DIR *dir=opendir(dir_path);
    if(!dir){perror("opendir infs"); return table;}
    struct dirent *ent;
    while((ent=readdir(dir))!=NULL){
        if(strncmp(ent->d_name,"inf_",4)!=0) continue;
        char *dot=strrchr(ent->d_name,'.');
        if(!dot) continue;
        // remove .csv
        char name[256]; strncpy(name, ent->d_name+4, sizeof(name));
        name[sizeof(name)-1]='\0';
        char *popj=strtok(name,".");
        char *popi=strtok(NULL,".");
        if(!popi||!popj) continue;
        int j=city_index(popj); int i=city_index(popi);
        if(i<0||j<0) continue;
        // construct path
        char path[512]; snprintf(path,sizeof(path),"%s/%s",dir_path,ent->d_name);
        FILE *fp=fopen(path,"r"); if(!fp){perror("open inf"); continue;}
        char line[1024];
        /* discard header line */
        if (!fgets(line, sizeof(line), fp)) { fclose(fp); continue; }
        int cap=1024; int nrow=0;
        uint64_t *positions=(uint64_t*)malloc(cap*sizeof(uint64_t));
        unsigned char *values=(unsigned char*)malloc(cap*sizeof(unsigned char));
        double *p_i=(double*)malloc(cap*sizeof(double));
        double *p_j=(double*)malloc(cap*sizeof(double));
        while(fgets(line,sizeof(line),fp)){
            char pos[64]; int val; double pi,pj;
            if(sscanf(line,"%63[^,],%d,%lf,%lf",pos,&val,&pi,&pj)!=4) continue;
            if(nrow==cap){cap*=2; positions=(uint64_t*)realloc(positions,cap*sizeof(uint64_t)); values=(unsigned char*)realloc(values,cap*sizeof(unsigned char)); p_i=(double*)realloc(p_i,cap*sizeof(double)); p_j=(double*)realloc(p_j,cap*sizeof(double));}
            char *chr=strtok(pos,"_");
            char *p=strtok(NULL,"_");
            if(!chr||!p) continue;
            positions[nrow]=encode_pos(chr,p);
            values[nrow]=(unsigned char)val;
            p_i[nrow]=pi; p_j[nrow]=pj;
            nrow++;
        }
        fclose(fp);
        InfoFile *info=(InfoFile*)malloc(sizeof(InfoFile));
        info->idx_i=i; info->idx_j=j; info->n=nrow;
        info->positions=positions; info->values=values; info->p_i=p_i; info->p_j=p_j;
        table[j][i]=info; // note orientation j vs i as in file name inf_{popj}.{popi}
    }
    closedir(dir);
    return table;
}

// processing per individual
static void process_individual(int idx, Sample *samples, GenotypeData *gd, InfoFile ***info_table, double **Q,
                               double *sum_res, double *count_res, double *sum_q_res, double *count_q_res,
                               int *has_match, int *correct_sum, int *correct_count, int *correct_sum_q, int *correct_count_q) {
    int ncity=num_cities;
    double *sum_cols = (double*)calloc(ncity,sizeof(double));
    double *count_cols = (double*)calloc(ncity,sizeof(double));
    double *sum_q_cols = (double*)calloc(ncity,sizeof(double));
    double *count_q_cols = (double*)calloc(ncity,sizeof(double));

    for(int i=0;i<ncity;i++){
        for(int j=0;j<ncity;j++){
            if(i==j) continue;
            InfoFile *info = info_table[i][j]; // note table[j][i] orientation earlier? check: we stored table[j][i] with i as city_i index (from file). But for pair (i,j), we need info either at [i][j] or [j][i].
            int swapped=0;
            if(!info){ info=info_table[j][i]; swapped=1; if(!info) continue; }
            double sum_i_for_j=0,sum_j_for_i=0; int count_i_for_j=0,count_j_for_i=0;
            for(int k=0;k<info->n;k++){
                uint64_t key=info->positions[k];
                int pos_idx=find_position(gd->positions, gd->n_positions, key);
                if(pos_idx<0) continue;
                unsigned char gt = gd->genotypes[idx][pos_idx];
                if(gt!=info->values[k]) continue;
                double pi = swapped ? info->p_j[k] : info->p_i[k];
                double pj = swapped ? info->p_i[k] : info->p_j[k];
                sum_i_for_j += pj; // evidence for j
                sum_j_for_i += pi; // evidence for i
                if (pj < pi) count_j_for_i++; else count_i_for_j++; // if equal -> i default in description
            }
            if(sum_i_for_j>sum_j_for_i){
                sum_cols[j]+=1; sum_q_cols[j]+= swapped?Q[j][i]:Q[i][j];
            } else if(sum_j_for_i>sum_i_for_j){
                sum_cols[i]+=1; sum_q_cols[i]+= swapped?Q[i][j]:Q[j][i];
            }
            if(count_i_for_j>count_j_for_i){
                count_cols[j]+=1; count_q_cols[j]+= swapped?Q[j][i]:Q[i][j];
            } else if(count_j_for_i>count_i_for_j){
                count_cols[i]+=1; count_q_cols[i]+= swapped?Q[i][j]:Q[j][i];
            } else { // tie use sums
                if(sum_j_for_i>sum_i_for_j){
                    count_cols[i]+=1; count_q_cols[i]+= swapped?Q[i][j]:Q[j][i];
                } else {
                    count_cols[j]+=1; count_q_cols[j]+= swapped?Q[j][i]:Q[i][j];
                }
            }
        }
    }
    double total=0; for(int c=0;c<ncity;c++) total+=sum_cols[c]; if(total>0) *has_match=1; else *has_match=0;
    // store results
    for(int c=0;c<ncity;c++){ sum_res[c]=sum_cols[c]; count_res[c]=count_cols[c]; sum_q_res[c]=sum_q_cols[c]; count_q_res[c]=count_q_cols[c]; }
    // max check
    double max_val=-1; int correct=0; for(int c=0;c<ncity;c++){ if(sum_cols[c]>max_val) {max_val=sum_cols[c]; correct=(c==samples[idx].city_idx); } else if(sum_cols[c]==max_val && c==samples[idx].city_idx) correct=1; }
    *correct_sum=correct;
    max_val=-1; correct=0; for(int c=0;c<ncity;c++){ if(count_cols[c]>max_val){max_val=count_cols[c]; correct=(c==samples[idx].city_idx);} else if(count_cols[c]==max_val && c==samples[idx].city_idx) correct=1; }
    *correct_count=correct;
    max_val=-1; correct=0; for(int c=0;c<ncity;c++){ if(sum_q_cols[c]>max_val){max_val=sum_q_cols[c]; correct=(c==samples[idx].city_idx);} else if(sum_q_cols[c]==max_val && c==samples[idx].city_idx) correct=1; }
    *correct_sum_q=correct;
    max_val=-1; correct=0; for(int c=0;c<ncity;c++){ if(count_q_cols[c]>max_val){max_val=count_q_cols[c]; correct=(c==samples[idx].city_idx);} else if(count_q_cols[c]==max_val && c==samples[idx].city_idx) correct=1; }
    *correct_count_q=correct;
    free(sum_cols); free(count_cols); free(sum_q_cols); free(count_q_cols);
}

static void write_results(const char *prefix, Sample *samples, int n_samples, double **sum_res, double **count_res, double **sum_q_res, double **count_q_res) {
    char path[256];
    int ncity=num_cities;
    // helper to write one file
    const char *names[4] = {"sum","count","sum_q","count_q"};
    double ***res_arr = (double***)malloc(4*sizeof(double**));
    res_arr[0]=sum_res; res_arr[1]=count_res; res_arr[2]=sum_q_res; res_arr[3]=count_q_res;
    for(int f=0;f<4;f++){
        snprintf(path,sizeof(path),"%s_%s.csv",prefix,names[f]);
        FILE *fp=fopen(path,"w"); if(!fp){perror("write result"); continue;}
        fprintf(fp,"iid"); for(int c=0;c<ncity;c++) fprintf(fp,",%s",cities[c]); fprintf(fp,",max,label\n");
        for(int s=0;s<n_samples;s++){
            fprintf(fp,"%s",samples[s].iid);
            double max_val=-1; char max_buf[512]={0}; int first=1;
            for(int c=0;c<ncity;c++){
                double v=res_arr[f][s][c]; fprintf(fp,",%g",v);
                if(v>max_val){ max_val=v; snprintf(max_buf,sizeof(max_buf),"%s",cities[c]); first=0; }
                else if(v==max_val){ if(first){snprintf(max_buf,sizeof(max_buf),"%s",cities[c]); first=0;} else {strcat(max_buf,";"); strcat(max_buf,cities[c]);}}
            }
            fprintf(fp,",%s,%s\n",max_buf,samples[s].city);
        }
        fclose(fp);
    }
}

static void write_stats(const char *prefix, int n_samples, int *has_match, int *correct_sum, int *correct_count, int *correct_sum_q, int *correct_count_q) {
    char path[256]; snprintf(path,sizeof(path),"%s_stats.txt",prefix);
    int contributors=0, correct=0; for(int i=0;i<n_samples;i++){ if(has_match[i]) contributors++; if(correct_sum[i]) correct++; }
    int contributors_count=0, correct_cnt=0; for(int i=0;i<n_samples;i++){ if(has_match[i]) contributors_count++; if(correct_count[i]) correct_cnt++; }
    int contributors_sum_q=0, correct_sq=0; for(int i=0;i<n_samples;i++){ if(has_match[i]) contributors_sum_q++; if(correct_sum_q[i]) correct_sq++; }
    int contributors_count_q=0, correct_cq=0; for(int i=0;i<n_samples;i++){ if(has_match[i]) contributors_count_q++; if(correct_count_q[i]) correct_cq++; }
    FILE *fp=fopen(path,"w"); if(!fp){perror("write stats"); return;}
    fprintf(fp,"sum_contributors %d\n",contributors);
    fprintf(fp,"sum_correct %d\n",correct);
    fprintf(fp,"count_contributors %d\n",contributors_count);
    fprintf(fp,"count_correct %d\n",correct_cnt);
    fprintf(fp,"sum_q_contributors %d\n",contributors_sum_q);
    fprintf(fp,"sum_q_correct %d\n",correct_sq);
    fprintf(fp,"count_q_contributors %d\n",contributors_count_q);
    fprintf(fp,"count_q_correct %d\n",correct_cq);
    fclose(fp);
}

int main(int argc, char **argv) {
    if(argc<7){
        fprintf(stderr,"Usage: %s prefix genotype_vcf labels informative_dir fst_file threads\n",argv[0]);
        return 1;
    }
    const char *prefix=argv[1];
    const char *vcf_path=argv[2];
    const char *labels_path=argv[3];
    const char *info_dir=argv[4];
    const char *fst_path=argv[5];
    int nthreads=atoi(argv[6]);

    int n_samples=0;
    Sample *samples=load_labels(labels_path,&n_samples);
    if(!samples){fprintf(stderr,"Failed to load labels\n");return 1;}
    GenotypeData *gd=load_vcf(vcf_path,samples,n_samples);
    if(!gd){fprintf(stderr,"Failed to load vcf\n");return 1;}
    double **Q=load_fst_Q(fst_path);
    if(!Q){fprintf(stderr,"Failed to load fst\n");return 1;}
    InfoFile ***info_table=load_info_files(info_dir);

    // allocate result arrays
    double **sum_res=(double**)malloc(n_samples*sizeof(double*));
    double **count_res=(double**)malloc(n_samples*sizeof(double*));
    double **sum_q_res=(double**)malloc(n_samples*sizeof(double*));
    double **count_q_res=(double**)malloc(n_samples*sizeof(double*));
    int *has_match=(int*)calloc(n_samples,sizeof(int));
    int *correct_sum=(int*)calloc(n_samples,sizeof(int));
    int *correct_count=(int*)calloc(n_samples,sizeof(int));
    int *correct_sum_q=(int*)calloc(n_samples,sizeof(int));
    int *correct_count_q=(int*)calloc(n_samples,sizeof(int));
    for(int i=0;i<n_samples;i++){
        sum_res[i]=(double*)calloc(num_cities,sizeof(double));
        count_res[i]=(double*)calloc(num_cities,sizeof(double));
        sum_q_res[i]=(double*)calloc(num_cities,sizeof(double));
        count_q_res[i]=(double*)calloc(num_cities,sizeof(double));
    }

    #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
    for(int i=0;i<n_samples;i++){
        process_individual(i,samples,gd,info_table,Q,sum_res[i],count_res[i],sum_q_res[i],count_q_res[i],
                           &has_match[i],&correct_sum[i],&correct_count[i],&correct_sum_q[i],&correct_count_q[i]);
    }

    write_results(prefix,samples,n_samples,sum_res,count_res,sum_q_res,count_q_res);
    write_stats(prefix,n_samples,has_match,correct_sum,correct_count,correct_sum_q,correct_count_q);

    return 0;
}
