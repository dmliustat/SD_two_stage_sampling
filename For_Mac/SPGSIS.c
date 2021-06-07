#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tables.h"
#include "common.h"
#include "header.h"
#include "gamma.h"


/*Claim Functions defined in this file*/
void input1(FILE *infile);
void input2(FILE *infile);


/*Claim External Functions*/
/*extern void exp_matrix(double **Q_out, double para); */ /*Q_out=exp{P*para}*/
extern int Build(long *TYPE, long *type2node, long *ntype, const double RHO, const double THETA, double *logw, tsk_table_collection_t *tables, long *num_nodes, long *num_rec, long *num_mut, long *num_coal, long* num_nonrec, long* mut_by_site, long* mut_site, long* mut_state, long* mut_node);


extern void FDupdate(long* TYPE, long* type2node, const double RHO, const double THETA, long k, long* ntype, long* nbranch, long* n, double* m, double* logw, tsk_table_collection_t* tables, int* rec_num, int* mut_num, int* coal_num, int* nonrecc, long* mut_by_site, long* mut_site, long* mut_state, long* mut_node, long* eligible);

extern double rexp(double rate);



/*Constants from input*/
long NRUN;
long INNER;
long Ne; /*effective population size*/
double NumLineage; /*the number of lineages to stop */
double THETA; /*mutation rate*/
double RHO; /*recombination rate*/
double SEQLEN; /*sequence length*/

long M; /*number of genes*/
long L; /*number of loci*/
double Mu; /*stationary distn of allele 0 in the samples*/
long *TYPE_copy; /*store observed haplotypes*/
long ntype_copy;
double P[2][2]; /*mutation matrix*/
//double empirical[2]; /*empirical distribution for p(N_tau) */
double *positions; /*genetic position for each site */

/*other constants*/
double adj_rate;

/*use a structure to store mutations*/
struct mutations
{
    long site_id;
    long state;
    long node_id;
};
int cmp(const void* a, const void* b) {
    struct mutations* x = (struct str*)a;
    struct mutations* y = (struct str*)b;
    if ((*x).site_id < (*y).site_id) return 1;
    else return 0;
}


int main(int argc, char *argv[]){
    /*declare the vatiables here*/
    FILE *in1, *in2, *out, *seed_file, *seed_file_w;

    double logw, *wt, wt_sum=0.0, *logwt;
    long *TYPE, *type2node, *type2node_copy, ntype, *mut_by_site;
    long *n_nodes, *n_rec, *n_mut, *n_nonrec, *n_coal, num_rec, num_mut, num_coal, num_nodes, num_nonrec, mut;
    char filename[200];

    long i,j, l, k;
    int flag, ret;
    tsk_table_collection_t tables, tables_copy;
	//tsk_treeseq_t ts;
    long* mut_by_site_copy;
    long* mut_site, * mut_state, * mut_node;

    
    if(argc>1) in1 = fopen(argv[1],"r");
    else in1 = fopen("infile1","r");
    
    input1(in1);
    fclose(in1);
    /*read input1*/
    
    
    if(argc>2) in2 = fopen(argv[2],"r");
    else in2 = fopen("infile2","r");
    
    input2(in2);
    fclose(in2);
    /*read input2*/
    
	//positions[L] = SEQLEN;
	/*scaled mutation and recombination rates*/
	THETA = 4.0 * Ne * THETA * SEQLEN / (double)L;
	RHO = 4.0 * Ne * RHO;
	//adj_rate = 0.9061782;
	adj_rate = 1;


    TYPE = long_vec_init((L+1)*TYPE_MAX);  
    type2node = long_vec_init(NODE_MAX);
    type2node_copy = long_vec_init(NODE_MAX);
    wt = dou_vec_init(NRUN);
	n_nodes = long_vec_init(NRUN);
	n_rec = long_vec_init(NRUN);
	n_mut = long_vec_init(NRUN);
    n_nonrec = long_vec_init(NRUN);
	n_coal = long_vec_init(NRUN);
	logwt = dou_vec_init(NRUN);
    ret = tsk_table_collection_init(&tables, 0);
    ret += tsk_table_collection_init(&tables_copy, 0);
    if(ret<0) fprintf(stderr, "Initialise table collection failed.\n");
    mut_by_site_copy = long_vec_init(L);
    mut_by_site = long_vec_init(L);
    mut_site = long_vec_init(MUT_MAX);
    mut_state = long_vec_init(MUT_MAX);
    mut_node = long_vec_init(MUT_MAX);
    struct mutations mut_record[MUT_MAX];

   

    /* initialise sequence length and sites table */
    tables_copy.sequence_length = SEQLEN;
    for(i=0; i<L; i++){
        ret = tsk_site_table_add_row(&tables_copy.sites, positions[i], 0, 0, NULL, 0);
        if(ret<0) fprintf(stderr, "Initialise position %f failed.\n", positions[i]);
    }

    /* initialize sample nodes, flag=1 */
    for(i=0; i<M; i++){
        ret = tsk_node_table_add_row(&tables_copy.nodes, 1, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
        if(ret<0) fprintf(stderr, "Initialise sample %ld failed.\n", i);
        type2node_copy[i] = ret;
    }

    /* initialize mutation number on each site*/
    for (i = 0; i < L; i++) {
        mut_by_site_copy[i] = 0;
        for (j = 0; j < ntype_copy; j++) {
            if(*(TYPE_copy + (L + 1) * j + i) == 1) mut_by_site_copy[i] += *(TYPE_copy + (L + 1) * j + L);
        }
    }

    for(i=0; i<NRUN; i++){
		fprintf(stderr, "\n--------------------------------------------------------------");
		fprintf(stderr, "\n This is run %ld\n", i + 1);

		ntype = ntype_copy;
        ret = tsk_table_collection_copy(&tables_copy, &tables, TSK_NO_INIT);
        if(ret<0) fprintf(stderr, "Copy table collection failed.\n");

        for(j=0; j<((L+1)*ntype); j++) *(TYPE+j) = *(TYPE_copy+j);
        for(j=0; j<M; j++) *(type2node+j) = *(type2node_copy+j);
        for(j=M; j<NODE_MAX; j++) *(type2node+j) = -1;
        for (j = 0; j < L; j++) *(mut_by_site + j) = *(mut_by_site_copy + j);

        flag = Build(TYPE, type2node, &ntype, RHO, THETA, &logw, &tables, &num_nodes, &num_rec, &num_mut, &num_coal, &num_nonrec, mut_by_site, mut_site, mut_state, mut_node);
		
        if(flag==0){
            /*sort mutation site id and record mutations*/
            for (j = 0; j < num_mut; j++) {
                mut_record[j].site_id = mut_site[j];
                mut_record[j].state = mut_state[j];
                mut_record[j].node_id = mut_node[j];
            }
            qsort(mut_record, num_mut, sizeof(mut_record[0]), cmp);
            for (j = 0; j < num_mut; j++) {
                mut = tsk_mutation_table_add_row(&tables.mutations, mut_record[j].site_id, mut_record[j].node_id, -1, mut_record[j].state, 0, NULL, 0);
                if (mut < 0) fprintf(stderr, "Adding mutation failed.\n");
            }

            /*.trees output */
			/*.trees are efficient and for analyzing; .txt are inefficient and for viewing only*/			
            if(argc>3){
                snprintf(filename, 200, "%s_%ld.trees", argv[3], i+1);
                ret = tsk_table_collection_dump(&tables, filename, 0);
				if(ret<0) fprintf(stderr, "Dumping .trees failed!\n");
            } 
            else{
                snprintf(filename, 200, "C:\\Users\\blabl\\Dropbox\\output\\out_%ld.trees", i+1);
                ret = tsk_table_collection_dump(&tables, filename, 0);
				if (ret < 0) fprintf(stderr, "Dumping .trees failed!\n");
            }

			if (logw < -312) {
				fprintf(stderr, "Unnormalized weight too small. logw = %lf.\n", logw);
			}
			else {
				fprintf(stderr, "Unnormalized weight logw = %lf.\n", logw);
			}
			logwt[i] = logw;
            wt[i] = exp(logw);
			n_nodes[i] = num_nodes;
			n_rec[i] = num_rec;
			n_mut[i] = num_mut;
            n_nonrec[i] = num_nonrec;
			n_coal[i] = num_coal;
			
			fprintf(stderr, "Building ARG succeeds.\n");

			//snprintf(filename, 200, "out_%d.txt", i+1);
			//test = fopen(filename, "w");

			//tsk_table_collection_print_state(&tables, test);
			//fclose(test);

        } 
        else{
            printf("Number of distinct haplotypes too large, ignoring trees\n");
            i--;
        }		


        ret = tsk_table_collection_clear(&tables);
        if(ret<0) fprintf(stderr, "Clear tree collection failed.\n");

        if(i%1000 ==999) fprintf(stderr, "End of run %d \n", (int)i+1);

    }

	
    /*calc the normalized SIS weights*/
	for (i = 0; i < NRUN; i++) {
		wt_sum += wt[i];
	}
	for (i = 0; i < NRUN; i++) {
		wt[i] /= wt_sum;
	}
	

    /*output the weights*/
    if(argc>3){
        snprintf(filename, 200, "%s_wt.txt", argv[3]);
        out = fopen(filename, "w");
    } 
    else{
        snprintf(filename, 200, "out.txt");
        out = fopen(filename, "w");
    }
    if(out){
		fprintf(out, "run \t normalized weight \t log(weight) \t num_nodes \t num_recombinations \t num_mutations \t num_nonrecurrent \t num_coalescence \n");
		for (i = 0; i < NRUN; i++) {
			fprintf(out, "%ld \t %lf \t %lf \t %ld \t %ld \t %ld \t %ld  \t %ld\n", i + 1, wt[i], logwt[i], n_nodes[i], n_rec[i], n_mut[i], n_nonrec[i], n_coal[i]);
		}
		
		fclose(out);
    }
	else {
		for (i = 0; i < NRUN; i++) fprintf(stderr, "Cannot open the weight file.");
	}
    

    
    /*change .seed file*/
    seed_file_w = fopen(".seed", "w");
    seed_put(seed_file_w);
    fclose(seed_file_w);

    
    /*free the variables*/
    ret = tsk_table_collection_free(&tables);
    ret = tsk_table_collection_free(&tables_copy);
    free(TYPE);
    free(TYPE_copy);
    free(type2node);
    free(type2node_copy);
    free(wt);
    free(logwt);
    free(positions);

    free(n_nodes);
    free(n_rec);
    free(n_mut);
    free(n_nonrec);
    free(n_coal);
    free(mut_by_site_copy);
    free(mut_by_site);
    

    return(0);

}

void input1(FILE* infile)
{
    char string[200];
    char iden[2];
    short lines;
    short flag = 1;
    short i;

    NRUN = -1; /*number of runs*/
    INNER = -1; /*number of inner SIS iteratons*/
    Ne = -1;
    RHO = -1; /*recombination rate*/
    THETA = -1; /*mutation rate*/
    SEQLEN = -1;
    NumLineage = -1;
    lines = 0;

    fprintf(stderr, "\n Data in file 1: \n");

    while (flag) {
        if (fgets(string, sizeof(string), infile) == NULL) {
            fprintf(stderr, "Error in file \n");
            abort();
        }
        else {
            i = 0;
            while (string[i] == ' ' || string[i] == '\t') i++;
            iden[0] = string[i];
            switch (iden[0]) {
            case EOF:
                flag = 0;
                break;
            case '#':
                flag = 0;
                break;
            case '\n':
                break;

            case 'N':
                if (NRUN != -1) {
                    fprintf(stderr, "number of runs redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf(string + i, "%ld", &NRUN);
                    fprintf(stderr, "number of runs %ld \n", NRUN);
                    break;
                }

            case 'I':
                if (INNER != -1) {
                    fprintf(stderr, "number of inner SIS iteratoins redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf(string + i, "%ld", &INNER);
                    fprintf(stderr, "number of inner SIS iterations %ld \n", INNER);
                    break;
                }

            case 'E':
                if (Ne != -1) {
                    fprintf(stderr, "Effective population size redefined on line %d \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf(string + i, "%ld", &Ne);
                    fprintf(stderr, "effective population size %ld \n", Ne);
                    break;
                }


            case'R':
                if (RHO != -1) {
                    fprintf(stderr, "Recombination rate redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf(string + i, "%lf", &RHO);
                    fprintf(stderr, "recombination rate %3.9lf per generation per unit physical distance\n", RHO);
                    break;
                }

            case'M':
                if (THETA != -1) {
                    fprintf(stderr, "Mutation rate redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf(string + i, "%lf", &THETA);
                    fprintf(stderr, "mutation rate %3.9lf per generation per base pair\n", THETA);
                    break;
                }

            case'L':
                if (SEQLEN != -1) {
                    fprintf(stderr, "Sequence length redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf(string + i, "%lf", &SEQLEN);
                    fprintf(stderr, "sequence length %lf \n", SEQLEN);
                    break;
                }

            case'S':
                if (NumLineage != -1) {
                    fprintf(stderr, "Number of lineages redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf(string + i, "%lf", &NumLineage);
                    fprintf(stderr, "number of lineages %lf \n", NumLineage);
                    break;
                }

            case'F':
                fprintf(stderr, "Random seed ");
                while (string[i] != '=') i++;
                i++;

                /*for (j = 0; j < 2; j++) {
                    if (fscanf(infile, "%ld", &seed[j]) != 1) {
                        fprintf(stderr, "Error in random seed %d \n", j + 1);
                        abort();
                    }
                    else fprintf(stderr, "%ld", seed[j]);
                }
                fprintf(stderr, "\n");*/
                sscanf(string + i, "%ld", &seed[0]);
                fprintf(stderr, "%ld ", seed[0]);
                while (string[i] != ',') i++;
                i++;
                sscanf(string + i, "%ld", &seed[1]);
                fprintf(stderr, "%ld\n", seed[1]);
                break;


            default:
                fprintf(stderr, "Non-standard input line %d in file 1. Ignored", lines);
                break;
            }
            lines++;
        }
    }

    flag = 1;

    if (NRUN <= 0) {
        fprintf(stderr, "Number of runs <= 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Number of runs = %ld checked\n", NRUN);

    if (INNER <= 0) {
        fprintf(stderr, "Number of inner SIS iterations <= 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Number of inner SIS iterations = %ld checked\n", INNER);

    if (Ne <= 0) {
        fprintf(stderr, "Effective population size <= 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Effective population size = %ld checked\n", Ne);

    if (THETA < 0) {
        fprintf(stderr, "Theta value < 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Theta value = %3.9lf checked\n", THETA);

    if (RHO < 0) {
        fprintf(stderr, "Rho value < 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Rho value = %3.9lf checked\n", RHO);

    if (SEQLEN < 0) {
        fprintf(stderr, "Sequence length < 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Sequence length = %lf checked\n", SEQLEN);

    if (NumLineage < 0) {
        fprintf(stderr, "Number of lineages < 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Number of lineages = %lf to stop\n", NumLineage);
    if (flag == 0) abort();
}

void input2(FILE* infile)
{
    char string[200];
    char iden[2];
    short lines, flag, count;
    long i, j;
    long temp;
    double tot;

    L = 0; /*number of loci*/
    M = 0; /*number of genes*/
    Mu = 0;
    flag = 1;
    lines = 0;
    count = 0;
    ntype_copy = 0;

    fprintf(stderr, "\n Data in file 2: \n");

    while (flag) {
        if (fgets(string, sizeof(string), infile) == NULL) {
            fprintf(stderr, "Error in file \n");
            abort();
        }
        else {
            i = 0;
            while (string[i] == ' ' || string[i] == '\t') i++;
            iden[0] = string[i];
            switch (iden[0]) {
            case EOF:
                flag = 0;
                break;
            case '#':
                flag = 0;
                break;
            case '\n':
                break;

            case 'L':
                if (L > 0) {
                    fprintf(stderr, "number of loci redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf(string + i, "%ld", &L);
                    fprintf(stderr, "number of loci %ld \n", L);
                    break;
                }

            case 'G':
                if (M > 0) {
                    fprintf(stderr, "number of genes redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf(string + i, "%ld", &M);
                    fprintf(stderr, "number of genes %ld \n", M);
                    break;
                }

            case 'D':
                if (ntype_copy > 0) {
                    fprintf(stderr, "number of distinct haplotypes on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf(string + i, "%ld", &ntype_copy);
                    fprintf(stderr, "number of distinct haplotypes %ld \n", ntype_copy);
                    break;
                }

            case 'S':
                if (Mu > 0) {
                    fprintf(stderr, "Stationary distn redefined on line %d - ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf(string + i, "%lf", &Mu);
                    fprintf(stderr, "Stationary distn of allele 0 is %lf \n", Mu);
                    break;
                }

            case 'H':
                if (L == 0 || ntype_copy == 0) {
                    fprintf(stderr, "Haplotypes redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    if ((count & 0x1) == 0x1) {
                        fprintf(stderr, "Haplotypes redefined on line %d. Ignored \n", lines);
                        break;
                    }
                    TYPE_copy = long_vec_init((L + 1) * ntype_copy);
                    fprintf(stderr, "Haplotypes:\n");
                    for (j = 0; j < (L + 1) * ntype_copy; j++) {
                        if (fscanf(infile, "%ld", &TYPE_copy[j]) != 1) {
                            fprintf(stderr, "Error in haplotypes for %ld th entry \n", j + 1);
                            abort();
                        }
                        else {
                            if (j % (L + 1) != L) fprintf(stderr, "%ld ", TYPE_copy[j]);
                            else fprintf(stderr, "%ld\n", TYPE_copy[j]);
                        }
                    }
                    fprintf(stderr, "\n");
                    count += 1;
                    break;
                }

            case 'M':
                if ((count & 0x2) == 0x2) {
                    fprintf(stderr, "Mutation matrix redefined on line %d. Ignore \n", lines);
                    break;
                }
                fprintf(stderr, "Mutation matrix: ");
                for (i = 0; i < K; i++) {
                    for (j = 0; j < K; j++) {
                        if (fscanf(infile, "%lf", &P[i][j]) != 1) {
                            fprintf(stderr, "Error in inputing mutation matrix element (%ld,%ld) \n", i + 1, j + 1);
                            abort();
                        }
                        else fprintf(stderr, "%lf ", P[i][j]);
                    }
                }
                fprintf(stderr, "\n");
                count += 2;
                break;

                /*case 'E':
                    if((count&0x4)==0x4){
                        fprintf(stderr,"Empirical distn redefined on line %d - ignored \n",lines);
                        break;
                    }
                    fprintf(stderr, "Empirical distn for H_tau:");
                    for(i=0;i<2;i++){
                        if(fscanf(infile,"%lf",&empirical[i])!=1){
                            fprintf(stderr,"Error in empirical ditn for %ld th entry \n",i+1);
                            abort();
                        }
                        else fprintf(stderr, "%lf ", empirical[i]);
                    }
                    fprintf(stderr, "\n");
                    count+=4;
                    break;*/

            case 'P':
                if (L == 0) {
                    fprintf(stderr, "Genetic positions defined before number of loci on line %d - ignore \n", lines);
                    break;
                }
                if ((count & 0x8) == 0x8) {
                    fprintf(stderr, "Genetic positions redefined on line %d \n", lines);
                    break;
                }
                positions = dou_vec_init(L);
                for (i = 0; i < L; i++) {
                    if (fscanf(infile, "%lf", &positions[i]) != 1) {
                        fprintf(stderr, "Error in inputting genetic positions on locus %ld \n", i + 1);
                        abort();
                    }
                }
                count += 8;
                break;

            default:
                fprintf(stderr, "Non-standard input line %d in file 2. Ignored \n", lines);
                break;
            }
            lines++;
        }
    }
    flag = 1;
    if (L <= 1) {
        fprintf(stderr, "Number of loci <= 1 or undefined \n");
        flag = 0;
    }
    if (ntype_copy <= 0) {
        fprintf(stderr, "Number of haplotypes <= 0 or undefined \n");
        flag = 0;
    }
    if (M <= 0) {
        fprintf(stderr, "Number of genes <= 0 or undefined \n");
        flag = 0;
    }
    if (Mu <= 0) {
        fprintf(stderr, "Stationary distn <= 0 or undefined \n");
        flag = 0;
    }
    if (flag == 0) abort();

    /*check empirical distn */
    /*if((count&0x4)!=0x4){
        fprintf(stderr, "Empirical distn for N_tau not defined \n");
        flag=0;
    }
    else{
        for(i=0;i<2;i++){
            if(empirical[i]<=0){
                fprintf(stderr, "Mean and std of empirical distn should be positive \n");
                flag=0;
            }
        }
    }*/
    /*check positions */
    if ((count & 0x8) != 0x8) {
        fprintf(stderr, "Genetic positions not defined \n");
        flag = 0;
    }
    else {
        for (i = 0; i < L; i++) {
            if (positions[i] < 0) {
                fprintf(stderr, "Genetic position %ld should be non-negative \n", i + 1);
                flag = 0;
            }
        }
    }

    if (flag == 0) abort();

    /*check TYPE*/
    temp = 0;
    for (i = 0; i < ntype_copy; i++) {
        for (j = 0; j < (L); j++) {
            if (TYPE_copy[(L + 1) * i + j]<0 || TYPE_copy[(L + 1) * i + j]>K - 1) {
                fprintf(stderr, "Type for haplotype %ld at site %ld is not an allele \n", j + 1, i + 1);
                flag = 0;
            }
        }
        if (TYPE_copy[(L + 1) * i + L] <= 0) {
            fprintf(stderr, "Number of %ld th haplotype is less than or equal to 0 \n", i);
            flag = 0;
        }
        temp += TYPE_copy[(L + 1) * i + L];
    }
    if (temp != M) {
        fprintf(stderr, "Number of genes from haplotypesdoes not equal the number of genes specified \n");
        flag = 0;
    }

    /*check mutation matrix*/
    if ((count & 0x2) != 0x2) {
        fprintf(stderr, "Mutation matrix unspecified \n");
        flag = 0;
    }
    else {
        for (i = 0; i < K; i++) {
            tot = 0.0;
            for (j = 0; j < K; j++) {
                if (P[i][j] < 0 || P[i][j]>1) {
                    fprintf(stderr, "(%ld,%ld)th entry of mutation matrix not a probability\n", i + 1, j + 1);
                    flag = 0;
                }
                tot += P[i][j];
            }
            if (tot >= 1.0000001 || tot <= 0.999999) fprintf(stderr, "%ldth row of mutation matrix does not sum to one (renormalising)\n", i + 1);
            for (j = 0; j < K; j++) P[i][j] /= tot;
        }
    }


    if (flag == 1) fprintf(stderr, "All checks done\n");
    if (flag == 0) abort();
}












