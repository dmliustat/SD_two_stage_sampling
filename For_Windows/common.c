#define _USE_MATH_DEFINES
#include "header.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


static double const normalization = 0.5 * M_SQRT1_2 * M_2_SQRTPI;

double *dou_vec_init(const size_t num);
long *long_vec_init(const size_t num);
short *sh_vec_init(const size_t num);
double runif(void);
double Normal_pdf(const double *empirical, double value);
void add_new(long *TYPE, long *type2node, long *ntype, long *newseq, long nnode, long nbranch);
void remove_old(long *TYPE, long *type2node, long *ntype, long c, long onode, long nbranch);
long add_s(long* TYPE, long* ntype, long* newseq, long nbranch);
int remove_s(long* TYPE, long* ntype, long c, long nbranch);
void reduce_seq(long *TYPE, long *ntype, long c);
long sample_node(long *TYPE, long *type2node, long k);
void shuffle(long* array, size_t n);


void seed_set(FILE *a);
void seed_put(FILE *a);

unsigned long seed[2];


extern long L;

double *dou_vec_init(const size_t num)
{
  double *p;
  if((p = (double *) calloc(num, sizeof(double))) == NULL)
    perror("initialise_array():calloc failed!");
  return(p);
}


long *long_vec_init(const size_t num)
{
  long *p;
  if((p = (long *) calloc(num, sizeof(long))) == NULL)
    perror("initialise_array():calloc failed!");
  return(p);
}


short *sh_vec_init(const size_t num)
{
  short *p;
  if((p = (short *) calloc(num, sizeof(short))) == NULL)
    perror("initialise_array():calloc failed!");
  return(p);
}

/*S+ runif routine*/
double runif(void)
{
	unsigned long n, lambda=69069;
	*seed = *seed * lambda;
	*(seed+1) ^= *(seed+1) >> 15;
	*(seed+1) ^= *(seed+1) << 17;
	n = *seed ^ *(seed+1);
	return( (((n>>1) & 017777777777)) / 2147483648. );
}


double Normal_pdf(const double *empirical, double value)
{
   return  normalization * exp( - 0.5 * (value-empirical[0]) * (value-empirical[0]) / (empirical[1] * empirical[1])) / empirical[1];
}


/*add new type to sample*/
void add_new(long *TYPE, long *type2node, long *ntype, long *newseq, long nnode, long nbranch)
{
	long i, j, c, flag, tot;
	
	c = -1;
	tot = 0;
	for(i=0; i<*ntype && c==-1; i++){
		flag = 1;
		for(j=0; j<L && flag==1; j++) flag = (*(TYPE+(L+1)*i+j) == *(newseq + j)); /*flag=1 the new type is identical to type i in H; flag=0 otherwise*/
		tot += *(TYPE + (L+1) * i + L);
		if (flag == 1) {
			c = i;
			break;
		}
	}
	if(c > -1){
		/*if the new type is identical to an existing type c, then add 1 to the number of this type*/
		*(TYPE + (L+1) * c + L) += 1;
		for (i = nbranch-1; i >= tot; i--) *(type2node + i + 1) = *(type2node + i);
		*(type2node + tot) = nnode;			
	} 
	else{ /*if the new type is not one of the existing types, add the new to H*/
		for(j=0; j<L; j++) *(TYPE + (L+1)* *ntype + j) = *(newseq + j);
		*(TYPE + (L+1)* *ntype + L) = 1;
		*ntype += 1;
		*(type2node + nbranch) = nnode;
	}
	return;
	
}


void remove_old(long *TYPE, long *type2node, long *ntype, long c, long onode, long nbranch)
{
	long i, j, tot;
	
	tot = 0;
	for(i=0; i<c; i++) tot += *(TYPE+(L+1)*i+L);

	if(*(TYPE + (L+1)*c + L) > 1){
		*(TYPE + (L+1)*c + L) -= 1;
		j = nbranch;
		for(i=tot; i<nbranch; i++){
			if(*(type2node+i) == onode) j = i;
			if(i>=j && i<nbranch-1) *(type2node+i) = *(type2node +(i+1));
		}		
		return;
	}
	else{
		for(i=c; i<(*ntype-1); i++){
			for(j=0; j<L+1; j++) *(TYPE + (L+1)*i + j) = *(TYPE + (L+1)*(i+1) + j);
		}
		j = nbranch;
		for(i=tot; i<nbranch; i++){
			if(*(type2node+i) == onode) j = i;
			if(i>=j && i<nbranch-1) *(type2node+i) = *(type2node +(i+1));
		}	

		*ntype -= 1;
		return;
	}
}

/*simplified version, without changing type2node. Return the type index of the new sequence*/
long add_s(long* TYPE, long* ntype, long* newseq, long nbranch)
{
	long i, j, c, flag;

	c = -1;
	for (i = 0; i < *ntype && c == -1; i++) {
		flag = 1;
		for (j = 0; j < L && flag == 1; j++) flag = (*(TYPE + (L + 1) * i + j) == *(newseq + j)); /*flag=1 the new type is identical to type i in H; flag=0 otherwise*/
		if (flag == 1) {
			c = i;
		}
	}
	if (c > -1) {
		/*if the new type is identical to an existing type c, then add 1 to the number of this type*/
		*(TYPE + (L + 1) * c + L) += 1;
	}
	else { /*if the new type is not one of the existing types, add the new to H*/
		for (j = 0; j < L; j++)* (TYPE + (L + 1) * *ntype + j) = *(newseq + j);
		*(TYPE + (L + 1) * *ntype + L) = 1;
		c = *ntype;
		*ntype += 1;
	}
	return c;

}

/*return 1 if ntype deduct 1; 0 otherwise*/
int remove_s(long* TYPE, long* ntype, long c, long nbranch)
{
	long i, j;

	if (*(TYPE + (L + 1) * c + L) > 1) {
		*(TYPE + (L + 1) * c + L) -= 1;
		return 0;
	}
	else {
		for (i = c; i < (*ntype - 1); i++) {
			for (j = 0; j < L + 1; j++)* (TYPE + (L + 1) * i + j) = *(TYPE + (L + 1) * (i + 1) + j);
		}

		*ntype -= 1;
		return 1;
	}
}

void reduce_seq(long* TYPE, long* ntype, long c) {
	long i;

	if (*(TYPE + (L + 1) * c + L) > 1) {
		*(TYPE + (L + 1) * c + L) -= 1;
		return;
	}
	else {
		for (i = 0; i < L+1; i++) *(TYPE + (L + 1) * c + i) = -2;
		*ntype -= 1;
		return;
	}
}


/*********************************
 *                               *
 * Routine to read seed for file *
 *                               *
 *********************************/

void seed_set(FILE *file)
{
  long i,temp;
  unsigned long max=4294967295;
  char ch;
 
  for(i=0;i<2;i++){
    /*ignore white space*/
    
    while(1){
      ch = fgetc(file);
      if(ch != '\n'&& ch != '\t') break;
    }
    
    temp= ch-'0';
    
    while(1){
      ch = fgetc(file);
      if(ch == '\n' || ch == '\t' || ch ==EOF) break;
      temp= (10*temp+ch-'0') % max;
    }
    *(seed+i)=temp;
  }
}


/***********************
 *                     *
 * put seed in file    *
 * file_r and file_w   *
 * point to same file  *
 * but for read (r)    *
 * and write (w)       *
 *                     *
 ***********************/

void seed_put(FILE *file_w)
{
 
  (void)fprintf(file_w,"%lu\n%lu", *seed,  *(seed+1));
}   


long sample_node(long *TYPE, long *type2node, long k)
{
	long i, num, startnum;

	num = *(TYPE+(L+1)*k+L);
	startnum = 0;

	if(k>0){
		for(i=0; i<k; i++){
			startnum += *(TYPE+(L+1)*i+L);
		}
	}

	i = rand() % num + startnum;
	
	return(*(type2node+i));

}


void shuffle(long* array, size_t n)
{
	if (n > 1) {
		size_t i;
		for (i = 0; i < n - 1; i++) {
			size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
			long t = array[j];
			array[j] = array[i];
			array[i] = t;
		}
	}
}

