#include <stdio.h>
#include <math.h>
#include "header.h"
#include "common.h"
#define ERROR_MIN (0.000000001)


extern long L;
extern long M;
extern double Mu;
extern double P[2][2];


/*calculate the proportion of haplotype 0 at each ancestral site */
void calcProp(double *prop, long *TYPE, long ntype, long *index, long index_len);



void calcProp(double *prop, long *TYPE, long ntype, long *index, long index_len)
{
  /*prop is a vector of length index_len */
  long i, j, p0, p1;

  for(i=0; i<index_len; i++)
  {
    p0 = 0;
    p1 = 0;
    for(j=0; j<ntype; j++)
    {
      if(*(TYPE+(L+1)*j+index[i]) == 0) p0 += *(TYPE+(L+1)*j+L);
      if(*(TYPE+(L+1)*j+index[i]) == 1) p1 += *(TYPE+(L+1)*j+L);
    }

    if(p0+p1 > 0) prop[i] = (double)p0/((double)p0+(double)p1);
    else prop[i] = 0.5;
  }
}


