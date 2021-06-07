# SamplingPartialGenealogies

The name of the software, SPGSIS, is short for ``sampling partial genealogies using sequential importance sampling". The SPGSIS.c source file consists of the main function and the input/output functions. Our software requires two input files and generates two kinds of output files, i.e., a bunch of .trees files that store the partial ARGs and a text file with the SIS weights for each partial ARG.


The LSprob.c consists of functions generating genetic events, updating partial ARGs and computing SIS weights. The three major functions are explained below.

- Build function: calls the FDupdate function to reconstruct the ARG sequentially, from present back to the time when the stopping criterion is met, and then computes the un-normalised SIS weights.

- FDupdate function: implements Fearnhead and Donnelly's sampling algorithm: a) choosing a chromosome out of the current configuration and b) sampling a genetic event to occur to the chosen chromosome.

- CSDapprx function: approximates the conditional sampling density using Li and Stephens' approach with variations.


The gamma.c consists of functions to calculate some required quantities used in the CSDapprx function. The common.c includes some small functions that are frequently called in other .c files, such as allocating memory for a vector, setting random seed, adding new and removing old haplotypes.



In input 1, the user should provide the following information:

- The number of runs;

- The number of inner SIS iterations;

- The effective population size;

- The recombination rate per base pair per unit physical distance; 

- The mutation rate per base pair per generation;

- The sequence length;

- The stopping criterion (i.e., the wanted number of lineages).


In input 2, the user should provide:

- The number of sites and the genetic position for each site;

- The number of haplotypes at present;

- The distinct haplotypes and their multiplicities;

- The stationary distribution of allele 0; If unknown, we use the proportion of allele 0 in the samples as an estimate of its stationary distribution;

- The mutation transition matrix in the form of a vector; The first two values denote the transition rate of allele 0 being mutated from allele 0 and 1, respectively; And the last two values denote the transition rate of allele 1 being mutated from allele 0 and 1, respectively. 



We use tskit, a C API, in our program to store the ARGs and write the output in the format of tree sequence. We use the Tables API only and thus the files tskit/core.[c,h] and tskit/tables.[c,h] are needed. The necessary files for using the Tables API, i.e., core.[c,h], tables.[c,h] and kastore.[c,h], are also uploaded. For details and source code of the tskit API, see [this link](https://tskit.readthedocs.io/en/latest/c-api.html#).


Our software generates one .trees file for each sampled partial ARG. Their corresponding SIS weights are stored in a text file. 

The user can use the following commands for compilation:

gcc -c common.c

gcc -c gamma.c

gcc -c LSprob.c

gcc -c tables.c

gcc -c kastore.c

gcc -c core.c

ar -rc liballfiles.a common.o gamma.o LSprob.o tables.o kastore.o core.o

gcc SPGSIS.c -L. -lallfiles -o SPGSIS

.\SPGSIS.exe infile1 infile2 out 