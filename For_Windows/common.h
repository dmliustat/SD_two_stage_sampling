extern unsigned long seed[2];

extern double *dou_vec_init(const size_t num);

extern long *long_vec_init(const size_t num);

extern short *sh_vec_init(const size_t num);

extern double runif(void);

extern double Normal_pdf(const double *empirical, double value);

extern void seed_set(FILE *a);

extern void seed_put(FILE *a);

extern void add_new(long* TYPE, long* type2node, long* ntype, long* newseq, long nnode, long nbranch);

extern void remove_old(long* TYPE, long* type2node, long* ntype, long c, long onode, long nbranch);

long add_s(long* TYPE, long* ntype, long* newseq, long nbranch);

int remove_s(long* TYPE, long* ntype, long c, long nbranch);

extern void reduce_seq(long* TYPE, long* ntype, long c);

extern long sample_node(long* TYPE, long* type2node, long k);

void shuffle(long* array, size_t n);
