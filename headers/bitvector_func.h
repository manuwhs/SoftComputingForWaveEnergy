#ifndef BITVECTORHEADER_H_
#define BITVECTORhEADER_H_

#define MAX_BITS 1000
#define MAX_CONF_SIZE 100
// -------------------> Bitvector Functions <----------------
int get_vector_weight(unsigned int *v, int length);
int get_random_bitarray(unsigned int *v, int max_bits);
int print_bitvector(unsigned int *v, int max_length);
int get_bitvector(unsigned int *v, int max_length, char *bitarray);
int or_bitvector(unsigned int *d, unsigned int *s, int num_bits);
int and_bitvector(unsigned int *d, unsigned int *s, int num_bits);
int bitvector_cpy(unsigned int *d, unsigned int *s, int num_bits);
unsigned int * get_1s_postions(unsigned int * vector, unsigned int n);
int zero_bitvector(unsigned int *v, int num_bits);
int fill_bitvector(unsigned int *v, int num_bits);
int mutate_bitvector(unsigned int *v, int n_bits, int n_mutations);
int bitvector_cmp(unsigned int *v, unsigned int *w);
unsigned int hamming_distance (unsigned int *v, unsigned int *w,int max_length);
unsigned int bitvector_setbit(unsigned int *v, int bit);
unsigned int bitvector_clearbit(unsigned int *v, int bit);
int get_ubits(unsigned int v, int max_length, char *bitarray);
unsigned int bin_to_gray (unsigned int v);

// Gets unsigned int bitvector from String of 0 and 1s
int get_u_vector(unsigned int *v, int max_length, char *bitarray);
#endif /* ELMHEADER_H_ */
