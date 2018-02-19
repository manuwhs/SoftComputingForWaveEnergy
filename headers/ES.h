
#ifndef ESHEADER_H_
#define ESHEADER_H_

typedef struct {
	int num_bits;		// Maximum number of bits (Total parameters)
	int num_1s_min;		// Minimum number of bits of the selected vectors
	int num_1s_max;		// Maximum number of bits of the selected vectors
						// Not counting the Es_fixed_1s
	int n_TOP;			// Top solutions that will be printed out.
	unsigned int * ES_selected;  //	Selected vector with the 1s we will do.
	unsigned int * ES_fixed_1s;  //	Vector with tha 1s that will always stay in the ES
	
	int n_ES;						// Number of total solutions we will do.
	unsigned int ** ES_bitvectors;	// Solutions of the bitvectors
	double *ES_ERMS;
}ES_data;

// ------------------> ES functions <-----------
int get_ES_selected(unsigned int *ones_pos, unsigned int *v, unsigned int n, int num_bits);	// FOR ES
unsigned int ** get_ES_vectors (unsigned int *ES_fixed_1s, unsigned int *ES_selected, unsigned int min1s,unsigned int n,unsigned int *num_vectors);
int ES_op ( ES_data * ES_p );
int plot_ES_Search_Space (ES_data * ES_p);

#endif /* ELMHEADER_H_ */
	
