
#ifndef SAHEADER_H_
#define SAHEADER_H_

#define  AUTO_TUNNING		0
#define  SET_TUNNING		1

// Simmulated anneling functions
typedef struct {
	unsigned int num_bits;		// Maximum number of bits (Total parameters)
	unsigned int n_ingots;		// Number of parallel tries.
	unsigned int n_freezings;	// Number of freezings done so far.
	int stop_cond;				// Stop condition
		
	double temp;		// Simulated anneling temperature
	double k;			// Simulated constant  temp[t+1] = temp[t]*k

	int max_mut;			// Number of maximum mutations per freezing
	unsigned int n_Selected;// Number of selected ingots for initial population
			
	unsigned int **Ingots;		// Ingots for the SA [n_organisms][sizeofvector]
	unsigned int **New_Ingots;	// New Ingots for the SA [n_organisms][sizeofvector]
	
	double *Energy;				// Energy of the current ingots
	double *New_Energy;			// Energy of the next freezing ingots
	
	unsigned int **Top_Ingots;	// Best ingots so far
	double *Top_Energy;			// Energy of the ingots
	
	unsigned int **Selected_Ingots;	// Ingots for the SA [n_organisms][sizeofvector]
	
}SA_data;

// ------------> SA functions

int SA_op(SA_data *SA_p, int FLAGS);
int get_Energy(SA_data *SA_p);
int get_New_Energy(SA_data *SA_p);
int freeze_ingots (SA_data *SA_p);
int apply_SA(SA_data *SA_p);
int auto_tune_k(SA_data *SA_p);
int Init_SA_Ingots(SA_data *SA_p);
int check_SA_convergence (SA_data *SA_p, double *ener, double *best_ener );

int print_SA_values(SA_data *SA_p);
int print_Ingots(SA_data *SA_p);
int print_SA_data(SA_data *SA_p);
void print_SA_stuff (SA_data *SA_p);
#endif /* ELMHEADER_H_ */
