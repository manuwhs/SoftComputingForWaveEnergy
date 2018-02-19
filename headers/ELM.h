
#ifndef ELMXHEADER_H_
#define ELMXHEADER_H_
#include "gsl_headers.h" 

//-------------> CONSTANTS <--------------------
#define NUM_CONFIG		20
#define SHOW_TIME_F   	1 << 0
#define CHECK_HINV_F	1 << 1
#define PLOT_DATA_F		1 << 2

//------------> NEW DATA TYPES <----------------

typedef struct {
	unsigned int n;				// Number of input parameters
	unsigned int n_chosen;		// Number of chosen parameters.
	unsigned int *chosen_vector;// Chosen vector from the parameters
	unsigned int n_ELMs;		// Number of times we will do the ELM to average
	unsigned int Nh;			// Number of hidden neurons of the ELM
	unsigned int activation_f; 	// Activation function used
	unsigned int FLAGS; 		// ELM FLAGS for differente pourposes
	
	gsl_vector ** x_train_param;// Array of gsl_vectors of parameters [n][Ntrain]
	gsl_vector ** x_test_param;
	gsl_vector *y_train;		// gsl_vector with the training outputs. [Ntrain]
	gsl_vector *y_test;
	
	int t_max;					// Normalization variables of the output.
	int t_min;
	double ERMS;				// RMSE error of the ELM
	
}ELM_params;

typedef struct {
	char train_input_dir[255];		// File with the training input vectors
	char train_output_dir[255];		// File with the training results
	char test_input_dir[255];		// File with the testing input vectors
	char test_output_dir[255];		// File with the testing results

	int Ntrain;		// Number of training vectors we will read
	int Ntest;		// Number of training vectors we will read
	
	int n_input;	// Number of parameters of the vectors
	int n_output;	// Number of output parameters (Always 1 so far)
	
	unsigned int n_ELMs;		// Number of times we will do the ELM to average
	int Nh;			// Number of hidden neurons
	unsigned int activation_f; 	// Activation function used
	int n_threads;
	}Input_params;

typedef struct {						// Number of parameters loaded
	unsigned int n_ELMs;
	unsigned int n_Selected;
	unsigned int **Selected;	// Chosen vector from the parameters
	unsigned int FLAGS;
	char get_distrib;		// Shows the random distribution of the ELM
	
	unsigned int n_ELMs_dis;	// Chosen vector from the parameters
	unsigned int n_divs;
	}ELMs_data;	
		
// ------------------> ELM functions <-----------
void set_WeightBias(gsl_matrix **W, gsl_vector **b, int Nh, int n);
gsl_matrix* get_H_matrix(gsl_matrix *Xtrain, int Nh,gsl_matrix *W,gsl_vector *b, unsigned int act_func);
gsl_matrix* get_Hinv_matrix(gsl_matrix *H);
gsl_vector* get_beta_vector(gsl_matrix *Hinv, gsl_vector *ytrain);
gsl_vector* test_ELM(gsl_matrix *Htest, gsl_vector *beta);
double check_psudoinverse(gsl_matrix *H, gsl_matrix *Hinv);
gsl_matrix* get_Hinv_matrix2(gsl_matrix *H);
double get_ELM_RMSE(ELM_params *p);
int ELMs_op (ELMs_data *ELM_p);

int get_ELM_random_distrib(unsigned int * selected, int n_ELMs, int n_divs);
int print_ELMs_values(ELMs_data *ELMs_p);
int print_ELM_values(Input_params *ELM_p);

#endif /* ELMHEADER_H_ */
