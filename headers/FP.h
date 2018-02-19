
#ifndef FPHEADER_H_
#define FPHEADER_H_


#define  FS_Fitness		1 << 0
#define  FS_Affinity	1 << 1
#define  SFS			1 << 2
#define  SBE			1 << 3

typedef struct {
	int n_F;				// Number of Features
	unsigned int Operations;	// Operation to do
	double *F_Fitness;		// Fitness values of every parameter
	double ** F_Pairs;		// Fitness of every 2 pairs of parameters
	double ** F_Affinity;	// Affinity between every 2 pairs of parameters
	
	
	int * order;	
	unsigned int **F_Selected; // Array with the Selected vectors to do ELM (n max)
	unsigned int **Aux_Selected; // Array with the Selected vectors to do ELM (n max)
	double *Aux_Fitness;		// Aux variable for fitness
}FP_data;

// ------------------> ES functions <-----------
int FP_op ( FP_data * FP_p );
int plot_FP_Fitness (FP_data * FP_p);
int plot_FP_Affinity(FP_data * FP_p);
int print_input_statistics (unsigned int ** input_vectors, int n, int Num);
int print_Fitness(FP_data * FP_p);
int print_Affinity(FP_data * FP_p, int num);
int print_vectors_data(unsigned int ** SV, int num_Vectors, char *file_name);

int get_possible_inclusions(unsigned int * VS, int n, unsigned int ** sol);
int get_possible_exclusions(unsigned int * VS, int n, unsigned int ** sol);

int FS_Evaluation_op(FP_data *FP_p);
int FS_Affinity_op(FP_data *FP_p);
int SFS_op(FP_data *FP_p);
int SBE_op(FP_data *FP_p);
#endif /* ELMHEADER_H_ */
	
