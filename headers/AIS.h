
#ifndef AISHEADER_H_
#define AISHEADER_H_


// Artificial inmune system funcions:
/*  We have 2 populations:  --> Any member of a population is a feasible solution
 * 	--> Antibodies: Random cells generated from the Antigens 
 * 				 -> Cloning has mutation 
 * 				 -> Hueeeesoss --> Amiiiigo
 * 	--> Antigens: They are the best possible solutions from where Antibodies are generated.
 * 
 * The system goes as follows:
 * - Initialize the antigen generation --> Generate a lot of solutions and get the best ones:
 * 		--> The antigens should differ from each other in a minimum number of bits -> Hamming distance
 * - EACH ITERATION:
 * 		Clone all (Several times) the antigens in the antibody population  --> With several mutations
 * 		Get their fitness
 * 		If (some fitness in the antibodies > antigens) -> Insert antibody into antigens from the tail
 * 			(antigens are from best to worst)
 * 			If (any existing antigen has a Hamming distance lower than specified, it is substituted with it,
 * 				if not, its just add it.
 * 		Eliminate the worst organisms, and keep the other best ones but mutating them a lot the next generation.
 * 
 * 		( For every antibody (nice solution) --> We generate X antibodies 
 * - 
 * - 
 */
			
typedef struct {
		unsigned int num_bits;		// Number of bits of the population	
		unsigned int n_cycles;		// Number of cycles this shit has reproduced
		
		unsigned int Ini_B;		// Number of initial B-cells
		unsigned int Max_B;		// Number of maximum B-cells
		unsigned int n_B;		// Current number of B-cells
		
		unsigned int n_clones;	// Number of clones (antibodies) of every B-cell
		unsigned int n_Ab;		// Number of random Ab created

		unsigned int B_domain;		// Minimum Hamming distance between 2 B-cells
									// The mutation of B-cells will be this size maximum
		unsigned int B_lifetime;	// Number of cycles an B is alive without getting upgraded
									// If it expires, it becomes a Master cell
										
		unsigned int Max_M;		// Maximum number of Master cells
		unsigned int n_M_cells;		// Current number of Master cells

		int stop_cond;				// Stop condition
				
		unsigned int **M_cells;// Population B-cells that could not be upgraded any further.
		unsigned int **B_cells;		// Population of B-cells [Max_Ag]
		unsigned int **Ab_cells;	// Population of Ab cells [Max_Ag*n_clones + n_Ab]
		
		double *Ab_Defense;			// How good the Ab-cells are (ERMS)
		double *B_Defense;			// How good the B-cells are (ERMS)
		double *M_Defense;		// How good the M-cells are (ERMS)
		
		unsigned int *Beato_cell;	// Best M-cell of all
		double Beato_Defense;		// How good the best M-cell is.
					
		unsigned int *B_rem_lifes;	// Lifes left for the B-cells
		int *order;					// Aux variable to keep order the Ab pool
		
	}AIS_data;

// ------------> SA functions

int AIS_op(AIS_data *AIS_p, int FLAGS);
int clone_Bs(AIS_data *AIS_p);
int mutate_Abs(AIS_data *AIS_p);
int Abs_defense(AIS_data *AIS_p);
int remove_B(AIS_data *AIS_p, unsigned int Ag_pos);
int	upgrade_Bs(AIS_data *AIS_p);
int Init_AIS (AIS_data *AIS_p);
int print_AIS_values(AIS_data *AIS_p);
int print_AIS_data(AIS_data *AIS_p);
void print_AIS_stuff (AIS_data *AIS_p);
#endif /* ELMHEADER_H_ */
