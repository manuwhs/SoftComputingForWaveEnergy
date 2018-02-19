
#ifndef GAHEADER_H_
#define GAHEADER_H_

#define  ROULETTE_WHEEL_SEL	0
#define  RANK_SEL			1
#define  TOURNAMENT_SEL		2

			
typedef struct {
		int n_organisms;	// Number of total organisms of a generation	Sould be 2*num_bits min
		int n_generations;	// Number of generations left to do.
		int num_bits;		// Number of bits of the population

		char sel_type;		// Type of selection algorithm.
		int stop_cond;		// Stop condition
		
		int max_mut;		// Number maximum of mutations
		int n_Selected;		// Number of Selected organisms for the initial population
		int Min_fighters;		// Number maximum of fighters per tournament		
		int Max_fighters;		// Number maximum of fighters per tournament
										
		unsigned int **Population;	// Population of organisms [n_organisms][sizeofvector]
		unsigned int **Aux_Population;	//Auxiliar population to do the selection.
		unsigned int **Selected_Population;	//Selected organisms for the initial population
				
		double *Evaluation;		// Evaluation values of the organisms
		double *Fitness;		// Fitness values of the organisms
		int *order;				// Fitness list will be ordered. This will have that order

		
		float Elit_rate;	// Proportion of the population that will surely be selected to the next generation	[0-1]
		float Xover_rate;	// Probability of crossover   Should be around 85% - 60% [0-1]
		float Mut_rate;		// Probability of mutation %	  Should be around 1-5% [0-1]
		
		unsigned int *fights_won;	// For the tournamente selection algorithm

	}GA_data;

// ------------> GAfunctions
void print_GA_data (GA_data *GA_p);
int GA_op(GA_data *GA_p, int FLAGS);
int crossover_1p(unsigned int *a, unsigned int *b, int num_bits);
int get_Next_Gen(GA_data *GA_p);
int get_Evaluation(GA_data *GA_p);
int get_Fitness(GA_data *GA_p);
int select_organisms(GA_data *GA_p);
int crossover_population(GA_data *GA_p);
int mutate_gen(GA_data *GA_p);

int Roulette_Wheel_Sel(GA_data *GA_p);
int Rank_Sel(GA_data *GA_p);
int Tournament_Sel(GA_data *GA_p);


int Init_GA_Population(GA_data *GA_p);
int print_organisms(GA_data *GA_p);
int print_GA_values(GA_data *GA_p);
int print_GA_generation(GA_data *GA_p);
#endif /* ELMHEADER_H_ */
