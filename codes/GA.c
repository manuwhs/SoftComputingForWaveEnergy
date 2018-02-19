#include "../ELMheader.h"

int GA_op(GA_data *GA_p, int FLAGS){
	int i,j;
	
	double *	Gen_best_ERMS;	// Array with the ERMS of the best organism of a generation
	double *	Gen_ave_ERMS;	// Array with the average ERMS  of a generation
	double *	Gen_worst_ERMS;	// Array with the ERMS of the worst organism of a generation
	double *	Gen_biodiversity;	// Array with the biodiversity of a generation
	multiple_graph mgra;
	gsl_vector **graph_GA;

	
	if (FLAGS & PLOT_DATA){
			//--------> Reserve memmory for the generation statistics <------------
		Gen_best_ERMS = (double *) malloc (sizeof(double)*GA_p->stop_cond);
		Gen_ave_ERMS = (double *) malloc (sizeof(double)*GA_p->stop_cond);	
		Gen_worst_ERMS = (double *) malloc (sizeof(double)*GA_p->stop_cond);
		Gen_biodiversity = (double *) malloc (sizeof(double)*GA_p->stop_cond);
		
		graph_GA = (gsl_vector **)malloc(4*sizeof(gsl_vector *));
		for (i = 0; i < 4; i++){
			graph_GA[i] = gsl_vector_alloc(GA_p->stop_cond);
		}
		for (i = 0; i < GA_p->stop_cond; i++){
			graph_GA[0]->data[i] = i;
		}	
	}	
	//--------> Reserve memmory for the operations <------------
	GA_p->Population = (unsigned int **) malloc (sizeof(unsigned int *)*GA_p->n_organisms);
	GA_p->Aux_Population = (unsigned int **) malloc (sizeof(unsigned int *)*GA_p->n_organisms);

	GA_p->Evaluation = (double *) malloc (sizeof(double)*GA_p->n_organisms);
	GA_p->Fitness= (double *) malloc (sizeof(double)*GA_p->n_organisms);
	GA_p->order = (int *) malloc (sizeof(int)*GA_p->n_organisms);
	
	GA_p->fights_won = (unsigned int *) malloc (GA_p->n_organisms*sizeof(unsigned int)); // Do it somewhere else
	
	for (i = 0; i < GA_p->n_organisms; i++){
		GA_p->Population[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
		GA_p->Aux_Population[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
	}
	
	//--------> Generate initial population <------------
	Init_GA_Population(GA_p);
	print_organisms(GA_p);
	print_GA_values(GA_p);
//--------> Evaluate Generations <------------	
	for (i = 0; i < GA_p->stop_cond; i++){

		printf("Generation %i \n", GA_p->n_generations );		
	
		get_Evaluation(GA_p);			// GET THE EVALUATION
		
		print_GA_data (GA_p);				// PRINT GENERATION		
		
		get_Next_Gen(GA_p);				// GETS THE NEW GENERATION FROM THE EVALUATION	
		
		//------> Get values for the statistics <------
		if (FLAGS & PLOT_DATA){
			order_double (GA_p->Evaluation, GA_p->order, GA_p->n_organisms);	
			Gen_best_ERMS[GA_p->n_generations] = GA_p->Evaluation[GA_p->n_organisms - 1];
			Gen_worst_ERMS[GA_p->n_generations] = GA_p->Evaluation[0] ;
			Gen_ave_ERMS[GA_p->n_generations] = 0;		
		
			for (j = 0; j < GA_p->n_organisms; j++) {
				Gen_ave_ERMS[GA_p->n_generations] +=  GA_p->Evaluation[j];
			}
			Gen_ave_ERMS[GA_p->n_generations]/= GA_p->n_organisms;
		}		
		GA_p->n_generations++;
	}
	
printf("Generations Done \n");

print_GA_generation(GA_p);

	if (FLAGS & PLOT_DATA){
		
		memcpy(graph_GA[1]->data, Gen_best_ERMS, sizeof(double)*GA_p->stop_cond);
		memcpy(graph_GA[2]->data,Gen_ave_ERMS, sizeof(double)*GA_p->stop_cond);
		memcpy(graph_GA[3]->data ,Gen_worst_ERMS, sizeof(double)*GA_p->stop_cond);
		
		mgra.y = graph_GA;
		mgra.n_curves = 2;
		mgra.graph_name = (char *)malloc(30*sizeof(char));
		mgra.curve_names= (char **)malloc(3*sizeof(char*));
		
		mgra.X_axis_n = (char *)malloc(30*sizeof(char));
		mgra.Y_axis_n = (char *)malloc(30*sizeof(char));
		
		for (i = 0; i < 3; i++){
			mgra.curve_names[i] = (char *)malloc(30*sizeof(char));
		}
		sprintf(mgra.curve_names[0], "Best");
		sprintf(mgra.curve_names[1], "Average");
		sprintf(mgra.curve_names[2], "Worst");
		
		sprintf(mgra.graph_name, "GA");
		
		sprintf(mgra.X_axis_n , "Generation");
		sprintf(mgra.Y_axis_n , "Evaluation (RMSE)");



		plot_multiple_graph(&mgra);
		
		
		free(Gen_best_ERMS);
		free(Gen_ave_ERMS);	
		free(Gen_worst_ERMS);
		free(Gen_biodiversity);
		
		for (i = 0; i < 4; i++){
			gsl_vector_free(graph_GA[i]);
		}
		free(graph_GA);
		 
		free(mgra.X_axis_n);
		free(mgra.Y_axis_n);
		free(mgra.graph_name);
		for (i = 0; i < 3; i++){
			free(mgra.curve_names[i]);
		}
		free(mgra.curve_names);
	}
	//--------------> Free SA data <---------------	
	for (i = 0; i < GA_p->n_organisms; i++){
		free(GA_p->Population[i]);
		free(GA_p->Aux_Population[i]);
	}
	free(GA_p->Population);
	free(GA_p->Aux_Population);

	free(GA_p->Evaluation);
	free(GA_p->Fitness);
	free(GA_p->order);	
	free(GA_p->fights_won); // Do it somewhere else
	if (GA_p->n_Selected > 0){
		for (i = 0; i < GA_p->n_Selected; i++){
			free(GA_p->Selected_Population[i]);
		}
		free(GA_p->Selected_Population);		
	}
	
	return 1;
}

int get_Next_Gen(GA_data *GA_p){
	
	get_Fitness(GA_p);				// Get the fitness of the organisms
/*	printf("Poblacion anterior \n");
	for (i = 0; i < GA_p->n_organisms; i++){
		print_bitvector(GA_p->Population[i],GA_p->num_bits );
	}
*/		
	select_organisms(GA_p);			// Select organisms for reproduction and shuffle
/*	printf("Poblacion seleccionada \n");
	for (i = 0; i < GA_p->n_organisms; i++){
		print_bitvector(GA_p->Population[i],GA_p->num_bits );
	}	
*/			
	crossover_population(GA_p);		// Sex them up
/*	printf("Poblacion reproducida\n");
	for (i = 0; i < GA_p->n_organisms; i++){
		print_bitvector(GA_p->Population[i],GA_p->num_bits );
	}
*/	
	mutate_gen(GA_p);				// Mutate them 
/*	printf("Mutados  \n");
	for (i = 0; i < GA_p->n_organisms; i++){
		print_bitvector(GA_p->Population[i],GA_p->num_bits );
	}
*/		
	return 1;
}	

int get_Evaluation(GA_data *GA_p){

	// WITH THREADS
	threads_WorkOut(GA_p->Population,GA_p->Evaluation, GA_p->n_organisms, &Thread_p);
	
	return 1;
}

int get_Fitness2(GA_data *GA_p){		// CHECKING FUNCTION TO TRY THE ALGORTHM
	int i;
	int n_organisms = GA_p->n_organisms;	
	double *Fitness = GA_p->Fitness;
	
	for(i = 0; i < n_organisms; i++){
		Fitness[i] = (GA_p->n_organisms - get_vector_weight(GA_p->Population[i], GA_p->num_bits));
// printf("Ev %i -> Fit %f \n",get_vector_weight(GA_p->Population[i], GA_p->num_bits), Fitness[i]);
	}
	return 1;
}

int get_Fitness3(GA_data *GA_p){
	int i;
	double *Evaluation = GA_p->Evaluation;
	int n_organisms = GA_p->n_organisms;	
	double *Fitness = GA_p->Fitness;
	
	double average_ev = 0;
	
	for(i = 0; i < n_organisms; i++){
		average_ev += 1/Evaluation[i];
	}
	average_ev = average_ev /n_organisms;
			
	
	for(i = 0; i < n_organisms; i++){
		Fitness[i] = (1/Evaluation[i])/average_ev;	// In our case the more EVAL, the worst so we do 1/Eval
// printf("Ev %lf -> Fit %f \n",Evaluation[i], Fitness[i]);
	}

	return 1;
}

int get_Fitness(GA_data *GA_p){
	int i;
	double max_Evaluation;
	double *Evaluation = GA_p->Evaluation;
	int n_organisms = GA_p->n_organisms;	
	double *Fitness = GA_p->Fitness;
	
	double average_ev = 0;
	max_Evaluation = 0;
	for(i = 0; i < n_organisms; i++){
		if (Evaluation[i] > max_Evaluation){
			max_Evaluation = Evaluation[i];
		}
	}
	for(i = 0; i < n_organisms; i++){
		Fitness[i]  = max_Evaluation - Evaluation[i];
	}
	for(i = 0; i < n_organisms; i++){
		average_ev += Fitness[i];
	}
	average_ev = average_ev /n_organisms;
			
	
	for(i = 0; i < n_organisms; i++){
		Fitness[i] = Fitness[i]/average_ev ;	
// printf("Ev %lf -> Fit %f \n",Evaluation[i], Fitness[i]);
	}

	return 1;
}

int select_organisms(GA_data *GA_p){
	// Dublicate the best ones and shuffel them
	unsigned int **aux;

	//----------> Roulette Wheel Selection <--------------- 
	if (GA_p->sel_type == ROULETTE_WHEEL_SEL){
		Roulette_Wheel_Sel(GA_p);
	}
	//----------> Rank Selection <--------------- 
	else if (GA_p->sel_type == RANK_SEL){
		Rank_Sel(GA_p);
	}
	else if (GA_p->sel_type == TOURNAMENT_SEL){
		Tournament_Sel(GA_p);
	}
	
	aux = GA_p->Aux_Population;		// Swap to put the next generation in Population.
	GA_p->Aux_Population = GA_p->Population;
	GA_p->Population = aux;
	return 1;
}

int crossover_population(GA_data *GA_p){
	int i;
	int n_cross = GA_p->n_organisms/2;	// Number of cross overs
	
	int father, mother;				// Index of the father and the mother.	
	for(i = 0; i < n_cross; i++){
//		printf("$$$$$$$$$$$$$$$$$$$$$$$ \n");
		if (rand() <= (unsigned int)((int)(GA_p->Xover_rate*RAND_MAX))){	// Prob of sex
//			printf("-------------------- \n");
			father = 2*i;
			mother = 2*i + 1;
			// We make a 2-point crossover from 2 1-point crossover
			crossover_1p(GA_p->Population[father], GA_p->Population[mother], GA_p->num_bits);
			crossover_1p(GA_p->Population[father], GA_p->Population[mother], GA_p->num_bits);
		}
	}	
	return 1;
}

int mutate_gen(GA_data *GA_p){
	int i;
	for(i = 0; i < GA_p->n_organisms; i++){
		if (rand() < RAND_MAX*GA_p->Mut_rate){
			mutate_bitvector(GA_p->Population[i], GA_p->num_bits, 1 + rand()% GA_p->max_mut);
		}
	}
	return 1;
}	

int crossover_1p(unsigned int *a, unsigned int *b, int num_bits){
	int i;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable
	
	unsigned int aux;
	
	int cross_point = rand()%num_bits;	// Change to rand()
	unsigned int last_mask = 0;			// When we swap bits, we use a mask with which we erase the swap bits from both vectors
							
	unsigned int swap_a = 0;		// Part of the first vector we are swapping
	unsigned int swap_b = 0;		
	
	var_chosen = cross_point/(sizeof(unsigned int)*8);
	bit_chosen = cross_point%(sizeof(unsigned int)*8);
	
	if (var_chosen > 0){
		for (i = 0; i < var_chosen; i++){   // Get full int variables interchanged
			aux = a[i];		// Get the swapped parts
			a[i] = b[i];
			b[i] = aux;	
		}
	}

// printf("Particion %i, %i, %i\n", cross_point, var_chosen,bit_chosen);
	
	for (i = 0; i < bit_chosen; i++){	// Get the mask
		last_mask |= 1 << i;
	}
	swap_a = a[var_chosen]&last_mask;		// Get the last part of the sapped parts
	swap_b = b[var_chosen]&last_mask;	

	a[var_chosen] &= ~last_mask;
	a[var_chosen] |= swap_b; 
	
	b[var_chosen] &= ~last_mask;
	b[var_chosen] |= swap_a; 	
	return 1;
}

int Roulette_Wheel_Sel(GA_data *GA_p){
	int i,j;
	double total_fitness = GA_p->n_organisms;
	double rem_fitness;
	int n_ELITISTS;
//	printf("Doinf roulete\n");
	for (i = 0; i < GA_p->n_organisms; i++){	// For every organism that will reproduce
		rem_fitness = rand()%((int)(1000*total_fitness));		// ROULETTE CHOSING !!!
		rem_fitness /= 1000;
// printf("Selected object  %f \n",rem_fitness);
		
		for (j = 0; j < GA_p->n_organisms; j++){  // We advance through the reoulette to get the selected organism
			rem_fitness -= GA_p->Fitness[j];
			if (rem_fitness <= 0){
				bitvector_cpy(GA_p->Aux_Population[i], GA_p->Population[j], GA_p->num_bits);
// printf("Rem_fitness %f -> selected %i\n",rem_fitness,j);
				break;
			}					
		}
	}
	
			//---> ********** ELITISTS ****************** <----
	//----------> DO SOMETHING FOR THIS NOT TO CONVERGE TO "ALL ELITISTS" COS THEY ARE ALWAYS COPIED
	order_double (GA_p->Fitness, GA_p->order, GA_p->n_organisms); // Order them to get the elitists
	n_ELITISTS = GA_p->n_organisms * GA_p->Elit_rate;
printf("Elitism: %i \n",n_ELITISTS);
	for (i = 0; i < n_ELITISTS; i++){	// For every organism that will reproduce
		bitvector_cpy(GA_p->Aux_Population[rand()%GA_p->n_organisms], GA_p->Population[GA_p->order[i]], GA_p->num_bits);
	}
	return 1;
}
	
int Rank_Sel(GA_data *GA_p){
	int i,j;
	int total_rank;
	int ran_sel_rank;

	order_double (GA_p->Fitness, GA_p->order, GA_p->n_organisms); // Order them to get the rank
	total_rank = (GA_p->n_organisms*(GA_p->n_organisms+1))/2;	// Sumatorio de 1 a n = n*(n+1)/2
	printf("Total rank is:	%i",total_rank);
	for (i = 0; i < GA_p->n_organisms; i++){	// For every organism that will reproduce
		ran_sel_rank = rand()%(total_rank + 1);		// GET the random value !!!
// printf("Selected object  %f \n",ran_sel_rank);
		
		for (j = 0; j < GA_p->n_organisms; j++){  // We advance through the ranked roulette to get the selected organism
			ran_sel_rank -= GA_p->n_organisms -j;//(GA_p->n_organisms - (1 + GA_p->order[j]));
			if (ran_sel_rank <= 0){
				bitvector_cpy(GA_p->Aux_Population[i], GA_p->Population[GA_p->order[j]], GA_p->num_bits);
//printf("Rem_fitness %f -> selected %i\n",rem_fitness,j);
				break;
			}					
		}
	}
	return 1;		
}
	
int Tournament_Sel(GA_data *GA_p){
	int i,j;
	int winners = 0;	// Number of organisms that got to the next generation
	int n_fighers;		// Number of fighters at every tournament 
	int n_champions;	// Number of champions we take to the next tournament
	int pos_fihgters[GA_p->n_organisms];	// Position of the fighters we will select
	
	// WE TAKE GROUPS OF ORGANISMS AND MAKE THE FIGHT WITH EACH OTHER, FROM EVERY GROUP WE TAKE THE BEST ONES
	// THE PARAMETERS OF THE ALGORITHM ARE:			SIZE OF EACH TOURNAMENT
	//												NUMBER OF WINNERS WE TAKE FROM EACH TOURNAMENT
	// AT THE END OF THE DAY WE MUST HAVE REPLACED THE WHOLE GENERATION
	
	while (winners < GA_p->n_organisms){
		n_fighers = GA_p->Min_fighters + (rand () % (GA_p->Max_fighters - GA_p->Min_fighters + 1));	// Multipli by a constant to do the selective pressure

// printf("Fighters %i  | Winners %i \n",n_fighers, winners );
		for (i = 0 ; i < n_fighers; i++){
			pos_fihgters[i] = rand () % GA_p->n_organisms;	// Select random fighters
			GA_p->fights_won[i] = 0;						// Make them initial loosers
		}
		for (i = 0 ; i < n_fighers; i++){			// Make it fight with each other
			for (j = 0 ; j < n_fighers; j++){			// Make it fight with each other
				if (GA_p->Fitness[pos_fihgters[i]] > GA_p->Fitness[pos_fihgters[j]]){
					GA_p->fights_won[i]++;
				}
			}
		}
		ordenar_int ((int *)GA_p->fights_won, GA_p->order, n_fighers);	
		
/* for (i = 0; i < n_fighers; i++) {
	printf("fighter %i with Fitness %f won -> %i \n", pos_fihgters[GA_p->order[i]], GA_p->Fitness[pos_fihgters[GA_p->order[i]]], GA_p->fights_won[i]);
} */
		
		n_champions = 1 + (rand ()% n_fighers)/2;

		if (n_champions + winners >= GA_p->n_organisms){
			n_champions = GA_p->n_organisms - ( winners);
		}
		
		for (i = 0 ; i < n_champions; i++){
			bitvector_cpy(GA_p->Aux_Population[winners], GA_p->Population[pos_fihgters[GA_p->order[i]]], GA_p->num_bits);
			winners++;
		}
	}
	return 1;
}
	//	prob = exp(delta/T);	
void print_GA_data (GA_data *GA_p){
 	 	FILE *pf;
	int nbytes = 0;
	char buffer[150];
	char command[100];
	char result_folder[MAX_CHAR] = "Results";
	char subfolder[MAX_CHAR] = "GA";
	char data_file[MAX_CHAR] = "gen";
	char base_folder[300];
	char mkdir[MAX_CHAR] = "mkdir ";
	char ls[MAX_CHAR] = "ls";
	char folder_name[50];		// Name of the folder with the files
	
	char aux_bitvector[n];
	char date[60];
	char dir[60];
	int i,k;
	time_t fechaActual ;
	struct tm * fechaPtr;
	
	sprintf(data_file + strlen(data_file), "%i", GA_p->n_generations);
	// Get the date.
	fechaActual = time(0) ;		// Get the time and print it
	fechaPtr = gmtime(&fechaActual) ;
	nbytes = sprintf(date, "%i-%i-%i(%i:%i:%i)", fechaPtr->tm_mday, fechaPtr->tm_mon + 1, 
		fechaPtr->tm_year + 1900, fechaPtr->tm_hour, fechaPtr->tm_min, fechaPtr->tm_sec);	
		
	//---------> Make the main folder <------------- 
	if (check_file("",result_folder) == 0){
		sprintf(command,"%s ", mkdir);
		sprintf(command + strlen(command),"%s", result_folder);
		system(command);
	}

//printf("%s \n", command);

	//-------------> Make the base folder <-----------------
	
	if (check_file(result_folder,subfolder) == 0){	
		sprintf(base_folder,"%s", result_folder);
		sprintf(base_folder + strlen(base_folder),"/%s", subfolder);	
		sprintf(command,"%s", mkdir);
		sprintf(command + strlen(command),"%s", base_folder);
		system(command);
	}

	//-----------> Write data into file <-----------------
	sprintf(dir,"%s/%s/%s", result_folder,subfolder,data_file);	
	pf = fopen(dir,"w");
	 if (pf == NULL){
		perror ("fopen");
	}

	nbytes = sprintf(buffer,"Date: %s \n", date);		// Date
	fwrite(buffer,sizeof(char),nbytes,pf);
  
 	nbytes = sprintf(buffer,"N_organisms: %u \n", GA_p->n_organisms);	// Number of organisms
	fwrite(buffer,sizeof(char),nbytes,pf);
	
 	nbytes = sprintf(buffer,"N_generations: %u \n", GA_p->n_generations);	// Number of generations
	fwrite(buffer,sizeof(char),nbytes,pf);
		
 	nbytes = sprintf(buffer,"Bits of the population %i\n\n", GA_p->num_bits);	// Num of parameters
	fwrite(buffer,sizeof(char),nbytes,pf);	

 	nbytes = sprintf(buffer,"Population: \n");	// Generation in "uint" form
	fwrite(buffer,sizeof(char),nbytes,pf);	
	for (k = 0; k < GA_p->n_organisms; k++){		
		get_bitvector(GA_p->Population[k], n, aux_bitvector); 
		nbytes = sprintf(buffer,"\n %s",aux_bitvector);  // Selected parameters in "bitvector" form
		fwrite(buffer,sizeof(char),nbytes,pf);	   

	}

    fclose(pf);
}

int Init_GA_Population(GA_data *GA_p){
	int i;
	if (GA_p->n_Selected == 0) {
		for (i = 0; i < GA_p->n_organisms; i++){
			zero_bitvector(GA_p->Population[i], GA_p->num_bits);
			get_random_bitarray(GA_p->Population[i], GA_p->num_bits);			// Initialize values of population
		}	
	}	
	else {		// If there is initial selected population
		for (i = 0; i < GA_p->n_Selected; i++){
			bitvector_cpy(GA_p->Population[i], GA_p->Selected_Population[i], GA_p->num_bits);
		}			
		for (i = i; i < GA_p->n_organisms; i++){
			zero_bitvector(GA_p->Population[i], GA_p->num_bits);
			get_random_bitarray(GA_p->Population[i], GA_p->num_bits);			// Initialize values of population
		}			
		
	}
}

int print_organisms(GA_data *GA_p){
	unsigned int j;
	for (j = 0; j < GA_p->n_organisms; j++) {
		print_bitvector(GA_p->Population[j],GA_p->num_bits);	
		printf("Energy: %f \n", GA_p->Evaluation[j]);
	}
}

int print_GA_values(GA_data *GA_p){
		printf("**********  Genetic Algorithm values ***************** \n");
		printf("N organism: %i \n",GA_p->n_organisms);		// Number of parallel tries Lingotes
		printf("N generatons: %i \n",GA_p->n_generations);	// Number of coldings.
		printf("N bits: %i \n",GA_p->num_bits);		// Number of bits of the population	
		printf("Stop cond: %i \n",GA_p->stop_cond);	// Number of coldings.
				
		printf("Mut rate %f \n",GA_p->Mut_rate);		// Simulated anneling temperature
		printf("Crossover rate  : %f \n",GA_p->Xover_rate);			// Simulated constant  temp[t+1] = temp[t]*k
		printf("Elitism rate  : %f \n",GA_p->Elit_rate);			// Simulated constant  temp[t+1] = temp[t]*k
		
		printf("Max mutations: %i \n",GA_p->max_mut);			// Number of maximum mutations per freezing
		printf("N Selected  : %i \n",GA_p->n_Selected);		// Number of selected ingots for initial population
		printf("********** ****************************************** \n");
		return 1;
}

int print_GA_generation(GA_data *GA_p){
		ELM_p.FLAGS = PLOT_DATA_F;
		init_threads(&ELM_p, &Thread_p);
		threads_WorkOut(GA_p->Population, GA_p->Evaluation, GA_p->n_organisms, &Thread_p);		
		return 1;
}


