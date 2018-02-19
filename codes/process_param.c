

#include "../ELMheader.h"


int get_ELM_params (char *config_file, Input_params *Input_p){
	FILE *pf;
	char conf_word[MAX_CONF_SIZE];
	char value[MAX_BITS];

	pf = fopen(config_file,"r");	// Open configuration file "read" mode
	if (pf == NULL){
		perror("fopen");
		printf("Fichero de configuracion invalido \n");
		exit(-1);
	}

	// File with the testing input vectors
	// File with the testing results

	// File with the training input vectors
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	strcpy(Input_p->train_input_dir, value);
	
	// File with the training output vectors
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	strcpy(Input_p->train_output_dir, value);
	
	// File with the testing input vectors
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	
	if (strcmp(value,"NONE") != 0){		// If we are given testing vectors
		strcpy(Input_p->test_input_dir, value);	
		fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
		fscanf(pf,"%s",value);		// Read the value
		strcpy(Input_p->test_output_dir, value);
	}
	else {
		Input_p->test_input_dir[0] = 0;
		fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
		fscanf(pf,"%s",value);		// Read the value	
		Input_p->test_output_dir[0] = 0;
	}
	
	// Get Number of training vectors
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	Input_p->Ntrain = atoi(value);
	
	// Get Number of testing vectors
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	Input_p->Ntest = atoi(value);

	// Get Number of input parameters
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	Input_p->n_input = atoi(value);
	
	// Get Number of output parameters
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	Input_p->n_output = atoi(value);
	
	// Number of ELM to do to get the average RMSE
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	Input_p->n_ELMs = atoi(value);
		
	// Number of Neurons
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	Input_p->Nh = atoi(value);

	// Get the Activation function.
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	Input_p->activation_f = atoi(value);
	
	// Number of Threads
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	Input_p->n_threads = atoi(value);
		
	fclose(pf);
	return 1;
}

int get_ELMs_params (char *config_file, ELMs_data *ELMs_p){
	int i,j;
	FILE *pf;
	char conf_word[MAX_CONF_SIZE];
	char value[MAX_BITS];

	pf = fopen(config_file,"r");	// Open configuration file "read" mode
	if (pf == NULL){
		perror("fopen");
		printf("Fichero de configuracion invalido \n");
		exit(-1);
	}
	
	// Number of Selected vectors  to do.
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	ELMs_p->n_Selected = atoi(value);


	// Get selected ELM value to do
	ELMs_p->Selected = (unsigned int **) malloc (sizeof(unsigned int *)*ELMs_p->n_Selected);
	for (i = 0; i < ELMs_p->n_Selected ; i++){		// $$$$$$$$$$ EXTERN VARIABLE $$$$$$$$$$$$
		ELMs_p->Selected[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
	}	
	for (i = 0; i < ELMs_p->n_Selected; i++){	
		fscanf(pf,"%s",value);
		get_u_vector(ELMs_p->Selected[i],n ,value);
	}	
	
	// Number of ELM to do.
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	ELMs_p->n_ELMs = atoi(value);


	
	// Get FLAGS
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	ELMs_p->FLAGS = atoi(value);
	
	// Will we print the random distribution 
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	
	if (strcmp(value,"yes") == 0){		// If we are given testing vectors
		ELMs_p->get_distrib = 1;
	} else {
		ELMs_p->get_distrib = 0;
	}
	// Get number of distributions for making the statistic
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	ELMs_p->n_ELMs_dis = atoi(value);

	// Get the number of divs
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	ELMs_p->n_divs = atoi(value);

	fclose(pf);
	return 1;
}

int get_ES_params (char *config_file, ES_data *ES_p){
	int i,j;
	FILE *pf;
	char conf_word[MAX_CONF_SIZE];
	char value[MAX_BITS];

	pf = fopen(config_file,"r");	// Open configuration file "read" mode
	if (pf == NULL){
		perror("fopen");
		printf("Fichero de configuracion invalido \n");
		exit(-1);
	}
	
	// Minimum number of 1s
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	ES_p->num_1s_min = atoi(value);

	// Maximum number of 1s
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	ES_p->num_1s_max = atoi(value);

	// Number of Top Solutions
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	ES_p->n_TOP = atoi(value);

	// Selected Es bitarray
	ES_p->ES_selected = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label		
	fscanf(pf,"%s",value);
	get_u_vector(ES_p->ES_selected, n ,value);
	
	// Fixed 1s bitarray
	ES_p->ES_fixed_1s = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label		
	fscanf(pf,"%s",value);
	get_u_vector(ES_p->ES_fixed_1s, n ,value);
		
	fclose(pf);
	return 1;
}

int get_GA_params (char *config_file, GA_data *GA_p){
	int i,j;
	FILE *pf;
	char conf_word[MAX_CONF_SIZE];
	char value[MAX_BITS];
	pf = fopen(config_file,"r");	// Open configuration file "read" mode
	if (pf == NULL){
		perror("fopen");
		printf("Fichero de configuracion invalido \n");
		exit(-1);
	}
	
	// Get the Number of organisms
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	GA_p->n_organisms = atoi(value);
	
	// Get the Number of genrations
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	GA_p->stop_cond = atoi(value);
	GA_p->n_generations = 0;
	
	// Probability of mutation
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	GA_p->Mut_rate = atof(value); 
	
	// Probability of crossover
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	GA_p->Xover_rate = atof(value); 	
	
	// Elistims proportionality
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	GA_p->Elit_rate = atof(value); 	

	// Max number of mutations
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	GA_p->max_mut = atoi(value); 	
		
	// Selection algorithm
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	GA_p->sel_type = atoi(value); 	

	// Min number of fighters
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	GA_p->Min_fighters = atoi(value); 	
		
	// Max number of fighters
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	GA_p->Max_fighters = atoi(value); 	

	// Initial population
	fscanf(pf,"%s",conf_word);		// This first read is to eliminate the label
	fscanf(pf,"%s",value);			
	GA_p->n_Selected = atoi(value);
	if (GA_p->n_Selected  == 0) {
	}
	else{
		GA_p->Selected_Population = (unsigned int **) malloc (sizeof(unsigned int *)*GA_p->n_Selected );
		for (i = 0; i < GA_p->n_Selected; i++){		// $$$$$$$$$$ EXTERN VARIABLE $$$$$$$$$$$$
			GA_p->Selected_Population[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
		}	
		for (i = 0; i < GA_p->n_Selected; i++){		
			fscanf(pf,"%s",value);
			get_u_vector(GA_p->Selected_Population[i],n , value);
		}	
	}
	
	GA_p->num_bits = n;
	fclose(pf);
	return 1;
}

int get_Graph_params (char *config_file, linear_graph *l_graph, linear_graph *l_graph2, unsigned int * Graph_selected){
	int i;
	FILE *pf;
	char conf_word[MAX_CONF_SIZE];
	char value[MAX_BITS];

	pf = fopen(config_file,"r");	// Open configuration file "read" mode
	if (pf == NULL){
		perror("fopen");
		printf("Fichero de configuracion invalido \n");
		exit(-1);
	}
	
	// ----------> FIRST GRAPH <----------
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	
	if (strcmp(value,"yes") == 0){}
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label	
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph->init_Xvalue = 	atoi(value);	// X-Axis parameters
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph->end_Xvalue = 	atoi(value);
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph->XstepSize = 	atoi(value);		
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph->init_Yvalue =	atoi(value);	// Y-Axis parameters
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph->end_Yvalue = 	atoi(value);		
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph->Yn_curves = 	atoi(value);	
	
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);	
	get_u_vector(Graph_selected,n , value);
	
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph->n_average = 	atoi(value);
	

// ----------> SECOND GRAPH <----------
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	
	if (strcmp(value,"yes") == 0){}
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph2->init_Xvalue = 	atoi(value);	// X-Axis parameters
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph2->end_Xvalue = 	atoi(value);
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph2->XstepSize = 	atoi(value);		
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph2->init_Yvalue =	atoi(value);	// Y-Axis parameters
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph2->end_Yvalue = 	atoi(value);		
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph2->Yn_curves = 	atoi(value);	
	
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value

	
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	l_graph2->n_average = 	atoi(value);
	
	fclose(pf);
	return 1;
}

int get_SA_params (char *config_file, SA_data *SA_p){
	int i,j;
	FILE *pf;
	char conf_word[MAX_CONF_SIZE];
	char value[MAX_BITS];

	pf = fopen(config_file,"r");	// Open configuration file "read" mode
	if (pf == NULL){
		perror("fopen");
		printf("Fichero de configuracion invalido \n");
		exit(-1);
	}
	
	// Get the Number of Ingots
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	SA_p->n_ingots = atoi(value);

	// Get the Number of freezings
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	SA_p->stop_cond = atoi(value);
	SA_p->n_freezings = 0;

	// Get the Number of max mutations
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	SA_p->max_mut = atoi(value);
			
	// Initial temperature
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);	
	SA_p->temp = atof(value);	
	
	// Decreasing konstant
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);	
	SA_p->k = atof(value);
	
	if ((SA_p->k <= 0.0) ||(SA_p->k > 1.0)){	// If there are too many threads they will require too much memmory
		printf("Decreasing constant must be [0 - 1] \n");
		exit(-1);	
	}	

	// Initial population
	fscanf(pf,"%s",conf_word);		// This first read is to eliminate the label
	fscanf(pf,"%s",value);			
	SA_p->n_Selected = atoi(value);
	if (SA_p->n_Selected  == 0) {
	}
	else{
		SA_p->Selected_Ingots = (unsigned int **) malloc (sizeof(unsigned int *)*SA_p->n_Selected );
		for (i = 0; i < SA_p->n_Selected; i++){		
			SA_p->Selected_Ingots[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
		}	
		for (i = 0; i < SA_p->n_Selected; i++){		
			fscanf(pf,"%s",value);
			get_u_vector(SA_p->Selected_Ingots[i], n ,value);
		}	
	}

	SA_p->num_bits = n;

	fclose(pf);
	return 1;
}

int get_AIS_params (char *config_file, AIS_data *AIS_p){
	FILE *pf;
	char conf_word[MAX_CONF_SIZE];
	char value[MAX_BITS];

	pf = fopen(config_file,"r");	// Open configuration file "read" mode
	if (pf == NULL){
		perror("fopen");
		printf("Fichero de configuracion invalido \n");
		exit(-1);
	}
	
	// Number of clonning cycles
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	AIS_p->stop_cond = atoi(value);
	
	// Initial number of B cells
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	AIS_p->Ini_B = atoi(value);
	
	// Maximum number of B cells
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	AIS_p->Max_B = atoi(value);	
	
	// Number of B clones
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	AIS_p->n_clones = atoi(value);		
	
	// Number of random Ab to get new B cells
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	AIS_p->n_Ab = atoi(value);		
	
	// B domain. Minimum hamming distance between B-cells
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	AIS_p->B_domain = atoi(value);	
	
	// Number of cycles a B cell getting upgraded till it becomes a M cell
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	AIS_p->B_lifetime = atoi(value);	

	// Maximum number of Master cells
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	AIS_p->Max_M = atoi(value);	
				
	AIS_p->n_cycles = 0;	
	AIS_p->n_B = 0;
	AIS_p->n_M_cells = 0;
	AIS_p->Beato_Defense = 1000000;	
	AIS_p->num_bits = n;

	fclose(pf);
	return 1;
}

int get_FP_params (char *config_file, FP_data *FP_p){
	int i,j;
	FILE *pf;
	char conf_word[MAX_CONF_SIZE];
	char value[MAX_BITS];

	pf = fopen(config_file,"r");	// Open configuration file "read" mode
	if (pf == NULL){
		perror("fopen");
		printf("Fichero de configuracion invalido \n");
		exit(-1);
	}
	FP_p->Operations = 0;
	FP_p->n_F = n;
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value

	if (strcmp(value,"yes") == 0){		// If we are given testing vectors
		FP_p->Operations |= FS_Fitness	;
	} 
	
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	
	if (strcmp(value,"yes") == 0){		// If we are given testing vectors
		FP_p->Operations |= FS_Affinity;
	} 
	
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	
	if (strcmp(value,"yes") == 0){		// If we are given testing vectors
		FP_p->Operations |= SFS;
	} 
	
	fscanf(pf,"%s",conf_word);	// This first read is to eliminate the label
	fscanf(pf,"%s",value);		// Read the value
	
	if (strcmp(value,"yes") == 0){		// If we are given testing vectors
		FP_p->Operations |= SBE;
	} 
	fclose(pf);
	return 1;
}


