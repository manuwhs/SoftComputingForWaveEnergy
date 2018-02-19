#include "../ELMheader.h"
// In the Simmulad anneling --> Allow to specify a factor to multiply the auto temperature value

int AIS_op(AIS_data *AIS_p, int FLAGS){
	unsigned int i,j;
	int non_stop = 1;
	unsigned int max_Ab = AIS_p->Max_B * AIS_p->n_clones + AIS_p->n_Ab;
	//------------> Variables for statistics <------------
	double *	Gen_best_ERMS;	// Array with the averBe ERMS of the best ingots so far
	double *	Gen_ave_ERMS;	// Array with the averBe ERMS  of the ingots

	multiple_graph mgra;
	gsl_vector **graph_GA;

	if (FLAGS & PLOT_DATA){
	//--------> Reserve memmory for the generation statistics <------------
		Gen_best_ERMS = (double *) malloc (sizeof(double)*AIS_p->stop_cond);
		Gen_ave_ERMS = (double *) malloc (sizeof(double)*AIS_p->stop_cond);	
		
	}	
	//--------> Reserve memmory for the operations <------------
	AIS_p->B_cells = (unsigned int **) malloc (sizeof(unsigned int *)*AIS_p->Max_B);
	for (i = 0; i < AIS_p->Max_B; i++){
		AIS_p->B_cells[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
	}
	
	AIS_p->Ab_cells = (unsigned int **) malloc (sizeof(unsigned int *)*max_Ab);
	for (i = 0; i < max_Ab; i++){
		AIS_p->Ab_cells[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
	}
	
	AIS_p->M_cells = (unsigned int **) malloc (sizeof(unsigned int *)*AIS_p->Max_M);
	for (i = 0; i < AIS_p->Max_M ; i++){
		AIS_p->M_cells[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
	}	
	AIS_p->Beato_cell = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
	
	AIS_p->B_rem_lifes = (unsigned int *) malloc (sizeof(unsigned int )*AIS_p->Max_B);
	AIS_p->order  = (int *) malloc (sizeof(int )*max_Ab);
	
	AIS_p->B_Defense = (double *) malloc (sizeof(double)*AIS_p->Max_B);
	AIS_p->Ab_Defense = (double *) malloc (sizeof(double)*max_Ab);	
	AIS_p->M_Defense = (double *) malloc (sizeof(double)*AIS_p->Max_M);	
	//--------> Generate initial set of B and Ab cells<------------
	
	Init_AIS (AIS_p);	
	print_AIS_values(AIS_p);
//--------> Start Cloning  <------------	
	for (i = 0; i < AIS_p->stop_cond; i++){
printf("Cycle %i, Num_B: %i, Num_M %i \n", i,AIS_p->n_B, AIS_p->n_M_cells);

print_AIS_cells(AIS_p);

		clone_Bs(AIS_p);
//printf("Cloned \n");
		mutate_Abs(AIS_p);
//printf("Mutated \n");	
		Abs_defense(AIS_p);
//printf("Got defense \n");	
		non_stop = upgrade_Bs(AIS_p);

		//------> Get values for the statistics <------
		if (FLAGS & PLOT_DATA){
			Gen_best_ERMS[AIS_p->n_cycles]= 0;
			Gen_ave_ERMS[AIS_p->n_cycles] = 0;	
			if (AIS_p->n_M_cells > 0){
				Gen_best_ERMS[AIS_p->n_cycles]= AIS_p->Beato_Defense;
				for (j = 0; j < AIS_p->n_M_cells; j++) {
					Gen_ave_ERMS[AIS_p->n_cycles] +=  AIS_p->M_Defense[j];
				}
				Gen_ave_ERMS[AIS_p->n_cycles]/= AIS_p->n_M_cells;
			}
			
		}		
		AIS_p->n_cycles++;
		
		if (non_stop == 0){
			break;
		}
	}
	
	if (AIS_p->n_cycles == AIS_p->stop_cond) {	// If it ended coz stop cpnd num of cycles
		for (i = AIS_p->n_M_cells; i = AIS_p->Max_M; i++ ){
			bitvector_cpy(AIS_p->M_cells[AIS_p->n_M_cells],AIS_p->B_cells[i], AIS_p->num_bits);
			AIS_p->M_Defense[AIS_p->n_M_cells] = AIS_p->B_Defense[i];		
		}
	}
	
printf("Generaciones acabadas \n");
	print_AIS_stuff (AIS_p);
	print_AIS_data(AIS_p);
	print_input_statistics (AIS_p->M_cells, n, AIS_p->n_M_cells);
	
	if (FLAGS & PLOT_DATA){
		
		graph_GA = (gsl_vector **)malloc(3*sizeof(gsl_vector *));
		for (i = 0; i < 3; i++){
			graph_GA[i] = gsl_vector_alloc(AIS_p->n_cycles);
		}
		for (i = 0; i < AIS_p->n_cycles; i++){
			graph_GA[0]->data[i] = i;
		}		
		
		memcpy(graph_GA[1]->data, Gen_best_ERMS, sizeof(double)*AIS_p->n_cycles);
		memcpy(graph_GA[2]->data,Gen_ave_ERMS, sizeof(double)*AIS_p->n_cycles);
		
		mgra.y = graph_GA;
		mgra.n_curves = 2;
		mgra.graph_name = (char *)malloc(30*sizeof(char));
		mgra.curve_names= (char **)malloc(2*sizeof(char*));
		
		mgra.X_axis_n = (char *)malloc(30*sizeof(char));
		mgra.Y_axis_n = (char *)malloc(30*sizeof(char));
		
		for (i = 0; i < 2; i++){
			mgra.curve_names[i] = (char *)malloc(30*sizeof(char));
		}
		sprintf(mgra.curve_names[0], "Best M-Cell (Beato-Cell)");
		sprintf(mgra.curve_names[1], "Average M-Cell");
		
		sprintf(mgra.graph_name, "AIS");
		
		sprintf(mgra.X_axis_n , "Clonnings");
		sprintf(mgra.Y_axis_n , "Evaluation");


		plot_multiple_graph(&mgra);
		
		free(Gen_best_ERMS);
		free(Gen_ave_ERMS);	
		
		for (i = 0; i < 3; i++){
			gsl_vector_free(graph_GA[i]);
		}
		free(graph_GA);
		 
		free(mgra.X_axis_n);
		free(mgra.Y_axis_n);
		free(mgra.graph_name);
		for (i = 0; i < 2; i++){
			free(mgra.curve_names[i]);
		}
		free(mgra.curve_names);
	}
	
	for (i = 0; i < max_Ab ; i++){
		free(AIS_p->Ab_cells[i]);
	}
	free(AIS_p->Ab_cells);
	
	for (i = 0; i < AIS_p->Max_B ; i++){
		free(AIS_p->B_cells[i] );
	}
	free(AIS_p->B_cells );
	
	for (i = 0; i < AIS_p->Max_M ; i++){
		free(AIS_p->M_cells[i]);
	}	
	free(AIS_p->M_cells);
	
	
	free(AIS_p->B_rem_lifes);  // If we let any of these we die
	free(AIS_p->order);
	
	free(AIS_p->B_Defense);
	free(AIS_p->Ab_Defense);
	free(AIS_p->M_Defense);	
	
	free(AIS_p->Beato_cell);
	
	return 1;
}

int Init_AIS (AIS_data *AIS_p){
	unsigned int i;
	unsigned int max_Ab = AIS_p->Max_B*AIS_p->n_clones + AIS_p->n_Ab;
		
	// Get initial random values for the Abs	
	for (i = 0; i < max_Ab; i++){
		zero_bitvector(AIS_p->Ab_cells[i], AIS_p->num_bits);
		get_random_bitarray(AIS_p->Ab_cells[i], AIS_p->num_bits);		// Initialize values of population
	}		
	
	// Get initial defenses
	threads_WorkOut(AIS_p->Ab_cells,AIS_p->Ab_Defense, max_Ab, &Thread_p);
	
	order_double (AIS_p->Ab_Defense,  AIS_p->order,  max_Ab);	// Order Abs to get the best ones
	AIS_p->n_B = 0;
	for (i = 0; i < AIS_p->Ini_B; i++) {		// Select as B the best Ab values
		bitvector_cpy(AIS_p->B_cells[i],AIS_p->Ab_cells[AIS_p->order[(max_Ab -1) - i]], AIS_p->num_bits);
		AIS_p->B_Defense[i] = AIS_p->Ab_Defense[AIS_p->order[(max_Ab -1) - i]];
//	printf("Mierda !! %f %f \n",AIS_p->B_Defense[i],AIS_p->Ab_Defense[AIS_p->order[(max_Ab -1) - i]] );
		AIS_p->B_rem_lifes[i] = AIS_p->B_lifetime;
		AIS_p->n_B++;
	}
	// PORQUE AQUI n_B no es igual a i ???????????????
	return 1;
}

int clone_Bs(AIS_data *AIS_p){
	unsigned int i,j;
	for (i = 0; i < AIS_p->n_B; i++){
		for (j = 0; j < AIS_p->n_clones; j++){
			bitvector_cpy(AIS_p->Ab_cells[i*AIS_p->n_clones + j],AIS_p->B_cells[i], AIS_p->num_bits);
		}
	}
	return 1;
}

int mutate_Abs(AIS_data *AIS_p){
	unsigned int i,j;
	int Ab_pool = AIS_p->n_B*AIS_p->n_clones;
	for (i = 0; i < AIS_p->n_B; i++){											// Mutate B clones
		for (j = 0; j < AIS_p->n_clones; j++){
			mutate_bitvector(AIS_p->Ab_cells[i*AIS_p->n_clones + j], AIS_p->num_bits, 1 + rand()%AIS_p->B_domain);
		}
	}
	// This second one is in order to get new B cells, we generate random Abs -> New Defenses
	for (i = 0; i < AIS_p->n_Ab; i++){	
		get_random_bitarray(AIS_p->Ab_cells[i + Ab_pool], AIS_p->num_bits);
	}
	return 1;		
}

int Abs_defense(AIS_data *AIS_p){
	threads_WorkOut(AIS_p->Ab_cells, AIS_p->Ab_Defense, AIS_p->n_B*AIS_p->n_clones + AIS_p->n_Ab, &Thread_p);
	
	return 1;
}

int remove_B(AIS_data *AIS_p, unsigned int B_pos){
	unsigned int i;
	if (B_pos < AIS_p->n_B){
		for (i = B_pos; i < AIS_p->n_B - 1; i++){	// Bitvectores 
			bitvector_cpy(AIS_p->B_cells[i],AIS_p->B_cells[i+1], AIS_p->num_bits);
		}
		for (i = B_pos; i < AIS_p->n_B - 1; i++){	// Lifes
			AIS_p->B_rem_lifes[i] = AIS_p->B_rem_lifes[i+1];
		}	
		for (i = B_pos; i < AIS_p->n_B - 1; i++){	// Defense
			AIS_p->B_Defense[i] = AIS_p->B_Defense[i+1];
		}		
		
	}
	AIS_p->n_B--;
	return 1;
}

int	upgrade_Bs(AIS_data *AIS_p){
	int i,j;
	unsigned int Ab_pool = AIS_p->n_B* AIS_p->n_clones; // Where the pool of Abs start
	unsigned int Init_B = AIS_p->n_B;
	int flB;
	unsigned int aux_B;
	int deleted[AIS_p->n_B];  // This has the list of deleted objects for B coincidences

	for (i = 0; i < AIS_p->n_B; i++){
		deleted[i] = 0;
	}
	
	// First we check if any of the clones is better than the father
	for (i = 0; i < AIS_p->n_B; i++){
		AIS_p->B_rem_lifes[i]--;			// Probably change place !!!
		for (j = 0; j < AIS_p->n_clones; j++){
			if(AIS_p->Ab_Defense[i*AIS_p->n_clones + j] < AIS_p->B_Defense[i]){	// Upgrade B_cell
		
//printf("************************************* \n");	
//printf("B cell number %i has been upgraded %i \n",i,AIS_p->B_lifetime);
//printf("************************************* \n");	

				bitvector_cpy(AIS_p->B_cells[i],AIS_p->Ab_cells[i*AIS_p->n_clones + j], AIS_p->num_bits);
				AIS_p->B_Defense[i] = AIS_p->Ab_Defense[i*AIS_p->n_clones + j];
				AIS_p->B_rem_lifes[i] = AIS_p->B_lifetime;										// Restart B lifes
			}
		}
	}

	// Now we have to eliminate the B cells that are in the same domain
	// This might happen because two B seeds in different domains, if they upgrade, maybe the 2 upgraded B cells 
	// are in the same domain so we only keep the best one. THE WAY THIS IS DONE CAN BE IMPROVED
	
	for (i = 0; i < AIS_p->n_B; i++){		// For every B
//	printf("checking coincidence for %i \n",i);
		for (j = 0; j < AIS_p->n_B; j++){	// We check if it is close to other B
			if (j != i){					// If its not the same cell
			if (hamming_distance (AIS_p->B_cells[i], AIS_p->B_cells[j], AIS_p->num_bits) <= AIS_p->B_domain){
				if (AIS_p->B_Defense[i] < AIS_p->B_Defense[j]){
					deleted[j] = 1;
					AIS_p->B_rem_lifes[i] = AIS_p->B_lifetime;
/*					
printf("************************************* \n");	
printf("B cell number %i has been deleted for %i \n",j,i);
print_bitvector(AIS_p->B_cells[j],AIS_p->num_bits);
print_bitvector(AIS_p->B_cells[i],AIS_p->num_bits);
printf("************************************* \n");
*/
				}
				else{
					deleted[i] = 1;
					AIS_p->B_rem_lifes[j] = AIS_p->B_lifetime;
/*						
printf("************************************* \n");	
printf("B cell number %i has been deleted for %i \n",i,j);
print_bitvector(AIS_p->B_cells[i],AIS_p->num_bits);
print_bitvector(AIS_p->B_cells[j],AIS_p->num_bits);
printf("************************************* \n");
*/
				}				
			}
			}
		}
	}
	// Get them removed !!! We have to go from last to least because the list decreses otherwise
	aux_B = AIS_p->n_B;
	for (i = aux_B  - 1; i >= 0 ; i--){
		if (deleted[i] == 1) {
			remove_B(AIS_p, i);
//printf("************************************* \n");	
//printf("B cell number %i died \n",i);
//printf("************************************* \n");
		}
	}
/*	
	printf("Before adding new %i\n", AIS_p->n_B);
	for (i = 0; i < AIS_p->n_B; i++){
		print_bitvector(AIS_p->B_cells[i],AIS_p->num_bits);
		printf("-------------> Lives %i \n",AIS_p->B_rem_lifes[i]);
	}	
*/	
	// Now we have to eliminate the B cells whose lifes are extinct --> Presumably coz they cannot improve
	for (i = 0; i < AIS_p->n_B; i++){
		deleted[i] = 0;
	}
	for (i = 0; i < AIS_p->n_B; i++){	
		if (AIS_p->B_rem_lifes[i] <= 0) {
			bitvector_cpy(AIS_p->M_cells[AIS_p->n_M_cells],AIS_p->B_cells[i], AIS_p->num_bits);
			AIS_p->M_Defense[AIS_p->n_M_cells] = AIS_p->B_Defense[i];
			if (AIS_p->M_Defense[AIS_p->n_M_cells] < AIS_p->Beato_Defense){
				bitvector_cpy(AIS_p->Beato_cell,AIS_p->B_cells[i], AIS_p->num_bits);
				AIS_p->Beato_Defense = AIS_p->M_Defense[AIS_p->n_M_cells] ;
			}
			AIS_p->n_M_cells++;
			if (AIS_p->n_M_cells == AIS_p->Max_M){
				printf(" All M cells set \n");
				return 0;				// Indicate the algorithm is over.
			}
			deleted[i] = 1;
		}
	}
	
	aux_B = AIS_p->n_B;
	for (i = aux_B  - 1; i  >= 0 ; i--){
		if(deleted[i] == 1){
			remove_B(AIS_p, i);
//printf("************************************* \n");	
//printf("B cell number %i has been extinguished \n",i);
//printf("************************************* \n");
			
		}
	}
	// Now we try to get a new B cell so that we can choose new stuff
	// We will order the Ab (not B cloned) pool from best to worst and see if we can add a new B
	// For this porpuse:
	//	- The Ab value must be considered good enough
	//	- There must be B cells slots available
	//	- It must be separed from the already existing B cells

/*	
	printf("Before adding new %i\n", AIS_p->n_B);
	for (i = 0; i < AIS_p->n_B; i++){
		print_bitvector(AIS_p->B_cells[i],AIS_p->num_bits);
		printf("-------------> Lives %i \n",AIS_p->B_rem_lifes[i]);
	}	
*/	
/*
	if (AIS_p->n_B < AIS_p->Max_B){		// We add just one
		bitvector_cpy(AIS_p->B_cells[AIS_p->n_B],AIS_p->Ab_cells[Ab_pool + 0], AIS_p->num_bits);
		AIS_p->B_Defense[AIS_p->n_B] = AIS_p->Ab_Defense[Ab_pool + 0];
		AIS_p->B_rem_lifes[AIS_p->n_B] = AIS_p->B_lifetime;
		AIS_p->n_B++; 
	}
*/	
/*	
	printf("After adding new %i\n", AIS_p->n_B);
	for (i = 0; i < AIS_p->n_B; i++){
		print_bitvector(AIS_p->B_cells[i],AIS_p->num_bits);
		printf("-------------> Lives %i \n",AIS_p->B_rem_lifes[i]);
	}	
*/	
	
//	printf("B %i MaxB %i M %i MaxM %i \n",AIS_p->n_B, AIS_p->Max_B,   AIS_p->n_M_cells, AIS_p->Max_M);
	// We introduce the best species from the pool --> The bigger the pool --> The better initially 
	if ((AIS_p->n_B < AIS_p->Max_B)&&		// If we can have more B cells
	   (AIS_p->n_B + AIS_p->n_M_cells < AIS_p->Max_M))	// If we can have more M-Cells
	{		// We add just one
		order_double (&AIS_p->Ab_Defense[Ab_pool],  AIS_p->order,  AIS_p->n_Ab);	// Order Ab pool from best to worst
			bitvector_cpy(AIS_p->B_cells[AIS_p->n_B],AIS_p->Ab_cells[Ab_pool + AIS_p->order[(AIS_p->n_Ab -1)]], AIS_p->num_bits);
			AIS_p->B_Defense[AIS_p->n_B] = AIS_p->Ab_Defense[Ab_pool + AIS_p->order[(AIS_p->n_Ab -1)]];
			AIS_p->B_rem_lifes[AIS_p->n_B] = AIS_p->B_lifetime;
			AIS_p->n_B++; 
printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");	
printf("New B has been added\n");
printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n");
		
	}
	return 1;
}

int print_AIS_values(AIS_data *AIS_p){
		printf("**********  Artificial Inmune System values ***************** \n");
		
		printf("N cycles: %i \n",AIS_p->n_cycles);	
		printf("N bits: %i \n",AIS_p->num_bits);		
		printf("Stop cond: %i \n",AIS_p->stop_cond);	
				
		printf("B-cell Initial %i \n",AIS_p->Ini_B);	
		printf("B-cell Max %i \n",AIS_p->Max_B);	
		printf("B-cell Current %i \n",AIS_p->n_B);	
		printf("B-cell Domain %i \n",AIS_p->B_domain);	
		printf("B-cell Lifetime %i \n",AIS_p->B_lifetime);	
				
		printf("N clones: %i \n",AIS_p->n_clones);			
		printf("N Ab random pool : %i \n",AIS_p->n_Ab);	
		
		printf("Max M cells %i \n",AIS_p->Max_M);			
		printf("Current M cells %i \n",AIS_p->n_M_cells);	
		
		printf("********** ****************************************** \n");
		return 1;
}

int print_AIS_cells(AIS_data *AIS_p){
	unsigned int j;
	printf("**********  Current B-Cells***************** \n");
	for (j = 0; j < AIS_p->n_B; j++) {
		print_bitvector(AIS_p->B_cells[j],AIS_p->num_bits);	
		printf("Defense %f -> Lives %i\n", AIS_p->B_Defense[j],AIS_p->B_rem_lifes[j]);

	}
	/*
	printf("**********  Current M-Cells***************** \n");
	for (j = 0; j < AIS_p->n_M_cells; j++) {
		print_bitvector(AIS_p->M_cells[j],AIS_p->num_bits);	
		printf("Defense: %f\n", AIS_p->M_Defense[j]);
	}
	return 1;
	*/ 
}

int print_AIS_data(AIS_data *AIS_p){
		ELM_p.FLAGS = PLOT_DATA_F;
		init_threads(&ELM_p, &Thread_p);
		threads_WorkOut(AIS_p->M_cells, AIS_p->M_Defense, AIS_p->n_M_cells, &Thread_p);	
	return 1;
}

void print_AIS_stuff (AIS_data *AIS_p){
 	FILE *pf;
	int nbytes = 0;
	char buffer[150];
	char command[100];
	char result_folder[MAX_CHAR] = "Results";
	char subfolder[MAX_CHAR] = "AIS";
	char data_file[MAX_CHAR] = "data";
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

	nbytes = sprintf(buffer,"Date: %s \n", date);		
	fwrite(buffer,sizeof(char),nbytes,pf);
  
 	nbytes = sprintf(buffer,"N cycles: %i \n",AIS_p->n_cycles);	
	fwrite(buffer,sizeof(char),nbytes,pf);
	
 	nbytes = sprintf(buffer,"B-cell Initial %i \n",AIS_p->Ini_B);	
	fwrite(buffer,sizeof(char),nbytes,pf);
		
 	nbytes = sprintf(buffer,"B-cell Max %i \n",AIS_p->Max_B);	
	fwrite(buffer,sizeof(char),nbytes,pf);

 	nbytes = sprintf(buffer,"B-cell Domain %i \n",AIS_p->B_domain);	
	fwrite(buffer,sizeof(char),nbytes,pf);

 	nbytes = sprintf(buffer,"B-cell Lifetime %i \n",AIS_p->B_lifetime);	
	fwrite(buffer,sizeof(char),nbytes,pf);

 	nbytes = sprintf(buffer,"Current M cells %i \n",AIS_p->n_M_cells);	
	fwrite(buffer,sizeof(char),nbytes,pf);
		
 	nbytes = sprintf(buffer,"N clones: %i \n",AIS_p->n_clones);	
	fwrite(buffer,sizeof(char),nbytes,pf);	
	
 	nbytes = sprintf(buffer,"N Ab random pool : %i \n",AIS_p->n_Ab);	
	fwrite(buffer,sizeof(char),nbytes,pf);

 	nbytes = sprintf(buffer,"Max M cells %i \n",AIS_p->Max_M);	
	fwrite(buffer,sizeof(char),nbytes,pf);


	for (k = 0; k < AIS_p->n_M_cells; k++){		
//		printf("fol*&^%^*&^%$%*&^%$");
		get_bitvector(AIS_p->M_cells[k], n, aux_bitvector); 
		nbytes = sprintf(buffer,"\n %s  -> %f",aux_bitvector, AIS_p->M_Defense[k]);  // Selected parameters in "bitvector" form
		fwrite(buffer,sizeof(char),nbytes,pf);	   
	}
 	nbytes = sprintf(buffer,"\nBeato Cell \n");	// Generation in "uint" form
	fwrite(buffer,sizeof(char),nbytes,pf);	
	get_bitvector(AIS_p->Beato_cell, n, aux_bitvector); 
	nbytes = sprintf(buffer,"\nSelected_vector: %s  -> %f",aux_bitvector, AIS_p->Beato_Defense);  // Selected parameters in "bitvector" form
	fwrite(buffer,sizeof(char),nbytes,pf);	  
    fclose(pf);
}

