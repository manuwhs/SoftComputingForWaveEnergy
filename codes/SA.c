#include "../ELMheader.h"
// Simmulad anneling

int SA_op(SA_data *SA_p, int FLAGS){
	unsigned int i,j;
	
	//------------> Variables for statistics <------------
	double *	Gen_best_ERMS;	// Array with the average ERMS of the best ingots so far
	double *	Gen_ave_ERMS;	// Array with the average ERMS  of the ingots

	multiple_graph mgra;
	gsl_vector **graph_SA;

	if (FLAGS & PLOT_DATA){
	//--------> Reserve memmory for the generation statistics <------------
		Gen_best_ERMS = (double *) malloc (sizeof(double)*SA_p->stop_cond);
		Gen_ave_ERMS = (double *) malloc (sizeof(double)*SA_p->stop_cond);	
		
	}	
	
	//--------> Reserve memmory for the operations <------------
	SA_p->Ingots = (unsigned int **) malloc (sizeof(unsigned int *)*SA_p->n_ingots);
	SA_p->New_Ingots = (unsigned int **) malloc (sizeof(unsigned int *)*SA_p->n_ingots);
	for (i = 0; i < SA_p->n_ingots; i++){
		SA_p->Ingots[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
		SA_p->New_Ingots[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
	}
	
	SA_p->Energy = (double *) malloc (sizeof(double)*SA_p->n_ingots);
	SA_p->New_Energy = (double *) malloc (sizeof(double)*SA_p->n_ingots);
	
	SA_p->Top_Ingots = (unsigned int **) malloc (sizeof(unsigned int *)*SA_p->n_ingots);
	for (i = 0; i < SA_p->n_ingots; i++){
		SA_p->Top_Ingots[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
	}	
	SA_p->Top_Energy = 	(double *) malloc (sizeof(double)*SA_p->n_ingots);			


	Init_SA_Ingots(SA_p);

	get_Energy(SA_p);			// GET THE Initial energy
	
	for (i = 0; i < SA_p->n_ingots; i++){ // Initialice Top ingots
		SA_p->Top_Energy[i] = SA_p->Energy[i];
	}		
	
	if (SA_p->temp <= 0.0){		// Do the autotunning
		auto_tune_k(SA_p);
	}
	print_SA_values(SA_p);
	
//--------> Start Cooling process <------------	
	for (i = 0; i < SA_p->stop_cond; i++){
		printf("Freezing %i, Temperature %f \n", i, SA_p->temp);

		freeze_ingots (SA_p);		// Mutate ingots	
		get_New_Energy(SA_p);		// Get the new energy
		apply_SA(SA_p);				// Choose weather or not we keep them
		
	//	print_Ingots(SA_p);
		
		//------> Get values for the statistics <------
		if (FLAGS & PLOT_DATA){
			
			Gen_best_ERMS[SA_p->n_freezings] = 0;
			Gen_ave_ERMS[SA_p->n_freezings] = 0;		
		
			for (j = 0; j < SA_p->n_ingots; j++) {
				Gen_ave_ERMS[SA_p->n_freezings] +=  SA_p->Energy[j];
			}
			Gen_ave_ERMS[SA_p->n_freezings]/= SA_p->n_ingots;
			
			for (j = 0; j < SA_p->n_ingots; j++) {
				Gen_best_ERMS[SA_p->n_freezings] +=  SA_p->Top_Energy[j];
			}
			Gen_best_ERMS[SA_p->n_freezings]/= SA_p->n_ingots;
		}		
		SA_p->n_freezings++;
		
		if (check_SA_convergence (SA_p, Gen_ave_ERMS, Gen_best_ERMS)){
			break;
		}
	}
printf("*************************************************** \n");	
printf("************ FREEZING DONE *********************** \n");
printf("*************************************************** \n");	

//----------> Treat the best solutions found <-----------

	print_SA_data(SA_p);
	print_input_statistics (SA_p->Top_Ingots, n, SA_p->n_ingots);
	print_SA_stuff (SA_p);
	
	
	
	if (FLAGS & PLOT_DATA){
		//--------------> Plot data <---------------
		graph_SA = (gsl_vector **)malloc(3*sizeof(gsl_vector *));
		for (i = 0; i < 3; i++){
			graph_SA[i] = gsl_vector_alloc(SA_p->n_freezings);
		}
		for (i = 0; i < SA_p->n_freezings; i++){
			graph_SA[0]->data[i] = i;
		}	
		memcpy(graph_SA[1]->data, Gen_best_ERMS, sizeof(double)*SA_p->n_freezings);
		memcpy(graph_SA[2]->data,Gen_ave_ERMS, sizeof(double)*SA_p->n_freezings);
		
		mgra.y = graph_SA;
		mgra.n_curves = 2;
		mgra.graph_name = (char *)malloc(30*sizeof(char));
		mgra.curve_names= (char **)malloc(2*sizeof(char*));
		
		mgra.X_axis_n = (char *)malloc(30*sizeof(char));
		mgra.Y_axis_n = (char *)malloc(30*sizeof(char));
		
		for (i = 0; i < 2; i++){
			mgra.curve_names[i] = (char *)malloc(30*sizeof(char));
		}
		sprintf(mgra.curve_names[0], "Best Ingots average");
		sprintf(mgra.curve_names[1], "Ingots average");
		
		sprintf(mgra.graph_name, "SA");
		
		sprintf(mgra.X_axis_n , "Freezing");
		sprintf(mgra.Y_axis_n , "Energy (RMSE)");

		plot_multiple_graph(&mgra);
		
		//--------------> Free Graph data <---------------
		free(Gen_best_ERMS);
		free(Gen_ave_ERMS);	
		
		for (i = 0; i < 3; i++){
			gsl_vector_free(graph_SA[i]);
		}
		free(graph_SA);
		 
		free(mgra.X_axis_n);
		free(mgra.Y_axis_n);
		free(mgra.graph_name);
		for (i = 0; i < 2; i++){
			free(mgra.curve_names[i]);
		}
		free(mgra.curve_names);
	}
	
	//--------------> Free SA data <---------------
	for (i = 0; i < SA_p->n_ingots; i++){
		free(SA_p->Ingots[i]);
		free(SA_p->New_Ingots[i]);
	}
	free(SA_p->Ingots);
	free(SA_p->New_Ingots);

	free(SA_p->Energy);					// Energy of the current blocks
	free(SA_p->New_Energy);				// Energy of the next colded blocks
		
	for (i = 0; i < SA_p->n_ingots; i++){
		free(SA_p->Top_Ingots[i]);
	}
	free(SA_p->Top_Ingots);
	free(SA_p->Top_Energy);				// Energy of the ingots
	
	if (SA_p->n_Selected > 0){
		for (i = 0; i < SA_p->n_Selected; i++){
			free(SA_p->Selected_Ingots[i]);
		}
		free(SA_p->Selected_Ingots);		
	}
	return 1;
}

int get_Energy(SA_data *SA_p){

	threads_WorkOut(SA_p->Ingots,SA_p->Energy, SA_p->n_ingots, &Thread_p);
	return 1;
}

int get_New_Energy(SA_data *SA_p){

	threads_WorkOut(SA_p->New_Ingots, SA_p->New_Energy, SA_p->n_ingots, &Thread_p);
	return 1;
}

int freeze_ingots (SA_data *SA_p){
	unsigned int i;
	for (i = 0; i < SA_p->n_ingots; i++){
		bitvector_cpy(SA_p->New_Ingots[i], SA_p->Ingots[i], SA_p->n_ingots);
		mutate_bitvector(SA_p->New_Ingots[i], SA_p->num_bits, 1 + rand()% SA_p->max_mut);
	}
	
	return 1;
}

int apply_SA(SA_data *SA_p){
	int i;
	int flag;
	double delta;
//	printf("Antes \n");
//	print_Ingots(SA_p);
	for (i = 0; i < SA_p->n_ingots; i++){  	// Apply SA
		delta = SA_p->Energy[i] -  SA_p->New_Energy[i];
		if (delta >= 0){		// If its a better solution we swap
			bitvector_cpy(SA_p->Ingots[i], SA_p->New_Ingots[i], SA_p->num_bits);
			SA_p->Energy[i] = SA_p->New_Energy[i];
		}
		else if (exp(delta/SA_p->temp) > (float)rand ()/RAND_MAX ){				// Else --> SA
			bitvector_cpy(SA_p->Ingots[i], SA_p->New_Ingots[i], SA_p->num_bits);
			SA_p->Energy[i] = SA_p->New_Energy[i];		
		}
	}
	//printf("Despues \n");
	//print_Ingots(SA_p);
	for (i = 0; i < SA_p->n_ingots; i++){  	
		if (SA_p->Energy[i] < SA_p->Top_Energy[i]){	// If its better solution
			bitvector_cpy(SA_p->Top_Ingots[i], SA_p->Ingots[i], SA_p->num_bits);
			SA_p->Top_Energy[i] = SA_p->Energy[i];	
		}		
	}

	SA_p->temp*= SA_p->k;	// Cool temperature
	return 1;
}

int auto_tune_k(SA_data *SA_p){
	/* For calculating the k we will run the algortihm once to see what is the average delta
	 * (change in energy) that happens between iterations.
	 */
	unsigned int i,j;
	double delta = 0;
	int n = 2;
	
	for (i = 0; i < n; i++) {
		freeze_ingots (SA_p);		// Mutate ingots	
		get_New_Energy(SA_p);		// Get the new energy
		
		for (j = 0; j < SA_p->n_ingots; j++){  	
			delta += abs_d(SA_p->Energy[j] -  SA_p->New_Energy[j]);
		} 
	}
	delta/= n*SA_p->n_ingots;
	SA_p->temp = delta;	// Cool temperature
	return 1;
}

int Init_SA_Ingots(SA_data *SA_p){
	int i,j;
	int num_min;
	int rest;
	
	if (SA_p->n_Selected == 0) {		// Select Random generation 
		for (i = 0; i < SA_p->n_ingots; i++){
			zero_bitvector(SA_p->Ingots[i], SA_p->num_bits);
			get_random_bitarray(SA_p->Ingots[i], SA_p->num_bits);			// Initialize values of population
		}	
	}
	else {		// Put the given vector as generation values sappin all the ingots
		num_min = SA_p->n_ingots / SA_p->n_Selected;
		rest = SA_p->n_ingots % SA_p->n_Selected;
		
		for (i = 0; i < SA_p->n_Selected; i++){
			for (j = 0; j < num_min; j++){
				bitvector_cpy(SA_p->Ingots[i*(num_min)+j], SA_p->Selected_Ingots[i], SA_p->num_bits);
			}
		}
		for (i = 0; i < rest; i++){
			bitvector_cpy(SA_p->Ingots[SA_p->n_ingots - (rest) + i], SA_p->Selected_Ingots[i], SA_p->num_bits);
		}	
	}
}

int print_SA_values(SA_data *SA_p){
		printf("********** Simulated Anneling values ***************** \n");
		printf("N Ingots: %i \n",SA_p->n_ingots);		// Number of parallel tries Lingotes
		printf("N freezings: %i \n",SA_p->n_freezings);	// Number of coldings.
		printf("N bits: %i \n",SA_p->num_bits);		// Number of bits of the population	
		printf("Stop cond: %i \n",SA_p->stop_cond);	// Number of coldings.
				
		printf("Temperature: %f \n",SA_p->temp);		// Simulated anneling temperature
		printf("Decreasing constant  : %f \n",SA_p->k);			// Simulated constant  temp[t+1] = temp[t]*k

		printf("Max mutations: %i \n",SA_p->max_mut);			// Number of maximum mutations per freezing
		printf("N Selected  : %i \n",SA_p->n_Selected);		// Number of selected ingots for initial population
		printf("********** ****************************************** \n");
		return 1;
}

int print_Ingots(SA_data *SA_p){
	unsigned int j;
	printf("********** Ingots ***************** \n");
	for (j = 0; j < SA_p->n_ingots; j++) {
		print_bitvector(SA_p->Ingots[j],SA_p->num_bits);	
		printf("Energy: %f \n", SA_p->Energy[j]);
	}
	printf("********** Best Ingots ***************** \n");
	for (j = 0; j < SA_p->n_ingots; j++) {
		print_bitvector(SA_p->Top_Ingots[j],SA_p->num_bits);	
		printf("Top Energy: %f \n", SA_p->Top_Energy[j]);
	}
	printf("*********************************** \n");
}
	
int print_SA_data(SA_data *SA_p){
		ELM_p.FLAGS = PLOT_DATA_F;
		init_threads(&ELM_p, &Thread_p);
		threads_WorkOut(SA_p->Top_Ingots, SA_p->Top_Energy, SA_p->n_ingots, &Thread_p);		
		return 1;
}

void print_SA_stuff (SA_data *SA_p){
 	FILE *pf;
	int nbytes = 0;
	char buffer[150];
	char command[100];
	char result_folder[MAX_CHAR] = "Results";
	char subfolder[MAX_CHAR] = "SA";
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
  
 	nbytes = sprintf(buffer,"N_ingots: %u \n", SA_p->n_ingots);	
	fwrite(buffer,sizeof(char),nbytes,pf);
	
 	nbytes = sprintf(buffer,"N_freezings: %u \n", SA_p->n_freezings);	
	fwrite(buffer,sizeof(char),nbytes,pf);
		
 	nbytes = sprintf(buffer,"Temperature: %f \n", SA_p->temp);	
	fwrite(buffer,sizeof(char),nbytes,pf);

 	nbytes = sprintf(buffer,"Constant: %f \n", SA_p->k);	
	fwrite(buffer,sizeof(char),nbytes,pf);
		
 	nbytes = sprintf(buffer,"Bits of the population %i\n\n", SA_p->num_bits);	
	fwrite(buffer,sizeof(char),nbytes,pf);	

	for (k = 0; k <SA_p->n_ingots; k++){		
		get_bitvector(SA_p->Top_Ingots[k], n, aux_bitvector); 
		nbytes = sprintf(buffer,"%s -> %f \n",aux_bitvector,SA_p->Top_Energy[k]);  // Selected parameters in "bitvector" form
		fwrite(buffer,sizeof(char),nbytes,pf);	   

	}

    fclose(pf);
}

int check_SA_convergence (SA_data *SA_p, double *ener, double *best_ener ) {
	int n = 20;
	int i;
	double ave = 0;
	double top_ave = 0;
	double diff;
	int flag;
	
	if (SA_p->n_freezings < n) {
		return 0;
	}
	/*
	for (i = 0; i < n; i++) {
		ave += ener[(SA_p->n_freezings - 1) - i];
		top_ave += best_ener[(SA_p->n_freezings - 1) - i];
	}
	
	diff = 2*(ave - top_ave)/(ave + top_ave);
	
	if ( diff < 0.1 ) {	// It has converged
		printf("Ended coz of convergence with the best \n");
		return 1;
	}
	
	flag = 0;
	for (i = 0; i < n - 1; i++) {
		if (ener[(SA_p->n_freezings-1)-i] < ener[(SA_p->n_freezings-2)-i]){	// Normaliced
			flag = 1;
			break;
		}
	}
	if ( flag == 0 ) {	// It has converged
		printf("Ended coz of variance \n");
		return 1;
	}
	*/
	return 0;
}

