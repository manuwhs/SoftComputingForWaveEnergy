
// Feature properties, this will calculate the result for every variable.


#include "../ELMheader.h"

int FP_op(FP_data *FP_p){
	int i,j;
	FP_p->order = (int *) malloc (sizeof(int)*FP_p->n_F);
	FP_p->Aux_Fitness = (double *) malloc (sizeof(double)*FP_p->n_F);	
	FP_p->F_Fitness = (double *) malloc (sizeof(double)*FP_p->n_F);
	FP_p->F_Pairs = (double **) malloc (sizeof(double*)*FP_p->n_F);
	
	FP_p->Aux_Selected = (unsigned int **) malloc (sizeof(unsigned int *)*FP_p->n_F);
	for (i = 0; i < FP_p->n_F; i++){
		FP_p->Aux_Selected[i] = (unsigned int *) malloc (sizeof(unsigned int)*n_var_chosen);
	}	
	
	
	for (i = 0; i < FP_p->n_F; i++){
		FP_p->F_Pairs[i] = (double *) malloc (sizeof(double)*FP_p->n_F);
	}
	FP_p->F_Affinity = (double **) malloc (sizeof(double*)*FP_p->n_F);
	for (i = 0; i < FP_p->n_F; i++){
		FP_p->F_Affinity[i] = (double *) malloc (sizeof(double)*FP_p->n_F);
	}
	FP_p->F_Selected = (unsigned int **) malloc (sizeof(unsigned int *)*FP_p->n_F);
	for (i = 0; i < FP_p->n_F; i++){
		FP_p->F_Selected[i] = (unsigned int *) malloc (sizeof(unsigned int)*n_var_chosen);
	}	
	
	if (FP_p->Operations & FS_Fitness) {
		FS_Evaluation_op(FP_p);
	}
	if (FP_p->Operations & FS_Affinity) {
		FS_Affinity_op(FP_p);
	}
	if (FP_p->Operations & SFS) {
		SFS_op(FP_p);
	}
	if (FP_p->Operations & SBE) {
		SBE_op(FP_p);
	}
	
	//-------> Free Memmory  <-------

}

int plot_FP_Fitness (FP_data * FP_p){
	int i;
	multiple_graph mgra;
	gsl_vector **graph_ES;
	int n_curves = 1;
	
	graph_ES = (gsl_vector **)malloc((n_curves + 1)*sizeof(gsl_vector *));
	for (i = 0; i < n_curves + 1; i++){
		graph_ES[i] = gsl_vector_alloc(FP_p->n_F);
	}
	for (i = 0; i < FP_p->n_F; i++){
		graph_ES[0]->data[i] = i;
	}	
	memcpy(graph_ES[1]->data, FP_p->F_Fitness, FP_p->n_F*sizeof(double));
	
	mgra.y = graph_ES;
	mgra.n_curves = n_curves;
	mgra.graph_name = (char *)malloc(30*sizeof(char));
	mgra.curve_names= (char **)malloc(n_curves*sizeof(char*));
	
	mgra.X_axis_n = (char *)malloc(30*sizeof(char));
	mgra.Y_axis_n = (char *)malloc(30*sizeof(char));
	
	for (i = 0; i < n_curves; i++){
		mgra.curve_names[i] = (char *)malloc(30*sizeof(char));
	}
	sprintf(mgra.curve_names[0], "Features");
	
	sprintf(mgra.graph_name, "FP");
	
	sprintf(mgra.X_axis_n , "Feature");
	sprintf(mgra.Y_axis_n , "Fitness");

	plot_multiple_graph(&mgra);
	print_Fitness(FP_p);
}

int print_Fitness(FP_data * FP_p){
	FILE *pf;
	int nbytes = 0;
	char buffer[150];
	char command[100];
	char result_folder[MAX_CHAR] = "Results";
	char subfolder[MAX_CHAR] = "FP";
	char data_file[MAX_CHAR] = "a_Feature_Fitness";
	char base_folder[300];
	char mkdir[MAX_CHAR] = "mkdir ";
	char ls[MAX_CHAR] = "ls";
	char folder_name[50];		// Name of the folder with the files
	
	char aux_bitvector[n];
	char date[60];
	char dir[60];
	time_t fechaActual ;
	struct tm * fechaPtr;
	
	int i,j, k;
	unsigned int * times_used;
	unsigned int * order;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	

	order = (unsigned int *) malloc (sizeof(unsigned int)*FP_p->n_F);
	order_double(FP_p->F_Fitness, order, FP_p->n_F);
	
	
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
	
	nbytes = sprintf(buffer,"Feature Fitness: %s \n", date);		
	fwrite(buffer,sizeof(char),nbytes,pf);
	
	for (i = 0; i < FP_p->n_F; i++){
		nbytes = sprintf(buffer,"%i -> %f \n", order[(FP_p->n_F - 1 ) -i], FP_p->F_Fitness[(FP_p->n_F - 1 ) - i]);
		fwrite(buffer,sizeof(char),nbytes,pf);
	}

    fclose(pf);		
	return 1;
}

int plot_FP_Affinity (FP_data * FP_p){
	int i;
	multiple_graph mgra;
	gsl_vector **graph_ES;
	int n_curves = 1;
	
	graph_ES = (gsl_vector **)malloc((n_curves + 1)*sizeof(gsl_vector *));
	for (i = 0; i < n_curves + 1; i++){
		graph_ES[i] = gsl_vector_alloc(FP_p->n_F);
	}
	for (i = 0; i < FP_p->n_F; i++){
		graph_ES[0]->data[i] = i;
	}	
	memcpy(graph_ES[1]->data, FP_p->F_Fitness, FP_p->n_F*sizeof(double));
	
	mgra.y = graph_ES;
	mgra.n_curves = n_curves;
	mgra.graph_name = (char *)malloc(30*sizeof(char));
	mgra.curve_names= (char **)malloc(n_curves*sizeof(char*));
	
	mgra.X_axis_n = (char *)malloc(30*sizeof(char));
	mgra.Y_axis_n = (char *)malloc(30*sizeof(char));
	
	for (i = 0; i < n_curves; i++){
		mgra.curve_names[i] = (char *)malloc(30*sizeof(char));
	}
	sprintf(mgra.curve_names[0], "Features");
	
	sprintf(mgra.graph_name, "FP");
	
	sprintf(mgra.X_axis_n , "Feature");
	sprintf(mgra.Y_axis_n , "Fitness");
	
	
	for (i = 0; i < FP_p->n_F; i++){ 
		sprintf(mgra.curve_names[0], "Features %i", i);
		sprintf(mgra.graph_name, "FP %i", i);	
		memcpy(graph_ES[1]->data, FP_p->F_Affinity[i], FP_p->n_F*sizeof(double));
		plot_multiple_graph(&mgra);
		print_Affinity(FP_p, i);
	}

}

int print_input_statistics (unsigned int ** input_v, int n, int Num){
	FILE *pf;
	int nbytes = 0;
	char buffer[150];
	char command[100];
	char result_folder[MAX_CHAR] = "Results";
	char subfolder[MAX_CHAR] = "FP";
	char data_file[MAX_CHAR] = "statistics";
	char base_folder[300];
	char mkdir[MAX_CHAR] = "mkdir ";
	char ls[MAX_CHAR] = "ls";
	char folder_name[50];		// Name of the folder with the files
	
	char aux_bitvector[n];
	char date[60];
	char dir[60];
	time_t fechaActual ;
	struct tm * fechaPtr;
	
	int i,j, k;
	unsigned int * times_used;
	unsigned int * order;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	

	times_used = (unsigned int *) malloc (sizeof(unsigned int)*n);
	order = (unsigned int *) malloc (sizeof(unsigned int)*n);
	for (i = 0; i < n; i++){
		times_used[i] = 0;
	}
	
	for (i = 0; i < Num; i++){
		for (j = 0; j < n; j++){	
			var_chosen = j/(sizeof(unsigned int)*8);
			bit_chosen = j%(sizeof(unsigned int)*8);			
			if ((input_v[i][var_chosen] & (1 << bit_chosen)) > 0){
				times_used[j]++;
			}
		}
	}
	
	ordenar_int (times_used, order, n);

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
	
	nbytes = sprintf(buffer,"Feature Statistics: %s \n", date);		
	fwrite(buffer,sizeof(char),nbytes,pf);
	
	for (i = 0; i < n; i++){
		nbytes = sprintf(buffer,"%i -> %i \n", order[i], times_used[i]);
		fwrite(buffer,sizeof(char),nbytes,pf);
	}

    fclose(pf);		
	return 1;
}

int print_Affinity(FP_data * FP_p, int num){
	FILE *pf;
	int nbytes = 0;
	char buffer[150];
	char command[100];
	char result_folder[MAX_CHAR] = "Results";
	char subfolder[MAX_CHAR] = "FP";
	char data_file[MAX_CHAR] = "F_Corr";
	char base_folder[300];
	char mkdir[MAX_CHAR] = "mkdir ";
	char ls[MAX_CHAR] = "ls";
	char folder_name[50];		// Name of the folder with the files
	
	char aux_bitvector[n];
	char date[60];
	char dir[60];
	time_t fechaActual ;
	struct tm * fechaPtr;
	
	int i,j, k;
	unsigned int * times_used;
	unsigned int * order;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	
	
	sprintf(data_file + strlen(data_file), "%i", num);	
	
	order = (unsigned int *) malloc (sizeof(unsigned int)*FP_p->n_F);
	order_double(FP_p->F_Affinity[num], order, FP_p->n_F);
	
	
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
	
	nbytes = sprintf(buffer,"Feature Affinity %i : %s \n",num, date);		
	fwrite(buffer,sizeof(char),nbytes,pf);
	
	for (i = 0; i < FP_p->n_F; i++){
		nbytes = sprintf(buffer,"%i ->  Affinity %f \n",order[i], FP_p->F_Affinity[num][i]);
		fwrite(buffer,sizeof(char),nbytes,pf);
	}

    fclose(pf);		
	return 1;
}

int SFS_op(FP_data *FP_p){
	int i;
	int F_left;		// Number of features left
	int num_moves;
	unsigned int * Selected;
	int num_iterations = 0;
	
	Selected = (unsigned int *) malloc (sizeof(unsigned int)*n_var_chosen);
	Selected[0] = 0;
	
	while (1) {
		num_moves = get_possible_inclusions(Selected, FP_p->n_F ,FP_p->F_Selected);
		
		if (num_moves == 0) {
			break;
		}
		threads_WorkOut(FP_p->F_Selected, FP_p->F_Fitness, FP_p->n_F, &Thread_p);	

		order_double (FP_p->F_Fitness, FP_p->order, num_moves);
		
		bitvector_cpy(Selected, FP_p->F_Selected[FP_p->order[num_moves-1]], n);
		bitvector_cpy(FP_p->Aux_Selected[num_iterations], Selected, n);
		FP_p->Aux_Fitness[num_iterations] =  FP_p->F_Fitness[num_moves-1];
		
		print_bitvector(Selected, n);
		
		num_iterations++;
	}
	print_vectors_data(FP_p->Aux_Selected, num_iterations, "SFS");
}

int SBE_op(FP_data *FP_p){
	int i;
	int F_left;		// Number of features left
	int num_moves;
	unsigned int * Selected;
	int num_iterations = 0;
	
	Selected = (unsigned int *) malloc (sizeof(unsigned int)*n_var_chosen);
	for (i = 0; i < FP_p->n_F; i++) {
		bitvector_setbit(Selected, i);
		
	}
	
	while (1) {
		num_moves = get_possible_exclusions(Selected, FP_p->n_F ,FP_p->F_Selected);
		
		if (num_moves == 1) {
			break;
		}
		threads_WorkOut(FP_p->F_Selected, FP_p->F_Fitness, FP_p->n_F, &Thread_p);	

		order_double (FP_p->F_Fitness, FP_p->order, num_moves);
		
		bitvector_cpy(Selected, FP_p->F_Selected[FP_p->order[num_moves-1]], n);
		bitvector_cpy(FP_p->Aux_Selected[num_iterations], Selected, n);
		FP_p->Aux_Fitness[num_iterations] =  FP_p->F_Fitness[num_moves-1];
		
		print_bitvector(Selected, n);
		
		num_iterations++;
	}
	printf("This is over \n");
	print_vectors_data(FP_p->Aux_Selected, num_iterations, "SBE");
	exit(0);
}

int get_possible_inclusions(unsigned int * VS, int n, unsigned int ** sol){
	int i;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	
	int u_bits = sizeof(unsigned int)*8;
	int	num_var = roof_int (n, u_bits);	
	
	int n_sol = 0;
	
	for (i = 0; i  < n ; i++){
		var_chosen = i/u_bits;
		bit_chosen = i%u_bits;
		if (((VS[var_chosen] >> bit_chosen) & 1 ) == 0){
			bitvector_cpy(sol[n_sol], VS, n);
			bitvector_setbit(sol[n_sol], i);
			n_sol++;
		}
	}
	return n_sol;
	
}

int get_possible_exclusions(unsigned int * VS, int n, unsigned int ** sol){
	int i;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	
	int u_bits = sizeof(unsigned int)*8;
	int	num_var = roof_int (n, u_bits);	
	
	int n_sol = 0;
	
	for (i = 0; i  < n ; i++){
		var_chosen = i/u_bits;
		bit_chosen = i%u_bits;
		if (((VS[var_chosen] >> bit_chosen) & 1 ) == 1){
			bitvector_cpy(sol[n_sol], VS, n);
			bitvector_clearbit(sol[n_sol], i);
			n_sol++;
		}
	}
	return n_sol;
	
}

int print_vectors_data(unsigned int ** SV, int num_Vectors, char *file_name){
	FILE *pf;
	int nbytes = 0;
	char buffer[150];
	char command[100];
	char result_folder[MAX_CHAR] = "Results";
	char subfolder[MAX_CHAR] = "FP";
	char data_file[MAX_CHAR];
	char base_folder[300];
	char mkdir[MAX_CHAR] = "mkdir ";
	char ls[MAX_CHAR] = "ls";
	char folder_name[50];		// Name of the folder with the files
	
	char aux_bitvector[n];
	char date[60];
	char dir[60];
	time_t fechaActual ;
	struct tm * fechaPtr;
	
	int i,j, k;
	unsigned int * times_used;
	unsigned int * order;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	

	double RMSE[num_Vectors];	// FOR DOING THE ELMs
	
	sprintf(data_file, "%s", file_name);
	
	ELM_p.FLAGS = PLOT_DATA_F;
	init_threads(&ELM_p, &Thread_p);
	threads_WorkOut(SV, RMSE, num_Vectors, &Thread_p);			
	ELM_p.FLAGS = 0;
	init_threads(&ELM_p, &Thread_p);
	times_used = (unsigned int *) malloc (sizeof(unsigned int)*n);
	order = (unsigned int *) malloc (sizeof(unsigned int)*n);
	for (i = 0; i < n; i++){
		times_used[i] = 0;
	}
	
	for (i = 0; i < num_Vectors; i++){
		for (j = 0; j < n; j++){	
			var_chosen = j/(sizeof(unsigned int)*8);
			bit_chosen = j%(sizeof(unsigned int)*8);			
			if ((SV[i][var_chosen] & (1 << bit_chosen)) > 0){
				times_used[j]++;
			}
		}
	}
	
	ordenar_int (times_used, order, n);

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
	
	
 	nbytes = sprintf(buffer,"Selection Vectors: \n");	// Generation in "uint" form
	fwrite(buffer,sizeof(char),nbytes,pf);	
	for (k = 0; k < num_Vectors; k++){		
		get_bitvector(SV[k], n, aux_bitvector); 
		nbytes = sprintf(buffer,"%s -> %f \n",aux_bitvector, RMSE[k]);  // Selected parameters in "bitvector" form
		fwrite(buffer,sizeof(char),nbytes,pf);	   

	}

	nbytes = sprintf(buffer,"Feature Statistics: %s \n", date);		
	fwrite(buffer,sizeof(char),nbytes,pf);
	
	for (i = 0; i < n; i++){
		nbytes = sprintf(buffer,"%i -> %i \n", order[i], times_used[i]);
		fwrite(buffer,sizeof(char),nbytes,pf);
	}

    fclose(pf);		
	return 1;
}

int FS_Evaluation_op(FP_data *FP_p){
	int i;
	//-------> Get the Fitness of every Feature <-------
	for (i = 0; i < FP_p->n_F; i++){
		zero_bitvector(FP_p->F_Selected[i], FP_p->n_F);
	}	
	for (i = 0; i < FP_p->n_F; i++) {
		bitvector_setbit(FP_p->F_Selected[i], i);
	}
	threads_WorkOut(FP_p->F_Selected, FP_p->F_Fitness, FP_p->n_F, &Thread_p );
	plot_FP_Fitness (FP_p);

}

int FS_Affinity_op(FP_data *FP_p){
	int i, j;
		//-------> Get the Fitness of every Feature <-------
	for (i = 0; i < FP_p->n_F; i++){
		zero_bitvector(FP_p->F_Selected[i], FP_p->n_F);
	}	
	for (i = 0; i < FP_p->n_F; i++) {
		bitvector_setbit(FP_p->F_Selected[i], i);
	}
	threads_WorkOut(FP_p->F_Selected, FP_p->F_Fitness, FP_p->n_F, &Thread_p );
	//-------> Get the Pairs <-------
	
	for (i = 0; i < FP_p->n_F; i++) {
		for (j = 0; j < FP_p->n_F; j++) {
			zero_bitvector(FP_p->F_Selected[j], FP_p->n_F);
			bitvector_setbit(FP_p->F_Selected[j], i);
			bitvector_setbit(FP_p->F_Selected[j], j);
		}
		threads_WorkOut(FP_p->F_Selected, FP_p->F_Pairs[i], FP_p->n_F, &Thread_p );
	}
	// If there is an improvement respect to the best (lowest) parameter, the Pairs will be positive
	// This value is expected to be near 0 and when very positive should be taken seriously.
	for (i = 0; i < FP_p->n_F; i++) {
		for (j = 0; j < FP_p->n_F; j++) {
			FP_p->F_Affinity[i][j] = min_d(FP_p->F_Fitness[i], FP_p->F_Fitness[j]) - FP_p->F_Pairs[i][j]; 
		}
	}

	plot_FP_Affinity (FP_p); 
}
