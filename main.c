#include"ELMheader.h"

/* USE:
	./ELM input_file -option options_file
	* 
	*/
	
int	n;			// Number of parameters
int n_var_chosen;	// Number of unsigned int variables needed for storing the parametrs
Thread_params Thread_p;
ELM_params ELM_p;

int main(int argc,char *argv[]){
	
	//  ----------------> TRAINING VARIABLES <----------------------
	double **x_train;	// Array with the training vectors
	gsl_vector **x_train_param;	// Array of gsl_vectors. One vector per parameter

	double *t_train;			// Array with the training results
	gsl_vector *y_train;		// gsl_vector with all the training results

	//  ----------------> TEST VARIABLES <----------------------
	double **x_test;			// Same as the training variables
	gsl_vector **x_test_param;
	
	double *t_test;
	gsl_vector *y_test;	

	// ---------------->  NORMALIZATION VARIABLES <------------------
	double *min_x_param;	// Array of the minimum value of each parameter n
	double *max_x_param;	// Array of the maximum value of each parameter n
	double t_min, t_max;		// Min and max value of the solution
	
	//---------------> BASIC PARAM VARIABLES <------------------
	char ** operation_config;		// Array of string where the params will be allocated
	Input_params datap;	
	int ELM_Operation;	// Operation we will do:
						// - TIME and ERMS test
						// - Plain ELMs
						// - GA algorithm
						// - ES algorithm
	int Nh = 0;			// Number of hidden neurons
	int FLAGS = 0;
	//--------------> Timing variables <----------------- //
	struct timeval time_start,time_end;		// Time variables for general use
	float time_passed;

	//--------------> General variables <----------------- //
	int i;
	int errores;
		
	//--------------> Graph Variables <----------------- //
	linear_graph l_graph;		// For linear graph
	linear_graph l_graph2;		// For linear graph
	unsigned int *Graph_selected;

	//--------------> Multiple ELM variables <----------------- //
	ELMs_data ELMs_p;
	
	//--------------> Exahustive search variables <----------------- //
	ES_data ES_p;
	
	// ----------------> Genetic algorithm variables <--------------------
	GA_data		GA_p;
	
	// ----------------> Simmulated Anneling variables <--------------------
	SA_data		SA_p;
	
	// ----------------> Artificial Inmmune System Variables <--------------------
	AIS_data AIS_p;
	
	// ----------------> Feature Selection <--------------------
	FP_data FP_p;
	
// ****************************************************************************************
// ***************** PROCESS THE CONFIG FILE AND GET PARAMMETERS **************************
// ****************************************************************************************

	errores = gettimeofday( &time_start, NULL);
	srand(time_start.tv_usec);		// START RANDOM
	if (argc < 4){
		printf("Not enought parameters \n");
		exit(-1);
	}
	
	// Basic config parameters
	get_ELM_params(argv[1], &datap);
	print_ELM_values(&datap);
	n = datap.n_input;
	n_var_chosen = roof_int (n, sizeof(unsigned int)*8);
	Nh = datap.Nh;
	
	//-----------> EXAUSTIVE SEARCH <--------------
	if (strcmp(argv[2],"-ES") == 0){
printf("ES operation selected \n");
		ELM_Operation = EXHAUSTIVE_SEARCH;
		get_ES_params (argv[3], &ES_p);
	}
	//-----------> SIMPLE ELMs <--------------
	else if (strcmp(argv[2],"-ELM") == 0){
printf("ELM operation selected \n");
		ELM_Operation = DO_ELMS ;
		get_ELMs_params (argv[3], &ELMs_p);
		FLAGS = ELMs_p.FLAGS;
	}
	//-----------> GENETIC ALGORITHM <--------------
	else if (strcmp(argv[2],"-GA") == 0){
printf("GA operation selected \n");
		ELM_Operation = GA_ALGORITHM;
		get_GA_params (argv[3], &GA_p);
	}
	//-----------> Artificial Inmune Sytem <--------------	
	else if (strcmp(argv[2],"-AIS") == 0){
printf("AIS operation selected \n");
		ELM_Operation = AIS_ALGORITHM;
		get_AIS_params (argv[3], &AIS_p);
	}
	//-----------> SIMMULATED ANNELING <--------------
	else if (strcmp(argv[2],"-SA") == 0){
printf("SA operation selected \n");
		ELM_Operation = SA_ALGORITHM;
		get_SA_params (argv[3], &SA_p);
		
	}
	//-----------> GRAPHING <--------------
	else if (strcmp(argv[2],"-Graph") == 0){
printf("Graph operation selected \n");
		ELM_Operation = TIME_ERMS_GRAPH;
		Graph_selected = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen); // No podemos modificarlo desde la funcion ahora
		get_Graph_params (argv[3],&l_graph, &l_graph2,Graph_selected);

	}
	else if (strcmp(argv[2],"-FP") == 0){
printf("FP operation selected \n");
		ELM_Operation = FP_ALGORITHM;
		get_FP_params (argv[3], &FP_p);
	}
	else {
		printf("Opcion no valida \n");
		exit(-1);
	}

// ***************** READ AND BUILD THE GENERAL TRAINING AND TESTING MATRIX  ***************

	printf("Loading data.. ");
	errores = gettimeofday( &time_start, NULL);
	if (errores == -1){
		printf("No hemos podido saber la hora\n");
	}


	//load_data3(&x_train, &t_train, n, datap.Ntrain, datap.train_input_dir);
	//load_data3(&x_test, &t_test, n, datap.Ntest, datap.test_input_dir);

	Super_LoadX(&x_train, &t_train, &x_test, &t_test, &datap) ;
	shuffle_train_test(x_train, t_train, x_test, t_test, &datap );
	
//	print_matrix(x_train,datap.n_input,datap.Ntrain);
//	print_array(t_train,datap.Ntrain);
//	print_matrix(x_test,datap.n_input,datap.Ntest);
//	print_array(t_test,datap.Ntest);
	
	errores = gettimeofday( &time_end, NULL);
	time_passed = get_time_passed(time_start,time_end);
	printf(" %f s\n",time_passed);
	
// *********************************** NORMALIZE THE DATA *********************************
	min_x_param = (double *)malloc(datap.Ntrain*sizeof(double));		// Reserve memmory		
	max_x_param = (double *)malloc(datap.Ntrain*sizeof(double));
	
	get_matrix_row_minmax(x_train,min_x_param,max_x_param,n,datap.Ntrain); // Get normalization values
	get_array_minmax(t_train,datap.Ntrain, &t_min, &t_max);
	
	normalize_rows_matrix(x_train,min_x_param,max_x_param,n,datap.Ntrain);	// Normalize training
	normalize_array(t_train,  datap.Ntrain, t_min, t_max);

	normalize_rows_matrix(x_test ,min_x_param,max_x_param,n,datap.Ntest);	// Normalize testing with same norm
	normalize_array(t_test,  datap.Ntest, t_min, t_max);
	
// ***************** TRANSFORM ARRAY INTO GSL STRUCTURES ******************************	
	convert_data_to_gsl(x_train, t_train, n, datap.Ntrain, &x_train_param, &y_train);
	convert_data_to_gsl(x_test, t_test, n, datap.Ntest, &x_test_param, &y_test);	

//	print_gsl_vector(y_train);
//	print_gsl_vector(x_train_param[0]);
	//exit(0);
	
// ****************************************************************************************
// ******************************* INIT ELM DATA AND POSIX THREADS ************************
// ****************************************************************************************
	// Init threads
	set_up_threads( datap.n_threads, &Thread_p);
	
	ELM_p.Nh = 			Nh;
	ELM_p.activation_f = 	datap.activation_f;
	ELM_p.n_ELMs = 	datap.n_ELMs;
	ELM_p.FLAGS = 		FLAGS; 
	ELM_p.x_train_param = x_train_param;
	ELM_p.x_test_param =  x_test_param;
	ELM_p.y_train = 		y_train;
	ELM_p.y_test = 		y_test  ;
	ELM_p.t_max = 		t_max ;
	ELM_p.t_min = 		t_min;

	init_threads(&ELM_p, &Thread_p);
// ************************************************************
// ********************* DO OPERATION *************************
// ************************************************************

	if (ELM_Operation == TIME_ERMS_GRAPH) {
	  print_bitvector(Graph_selected, n);
	  Graph_op (&l_graph, &l_graph2, Graph_selected);
	  free(Graph_selected);
	}
	else if (ELM_Operation == EXHAUSTIVE_SEARCH){
	  ES_op (&ES_p);
	}
	else if (ELM_Operation == DO_ELMS){
		ELMs_op (&ELMs_p);
	}
	else if (ELM_Operation == GA_ALGORITHM) {	
		GA_op(&GA_p, PLOT_DATA);
	}
	else if (ELM_Operation == SA_ALGORITHM) {	
		SA_op(&SA_p, PLOT_DATA);
	}
	else if (ELM_Operation == AIS_ALGORITHM) {	
		AIS_op(&AIS_p, PLOT_DATA);
	}
	else if (ELM_Operation == FP_ALGORITHM) {	
		FP_op(&FP_p);
	}
// ************************************************************
// *****************  FREE MEMMORY  ************************
// ************************************************************

	for (i = 0; i < n; i++){
		free(x_train[i]);
		free(x_test[i]);
	}
	free (x_train);
	free (x_test);
	free(t_train);
	free(t_test);
	free(min_x_param);
	free(max_x_param);
	

	for (i = 0; i < n; i++){
		gsl_vector_free(x_train_param[i]);
		gsl_vector_free(x_test_param[i]);
	}
	free(x_train_param);
	free(x_test_param);

	gsl_vector_free (y_train);	
	gsl_vector_free (y_test);	

	destroy_threads(&Thread_p);
	return 0;
}


