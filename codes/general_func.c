#include "../ELMheader.h"

gsl_matrix* get_selected_input_matrix(gsl_vector **param, unsigned int *chosen_param, int n_chosen){
	int i, aux_i;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable
	gsl_matrix *final_matrix;		// Returned, not freed.
	// Allocate memory for the matrix
	final_matrix = gsl_matrix_alloc(param[0]->size,n_chosen);	// Use the size of the first vector for convenience
																// Note that the actual number of elements could be
																// greater than this value. This way is a fast methof
																// for using a changing number of TRAINING VECTORS samples
	aux_i= 0;
	i = 0;
	// Put the vectors on the matrix
	while(aux_i < n_chosen){		// For every parammeter selected
		var_chosen = i/(sizeof(int)*8);
		bit_chosen = i%(sizeof(int)*8);

		if((chosen_param[var_chosen] & (1 << bit_chosen)) > 0){	// If its a choosen param
			param[i]->size = param[0]->size; // ******* SET THE SIZE TO THE GIVEN BY THE FIRST VECTOR ****
			//*** THIS IS NOT NECESARY IF WE NEVER CHANGE THE LENGTH OF THE VECTORS ***
			
			
			gsl_matrix_set_col (final_matrix, aux_i, param[i]);
			// Decrement the number of chosen
			aux_i++;
		}
		i++;
	}
	
	return final_matrix;
	
}
// Pasamos X y t por referencia

void convert_data_to_gsl(double **X, double *t, int n, int Nsamples, gsl_vector ***Xgsl, gsl_vector **tgsl){
	int i;
	*Xgsl = (gsl_vector **)malloc(n*sizeof(gsl_vector *));
	if (*Xgsl == NULL){
		perror("Malloc");
	}	
	for (i = 0; i < n; i++){						// For every parammeter of the vector.
			(*Xgsl)[i] =(gsl_vector *) gsl_vector_alloc(Nsamples);		// Reserve memmory
			
			//*******  PORQUE HAY QUE HACER LO DE LOS 2 !!!!!!!!!!!! **********  ?????????????????
			
			memcpy((void *)((*Xgsl)[i]->data), &(X[i][0]),(Nsamples)*sizeof(double));	// Copy values to vector
			
	}
//print_array(x_train[0], Ntrain);
//print_gsl_vector(x_train_param[0]);
	
	*tgsl = gsl_vector_alloc(Nsamples);		// Reserve memmory	
	memcpy((*tgsl)->data, t,Nsamples*sizeof(double));	// Copy values to vector	
	
}

float get_time_passed(struct timeval time_start,struct timeval time_end){
	int dif_sec, dif_usec;
	float result = 0;
	
	dif_sec = time_end.tv_sec - time_start.tv_sec;
	dif_usec = time_end.tv_usec - time_start.tv_usec;
	if (dif_sec == 0) {
		result = (float)dif_usec/1000000;
	} else {
		result += dif_sec;
		result += ((float)dif_usec)/1000000;
	}
	return result;
}

double math_combination(int n, int r){
	// nCr = n! /(n - r)! 
    double cnm = 1.0;
    int aux = n-r;
    
    int i;

    for ( i = 0; i < n ; i++){
		cnm *= n;
		n--;
    }
    r = aux;
    for ( i = 0; i < aux ; i++){
		cnm /= r;
		r--;
    }
    return cnm;
}

void load_data(double ***X, double **t, int n, int Nsamples, char *Xdir, char *tdir){
	int i, j;
	FILE *pf;
//******** Allocate memmory for the general matrix *********** 
	*X = (double **) malloc(n*sizeof(double*));
	if (*X == NULL){
		perror("Malloc");
	}
	for (i = 0; i < n; i++){
		(*X)[i] = (double *)malloc(Nsamples*sizeof(double));
		if ((*X)[i] == NULL){
			perror("Malloc");
		}
	}
	
	//printf("Memoria para la matriz basica principal reservada\n");
	//printf("Dimensiones %i x %i \n",n,Nsamples); 
	
	// Read the vectors float by float accoding to the format spacified.
	pf = fopen(Xdir,"r");	// Open file "read" mode
	if (pf == NULL){
		perror("fopen");
	}
	
	for (i = 0; i < Nsamples; i++){	// For every samples vector
		for (j = 0; j < n; j++){	// For every parammeter of the vector.
			fscanf(pf,"%lf",&((*X)[j][i]));
		}
	}
	fclose(pf);


	//------------>Read and allocate memmory for the result vector <---------------/
	
	*t = (double *) malloc(Nsamples*sizeof(double));
	if (t == NULL){
		perror("Malloc");
	}	
	pf = fopen(tdir,"r");	// Open file "read" mode
	if (pf == NULL){
		perror("fopen");
	}
	for (i = 0; i < Nsamples; i++){	// For every training vector
			fscanf(pf,"%lf",&((*t)[i]));
	}
	fclose(pf);		
	// printf(" \nTrain datos leidos \n");
}

/* This function receives an input file with the following format !!
 * [Año Mes Día Hora Validez ALTURA Parámetro1 ... Parámetro14].
 * It will place the ALTURA into the double **t structure and the Parametro1-14 into the double ***X,
 * n must be 14 from the outside. 
 */

void Super_Load(double ***X_tr, double **t_tr, double ***X_te, double **t_te, Input_params *Input_p) {
	int i, j;
	FILE *pf;
	
//******** Allocate memmory *********** 
//  X_tr[param][muestras]

	*X_tr = (double **) malloc(Input_p->n_input*sizeof(double*));
	if (*X_tr == NULL){
		perror("Malloc");
	}
	for (i = 0; i < n; i++){
		(*X_tr)[i] = (double *)malloc(Input_p->Ntrain*sizeof(double));
		if ((*X_tr)[i] == NULL){
			perror("Malloc");
		}
	}
	*t_tr = (double *) malloc(Input_p->Ntrain*sizeof(double));
	if (t_tr == NULL){
		perror("Malloc");
	}	
	
	*X_te = (double **) malloc(Input_p->n_input*sizeof(double*));
	if (*X_te == NULL){
		perror("Malloc");
	}
	for (i = 0; i < n; i++){
		(*X_te)[i] = (double *)malloc(Input_p->Ntest*sizeof(double));
		if ((*X_te)[i] == NULL){
			perror("Malloc");
		}
	}
	*t_te = (double *) malloc(Input_p->Ntest*sizeof(double));
	if (t_te == NULL){
		perror("Malloc");
	}	
	
//printf("Memoria para la matriz basica principal reservada\n");
//printf("Dimensiones %i x %i \n",n,Nsamples); 
	
// Read the vectors float by float accoding to the format spacified.

	if (Input_p->test_input_dir[0] != 0 ) {	// If we have specified separete files for train and test
		// Read training Input vectors 
		pf = fopen(Input_p->train_input_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen");
		}
		for (i = 0; i < Input_p->Ntrain; i++){	// For every samples vector
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_tr)[j][i]));
			}
		}
		fclose(pf);
		
		// Read training Output vectors
		pf = fopen(Input_p->train_output_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen");
		}
		for (i = 0; i < Input_p->Ntrain; i++){	// For every samples vecto
			fscanf(pf,"%lf",&((*t_tr)[i]));
		}
		fclose(pf);	
		
		// Read testing Input vectors 
		pf = fopen(Input_p->test_input_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen");
		}
		for (i = 0; i < Input_p->Ntest; i++){	// For every samples vector
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_te)[j][i]));
			}
		}
		fclose(pf);
		
		// Read testing Output vectors
		pf = fopen(Input_p->test_output_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen");
		}
		for (i = 0; i < Input_p->Ntrain; i++){	// For every samples vector
			fscanf(pf,"%lf",&((*t_te)[i]));
		}
		fclose(pf);		
	}
	else {									// If we havent specified separate files
			// Read training Input vectors 
		pf = fopen(Input_p->train_input_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen");
		}
		for (i = 0; i < Input_p->Ntrain; i++){	// For every samples vector
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_tr)[j][i]));
			}
		}
		for (i = 0; i < Input_p->Ntest; i++){	// For every samples vector
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_te)[j][i]));
			}
		}	
		fclose(pf);
		
		// Read training Output vectors
		pf = fopen(Input_p->train_output_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen");
		}
		for (i = 0; i < Input_p->Ntrain; i++){	// For every samples vector
			fscanf(pf,"%lf",&((*t_tr)[i]));
		}
		for (i = 0; i < Input_p->Ntest; i++){	// For every samples vector
			fscanf(pf,"%lf",&((*t_te)[i]));
		}
		fclose(pf);	
	}
}


void load_data2(double ***X, double **t, int n, int Nsamples, char *dir){
	int i, j;
	FILE *pf;
	double aux_d;
//******** Allocate memmory *********** 
	*X = (double **) malloc(n*sizeof(double*));
	if (*X == NULL){
		perror("Malloc");
	}
	
	*t = (double *) malloc(Nsamples*sizeof(double));
	if (t == NULL){
		perror("Malloc");
	}	
	
	for (i = 0; i < n; i++){
		(*X)[i] = (double *)malloc(Nsamples*sizeof(double));
		if ((*X)[i] == NULL){
			perror("Malloc");
		}
	}
	
//printf("Memoria para la matriz basica principal reservada\n");
//printf("Dimensiones %i x %i \n",n,Nsamples); 
	
	// Read the vectors float by float accoding to the format spacified.
	pf = fopen(dir,"r");	// Open file "read" mode
	if (pf == NULL){
		perror("fopen");
	}
	
	for (i = 0; i < Nsamples; i++){	// For every samples vector
		fscanf(pf,"%lf",&aux_d);	// Año 
		fscanf(pf,"%lf",&aux_d);	// Mes 
		fscanf(pf,"%lf",&aux_d);	// Día 
		fscanf(pf,"%lf",&aux_d);	// Hora 
		fscanf(pf,"%lf",&aux_d);	// Validez
		fscanf(pf,"%lf",&((*t)[i]));	// Get the ALTURA into the output array
		for (j = 0; j < n; j++){	// For every parammeter of the vector.
			fscanf(pf,"%lf",&((*X)[j][i]));
		}
		if (aux_d <0) {		// aux_d = -1  -> Dato no valido
			i--;
		}
	}
	fclose(pf);
}

void load_data3(double ***X, double **t, int n, int Nsamples, char *dir){
	int i, j;
	FILE *pf;
	double aux_d;
//******** Allocate memmory *********** 
	*X = (double **) malloc(n*sizeof(double*));
	if (*X == NULL){
		perror("Malloc");
	}
	
	*t = (double *) malloc(Nsamples*sizeof(double));
	if (t == NULL){
		perror("Malloc");
	}	
	
	for (i = 0; i < n; i++){
		(*X)[i] = (double *)malloc(Nsamples*sizeof(double));
		if ((*X)[i] == NULL){
			perror("Malloc");
		}
	}
	
//printf("Memoria para la matriz basica principal reservada\n");
//printf("Dimensiones %i x %i \n",n,Nsamples); 
	
	// Read the vectors float by float accoding to the format spacified.
	pf = fopen(dir,"r");	// Open file "read" mode
	if (pf == NULL){
		perror("fopen");
	}
	
	for (i = 0; i < Nsamples; i++){	// For every samples vector
		fscanf(pf,"%lf",&aux_d);	// Año 
		fscanf(pf,"%lf",&aux_d);	// Mes 
		fscanf(pf,"%lf",&aux_d);	// Día 
//		fscanf(pf,"%lf",&aux_d);	// Hora 
		fscanf(pf,"%lf",&((*X)[0][i]));		// Pillamos la hora
		fscanf(pf,"%lf",&aux_d);	// Validez
		fscanf(pf,"%lf",&((*t)[i]));	// Get the ALTURA into the output array
		for (j = 1; j < n; j++){	// For every parammeter of the vector.
			fscanf(pf,"%lf",&((*X)[j][i]));
		}
		if (aux_d <0) {		// aux_d = -1  -> Dato no valido
			i--;
		}
	}
	fclose(pf);
}

// In this version, the output parameter is the last one of the first file.
void Super_Load2(double ***X_tr, double **t_tr, double ***X_te, double **t_te, Input_params *Input_p) {
	int i, j;
	FILE *pf;
	
//******** Allocate memmory *********** 
//  X_tr[param][muestras]

	*X_tr = (double **) malloc(Input_p->n_input*sizeof(double*));
	if (*X_tr == NULL){
		perror("Malloc");
	}
	for (i = 0; i < n; i++){
		(*X_tr)[i] = (double *)malloc(Input_p->Ntrain*sizeof(double));
		if ((*X_tr)[i] == NULL){
			perror("Malloc");
		}
	}
	*t_tr = (double *) malloc(Input_p->Ntrain*sizeof(double));
	if (t_tr == NULL){
		perror("Malloc");
	}	
	
	*X_te = (double **) malloc(Input_p->n_input*sizeof(double*));
	if (*X_te == NULL){
		perror("Malloc");
	}
	for (i = 0; i < n; i++){
		(*X_te)[i] = (double *)malloc(Input_p->Ntest*sizeof(double));
		if ((*X_te)[i] == NULL){
			perror("Malloc");
		}
	}
	*t_te = (double *) malloc(Input_p->Ntest*sizeof(double));
	if (t_te == NULL){
		perror("Malloc");
	}	
	
//printf("Memoria para la matriz basica principal reservada\n");
//printf("Dimensiones %i x %i \n",n,Nsamples); 
	
// Read the vectors float by float accoding to the format spacified.

	if (Input_p->test_input_dir[0] != 0 ) {	// If we have specified separete files for train and test
		// Read training Input vectors 
		pf = fopen(Input_p->train_input_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen");
		}
		for (i = 0; i < Input_p->Ntrain; i++){	// For every samples vector
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_tr)[j][i]));
			}
			fscanf(pf,"%lf",&((*t_tr)[i]));
		}
		fclose(pf);
		
		// Read testing Input vectors 
		pf = fopen(Input_p->test_input_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen");
		}
		for (i = 0; i < Input_p->Ntest; i++){	// For every samples vector
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_te)[j][i]));
			}
			fscanf(pf,"%lf",&((*t_te)[i]));
		}
		fclose(pf);	
	}
	else {									// If we havent specified separate files
			// Read training Input vectors 
		pf = fopen(Input_p->train_input_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen ");
			printf("Error oppening %s \n",Input_p->train_input_dir); 
		}
		for (i = 0; i < Input_p->Ntrain; i++){	// For every samples vector
			fscanf(pf,"%lf",&((*t_tr)[i]));
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_tr)[j][i]));
			}
			
		}
		for (i = 0; i < Input_p->Ntest; i++){	// For every samples vector
			fscanf(pf,"%lf",&((*t_te)[i]));
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_te)[j][i]));
			}
		
		}	
		
	}
}


void Super_LoadX(double ***X_tr, double **t_tr, double ***X_te, double **t_te, Input_params *Input_p) {
	int i, j;
	FILE *pf;
	double aux_d;
//******** Allocate memmory *********** 
//  X_tr[param][muestras]

	*X_tr = (double **) malloc(Input_p->n_input*sizeof(double*));
	if (*X_tr == NULL){
		perror("Malloc");
	}
	for (i = 0; i < n; i++){
		(*X_tr)[i] = (double *)malloc(Input_p->Ntrain*sizeof(double));
		if ((*X_tr)[i] == NULL){
			perror("Malloc");
		}
	}
	*t_tr = (double *) malloc(Input_p->Ntrain*sizeof(double));
	if (t_tr == NULL){
		perror("Malloc");
	}	
	
	*X_te = (double **) malloc(Input_p->n_input*sizeof(double*));
	if (*X_te == NULL){
		perror("Malloc");
	}
	for (i = 0; i < n; i++){
		(*X_te)[i] = (double *)malloc(Input_p->Ntest*sizeof(double));
		if ((*X_te)[i] == NULL){
			perror("Malloc");
		}
	}
	*t_te = (double *) malloc(Input_p->Ntest*sizeof(double));
	if (t_te == NULL){
		perror("Malloc");
	}	
	
//printf("Memoria para la matriz basica principal reservada\n");
//printf("Dimensiones %i x %i \n",n,Nsamples); 
	
// Read the vectors float by float accoding to the format spacified.

	if (Input_p->test_input_dir[0] != 0 ) {	// If we have specified separete files for train and test
		// Read training Input vectors 
		pf = fopen(Input_p->train_input_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen");
		}
		for (i = 0; i < Input_p->Ntrain; i++){	// For every samples vector
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_tr)[j][i]));
			}
			fscanf(pf,"%lf",&((*t_tr)[i]));
		}
		fclose(pf);
		
		// Read testing Input vectors 
		pf = fopen(Input_p->test_input_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen");
		}
		for (i = 0; i < Input_p->Ntest; i++){	// For every samples vector
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_te)[j][i]));
			}
			fscanf(pf,"%lf",&((*t_te)[i]));
		}
		fclose(pf);	
	}
	else {									// If we havent specified separate files
			// Read training Input vectors 
		pf = fopen(Input_p->train_input_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen "); 
			printf("Error oppening %s \n",Input_p->train_input_dir); 
		}
		for (i = 0; i < Input_p->Ntrain; i++){	// For every samples vector
			fscanf(pf,"%lf",&((*t_tr)[i]));
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_tr)[j][i]));
			}
			for (j = 0; j < 60; j++) {	// To avoid the parameters of that boyamake
				fscanf(pf,"%lf",&aux_d);
			}
			
		}
		for (i = 0; i < Input_p->Ntest; i++){	// For every samples vector
			fscanf(pf,"%lf",&((*t_te)[i]));
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_te)[j][i]));
			}
			for (j = 0; j < 60; j++) {	// To avoid the parameters of that boyamake
				fscanf(pf,"%lf",&aux_d);
			}
		}	
		
	}
}

void Super_LoadX2(double ***X_tr, double **t_tr, double ***X_te, double **t_te, Input_params *Input_p) {
	int i, j;
	FILE *pf;
	double aux_d;
//******** Allocate memmory *********** 
//  X_tr[param][muestras]

	*X_tr = (double **) malloc(Input_p->n_input*sizeof(double*));
	if (*X_tr == NULL){
		perror("Malloc");
	}
	for (i = 0; i < n; i++){
		(*X_tr)[i] = (double *)malloc(Input_p->Ntrain*sizeof(double));
		if ((*X_tr)[i] == NULL){
			perror("Malloc");
		}
	}
	*t_tr = (double *) malloc(Input_p->Ntrain*sizeof(double));
	if (t_tr == NULL){
		perror("Malloc");
	}	
	
	*X_te = (double **) malloc(Input_p->n_input*sizeof(double*));
	if (*X_te == NULL){
		perror("Malloc");
	}
	for (i = 0; i < n; i++){
		(*X_te)[i] = (double *)malloc(Input_p->Ntest*sizeof(double));
		if ((*X_te)[i] == NULL){
			perror("Malloc");
		}
	}
	*t_te = (double *) malloc(Input_p->Ntest*sizeof(double));
	if (t_te == NULL){
		perror("Malloc");
	}	
	
//printf("Memoria para la matriz basica principal reservada\n");
//printf("Dimensiones %i x %i \n",n,Nsamples); 
	
// Read the vectors float by float accoding to the format spacified.

	if (Input_p->test_input_dir[0] != 0 ) {	// If we have specified separete files for train and test
		// Read training Input vectors 
		pf = fopen(Input_p->train_input_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen");
		}
		for (i = 0; i < Input_p->Ntrain; i++){	// For every samples vector
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_tr)[j][i]));
			}
			fscanf(pf,"%lf",&((*t_tr)[i]));
		}
		fclose(pf);
		
		// Read testing Input vectors 
		pf = fopen(Input_p->test_input_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen");
		}
		for (i = 0; i < Input_p->Ntest; i++){	// For every samples vector
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_te)[j][i]));
			}
			fscanf(pf,"%lf",&((*t_te)[i]));
		}
		fclose(pf);	
	}
	else {									// If we havent specified separate files
			// Read training Input vectors 
		pf = fopen(Input_p->train_input_dir,"r");	// Open file "read" mode
		if (pf == NULL){
			perror("fopen "); 
			printf("Error oppening %s \n",Input_p->train_input_dir); 
		}
		for (i = 0; i < Input_p->Ntrain; i++){	// For every samples vector
			fscanf(pf,"%lf",&((*t_tr)[i]));
			
			for (j = 0; j < 14; j++) {	// To avoid the parameters of that boyamake
				fscanf(pf,"%lf",&aux_d);
			}
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_tr)[j][i]));
			}
			
		}
		for (i = 0; i < Input_p->Ntest; i++){	// For every samples vector
			fscanf(pf,"%lf",&((*t_te)[i]));
			
			for (j = 0; j < 14; j++) {
				fscanf(pf,"%lf",&aux_d);
			}
			for (j = 0; j < Input_p->n_input; j++){			// For every parammeter of the vector.
				fscanf(pf,"%lf",&((*X_te)[j][i]));
			}
		
		}	
		
	}
}

void shuffle_train_test(double **X_tr, double *t_tr, double **X_te, double *t_te, Input_params *Input_p ) {
	double aux_vector[Input_p->n_input];
	double aux_output;
	int itr, ite;
	
	int i,j,k;
	
	for (i = 0; i < 2*Input_p->Ntest; i++) {
		itr = rand() % Input_p->Ntrain;
		ite = rand() % Input_p->Ntest;
		
		for (k = 0; k < Input_p->n_input; k++) {	// Swap Input and output vectors
			aux_vector[k] = X_tr[k][itr];
			X_tr[k][itr] = X_te[k][ite];
			X_te[k][ite] = aux_vector[k];		
		}	
		aux_output = t_tr[itr];
		t_tr[itr] = t_te[ite];
		t_te[ite] = aux_output;	
	}
	
}
