/*
 * elm.c
 *
 *  Created on: Nov 3, 2013
 *      Author: alvaro
 */

#include "../ELMheader.h"
#define SIGMOID	0
#define SINE	1
#define COSINE	2
#define TANH	3
#define LINEAL	4

#include <math.h>

void set_WeightBias(gsl_matrix **W, gsl_vector **b, int Nh, int n){
	unsigned int i,j;
	struct timeval t1;	// Random number to feed the pseudorandom function
	gsl_rng *r ;		// Pseudo random number generated.
			
	r = gsl_rng_alloc (gsl_rng_ranlux);
	gettimeofday(&t1,NULL);
	gsl_rng_set(r,(unsigned long)t1.tv_usec*t1.tv_sec);
	
	*W = gsl_matrix_alloc((size_t)Nh,n);
	*b = gsl_vector_alloc((size_t)Nh);
	
	// -------> Give random values to weights and bias <-------- 
	for(i=0;i< (*W)->size1; i++){
		
		gsl_vector_set(*b,i,2*gsl_rng_uniform(r)-1);// Offset always positive ??
		
		for(j=0;j< (*W)->size2; j++){
			gsl_matrix_set(*W,i,j,2*gsl_rng_uniform(r)-1);
		}
	}		
	gsl_rng_free(r); 
}

gsl_matrix* get_H_matrix(gsl_matrix *Xtrain, int Nh,gsl_matrix *W,gsl_vector *b, unsigned int act_func){
	unsigned int i;
	unsigned int j;
	double aux_double;
	
	gsl_matrix *H;		// H matrix with the response of the hidden neurons	(Returned, not freed)

	gsl_vector *w;		// Aux variable that will hold each nwuron weight
	gsl_vector *x;		// Aux variable that will hold each input vector

	//Reserve memmory for vectors and matrix
	H = gsl_matrix_alloc((size_t)Xtrain->size1,Nh);	// Ntrain x Nh

	x = gsl_vector_alloc((size_t)Xtrain->size2);	// n elements, 1 per intput
	w = gsl_vector_alloc((size_t)Xtrain->size2);	// n elements, 1 per intput


	// The vector bias b has one bias for every neuron and we only have to 
	// calculate it once and use every differente value.
	// The vector w has to change from neuron to neuron ???
	
	for(i=0; i < H->size1; i++){		// For evey input vector
		gsl_matrix_get_row(x,Xtrain,i);	// Get the input vector
		
		for(j = 0; j < Nh; j++){	// solve g(w*h + b) for every neuron
			gsl_matrix_get_row(w,W,j);	// Get the weights vector
			gsl_vector_mul(w,x);		// W(1xn)X(nx1) = 1x1 Multiply weight by input 
			aux_double = sum_gsl_vector(w) + gsl_vector_get(b,j); 
			switch (act_func){
				case SIGMOID:
					aux_double = 1/(1+exp(-aux_double)); //sigmoid function
					break;
				case SINE:
					aux_double = sin(aux_double); //sine function
					break;
				case COSINE:
					aux_double = cos(aux_double); //cosine function
					break;
				case TANH:
					aux_double = (exp(aux_double) - exp(-aux_double)) / (exp(aux_double) + exp(-aux_double)); //cosine function
					break;
				case LINEAL:
					aux_double = aux_double; //cosine function
					break;
			}
			gsl_matrix_set(H,i,j,aux_double);
		}
	}

	gsl_vector_free(x);
	gsl_vector_free(w);

// print_gsl_matrix(H);
	
	return H;
}

gsl_matrix* get_Hinv_matrix(gsl_matrix *H){	// SVG metho
	
	unsigned int i;
	gsl_matrix *Hinv;		// Penrouse-Moore inverse
	gsl_matrix *V;			// 
	gsl_matrix *Sinv;
	gsl_matrix *Aux;
	gsl_vector *s;
	gsl_vector *work;
	gsl_matrix *Hcopy;	// Copy of H coz the SVG function changes it
	
	//------> Allocate memmory <-----------
	Hinv 	= gsl_matrix_alloc(H->size2,H->size1);
	Sinv 	= gsl_matrix_alloc(H->size2,H->size2);
	V 		= gsl_matrix_alloc(H->size2,H->size2);
	s 		= gsl_vector_alloc(H->size2);
	work 	= gsl_vector_alloc(H->size2);
	Aux 	= gsl_matrix_alloc(H->size2,H->size1);
	
	Hcopy = gsl_matrix_alloc(H->size1,H->size2);
	
	//---------> Algoritm <--------------
	gsl_matrix_memcpy (Hcopy, H);
	gsl_linalg_SV_decomp(Hcopy,V,s,work); //H is replaced by U

	for(i = 0;i < Sinv->size1; i++){
		gsl_matrix_set(Sinv,i,i,1/gsl_vector_get(s,i));
	}
	
	gsl_blas_dgemm(CblasNoTrans,CblasConjTrans,1,Sinv,Hcopy,0,Aux);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,Aux,0,Hinv);

	//------> Free memmory <-----------
	gsl_matrix_free(Sinv);
	gsl_matrix_free(V);
	gsl_vector_free(s);
	gsl_vector_free(work);
	gsl_matrix_free(Aux);
	
	gsl_matrix_free(Hcopy);
	// gsl_matrix_free(H);	// We wont need ir anymore
	return Hinv;
}

gsl_vector* get_beta_vector(gsl_matrix *Hinv, gsl_vector *ytrain){
	gsl_vector *beta;		// Beta of the hidden neurons
	beta = gsl_vector_alloc(Hinv->size1);
	gsl_blas_dgemv(CblasNoTrans,1,Hinv,ytrain,0,beta);
	
	// gsl_matrix_free(Hinv);
	return beta;
}

gsl_vector* test_ELM(gsl_matrix *Htest, gsl_vector *beta){
	
	gsl_vector *y_pred;	// Predicted results from the function	(returned, not freed)
	y_pred = gsl_vector_alloc(Htest->size1);
	
	//Compute neural network to get the ypred
	gsl_blas_dgemv(CblasNoTrans,1,Htest,beta,0,y_pred);

	return y_pred;
}

double check_psudoinverse(gsl_matrix *H, gsl_matrix *Hinv){
	unsigned int i;		// Checks for property:  H*Hinv*H = H
	double error = 0;

	// H = (Ntrain x n) 		Hinv = (n x Ntrain) 
	gsl_matrix *Aux1;		// Aux1 = H x Hinv = (Ntrain x Ntrain)
	gsl_matrix *Aux2;		// Aux2 = Aux1 x H = (Ntrain x n) 
	gsl_vector *aux_vector;	// Vector for calculating the error
	
	Aux1 = gsl_matrix_alloc(H->size1,Hinv->size2);
	Aux2 = gsl_matrix_alloc(Aux1->size1, H->size2);
	aux_vector = gsl_vector_alloc(Aux2->size2);
	
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,H,Hinv,0,Aux1);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Aux1,H,0,Aux2);
	
	gsl_matrix_sub (Aux2, H);	// Compare result to ideal
	
	for (i = 0; i < Aux2->size1; i++){
		gsl_matrix_get_row(aux_vector,Aux2,i);	// Get the row
		error += sum_abs_gsl_vector(aux_vector) ;
	}
	error = error/(Aux2->size1*Aux2->size2);

	//print_gsl_matrix(Aux2);	
	// printf("\n El error es %lf \n",error);
	 
	gsl_matrix_free(Aux1);
	gsl_matrix_free(Aux2);
	gsl_vector_free(aux_vector); 
	
	return error;
}

gsl_matrix* get_Hinv_matrix2(gsl_matrix *H){
	
	int k, transpose = 0;	// Dont transpose initially
	int r = 0;
	int m, n;
	double tol, aux_value;
	int np;			// np = n -1 to easyly adapt C to MATLAB
	
	gsl_matrix *Hinv;	// Penrouse-Moore inverse
	
	gsl_matrix *A;		// Matriz A for the algotihm
	gsl_vector *A_subvector;		// Subvector of matrix L
	gsl_vector *dA;		// Diagonal of matrix A
	
	gsl_matrix *L;		// Matriz L for the algotihm
	gsl_matrix *L_submatrix;		// Submatrix of matrix L
	gsl_vector *L_subvector;		// Subvector of matrix L
	gsl_vector *L_subproduct;		// Subvector of matrix L
	gsl_matrix *L_r;		// Subvector of matrix L

	gsl_matrix *Aux;
	gsl_matrix *Aux2;

	gsl_matrix *M;		// Matrix M
	gsl_matrix *Minv;		// Inverse of Matrix M	
			
	gsl_matrix *G;		// Matrix G
	
	m = H->size1;
	n = H->size2;

	Hinv 	= gsl_matrix_alloc(n,m);		// Returnes, not freed
	
	// We will reserve maximum memmory for all matrix and vectors and we will change
	// their attributes   size for vector  and  size1,size2 for matrix to change their size.
	
//------> Allocate general memmory <-----------

		
	if (m < n){			// 
		transpose = 1;
		n = m;			// n = min(m,n)  !!!!!!
	}	


	A			= gsl_matrix_alloc(n,n);		
	A_subvector = gsl_vector_alloc(n);		
	dA			= gsl_vector_alloc(n);		
	
	L			= gsl_matrix_alloc(n,n);			
	L_submatrix			= gsl_matrix_alloc(n,n);		
	L_subvector			= gsl_vector_alloc(n);		
	L_subproduct		= gsl_vector_alloc(n);		
	L_r			= gsl_matrix_alloc(n,n);		

	Aux			= gsl_matrix_alloc(n,n);
	Aux2		= gsl_matrix_alloc(n,n);

	M			= gsl_matrix_alloc(n,n);
	Minv		= gsl_matrix_alloc(n,n);		
	G			= gsl_matrix_alloc(n,n);	


//printf("Memory allocated \n");
	// ------------> Check if we have to transpose it <---------------
	if (transpose){											// m < n 
		gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,H,H,0,A);	// A = H*H' -> (mxn)x(nxm) = (mxm)
	}
	else {													// m > n 
		gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,H,H,0,A);	//A = H'*H (nxm)x(mxn) = (nxn)
	} 
	
//printf("Matriz H \n");
//print_gsl_matrix(H);

	// In both cases m or n will be the shortest size and for now on it will be  "n = shortest size"
	
	//Full rank Cholesky factorization of A
//	WE WILL HAVE GENERAL PORPOUSE MATRIX AND VECTOR WITH MEMMORY RESERVED FOR THE MAXIMUM SIZE AND WE WILL USE IT FOR DIFFERENT
// SIZES CHANGING THE VARIABLES -> SIZE 
	gsl_matrix_diag(A,dA);
//printf("Diagonal ");
//print_gsl_vector(dA);
	tol = abs(gsl_vector_get(dA,0))/1000000000;	// Tolerance
	gsl_matrix_set_zero(L);	// Set matrix to 0   NECESARY ?????????????????????????????????
	r = -1;		// So that at the beginning of the loop it is 0

//printf("A = H'*H \n");
//print_gsl_matrix(A);	
	//	k {0,n-1}	-> Longitud  n
	//	r {0,n-1}	-> Longitud  n
	np = n -1;
	for (k = 0; k < n; k++){	// For every row
		r++;
		L_subvector->size = n;			// Set its size to its maximum
		gsl_vector_set_zero(L_subvector);	// Set vector to 0	L = { 0 0 0 0 0 .. } E n 
		A_subvector->size = n;			// Set its size to its maximum
		gsl_vector_set_zero(A_subvector);	// Set vector to 0	L = { 0 0 0 0 0 .. } E n 
		L_subproduct->size = n;			// Set its size to its maximum
		gsl_vector_set_zero(L_subproduct);	// Set vector to 0	L = { 0 0 0 0 0 .. } E n 
		
//	printf("*********************** Iteracion %i **********************\n",k);

//printf("Matriz L\n");															
//print_gsl_matrix(L);	
	
		if (r > 0) {		// Because when r <= 0 -->   L(k:n,1:(r-1)) =  {0}
			gsl_matrix_submatrix_cpy(L,L_submatrix, k, 0, (np-k)+1,r);	// L_submatrix = L(k:n,1:(r-1))  E (n-k)xr	
			
//printf("L_submatrix = L(k:n,1:(r-1))  E (n-k)xr -> %ix%i\n",(np-k)+1,r);
//print_gsl_matrix(L_submatrix);			
			gsl_matrix_subvector_row_get (L,L_subvector, k, 0,r);				// L_subvector = L(k,1:(r-1))'  E r
															// Using vectors in GSL there is no difference between L and L'
															// Using matrices it does matter
//printf("L_subvector = L(k,1:(r-1))' E r -> %ix%i\n", r,1);																	
//print_gsl_vector(L_subvector);	

			
			L_subproduct->size = (n-k);					// Set proper size for the result of the product					
			gsl_blas_dgemv(CblasNoTrans,1,L_submatrix,L_subvector,0,L_subproduct); // ((n-k)xr)x(rx1) = (n-k)x1;
//printf("L_subprudct  (n-k)x1 -> %ix%i\n",(np-k)+1,1);	
//print_gsl_vector(L_subproduct);	

		}
	
		gsl_matrix_subvector_col_get(A, A_subvector,k,k,(np-k)+1);	// A_subvector = A(k:n,k)		Column vector	(n-k+1)x1
//printf("A_subvector = A(k:n,k) \n");

//printf("A_subvector=  (n-k)x1 -> %ix%i\n",(n-k),1);	
//print_gsl_vector(A_subvector);			
		gsl_vector_dotsub(A_subvector,L_subproduct);	//	L_vector  = A_subvector + L_subproduct 
//printf("Vectores restados \n");	
//print_gsl_vector(A_subvector);					
																// Reusing the L_subvector to store final value
		A_subvector->size = n;									// Put same size as matriz (added values will be 0 as desired)													
	
//printf("Colocando vector de  longitud %i  en posicion (%i,%i) \n",n-k,r,k);		
		gsl_matrix_subvector_col_set(L, A_subvector, r, k, (np-k)+1);	// Set column vector L(k:n,r) = L_subvector
		
//printf("Cogiendo valor en posicion %i,%i \n",k,r, );	

		aux_value = gsl_matrix_get(L, k, r);		// L(k,r)
		
		if (aux_value > tol){		// If the value is greater than the tolerance  	aux_value > tol	
									// For not to be dividing by 0 ???
		
			aux_value = sqrt(aux_value);
			gsl_matrix_set(L, k,r,aux_value);	// L(k,r)=sqrt(L(k,r));

			if (k < (np)){
//				printf("*********** SHIT ******* \n");
				
				gsl_matrix_subvector_col_get(L,L_subvector, r, k+1, (np-(k+1))+1);	// Column vector r of L((k+1):n,r) 
				gsl_vector_scale_div (L_subvector, aux_value);		// L((k+1):n,r) / L(k,r)
				gsl_matrix_subvector_col_set(L, L_subvector, r,k+1, (np-(k+1))+1);		// L((k+1):n,r)=L((k+1):n,r)/L(k,r);
				
//printf("Vectore columna %i, en posicion %i, de longitud %i \n", r,k+1,(np-(k+1)) +1);	
//print_gsl_vector(L_subvector);		
			}
		}
		else{	
			// printf(" *********** FUCK YOU ******* \n" );
			// break;		// Added coz since its fucked we exit
			/* If this happens the algorith is pretty much fucked up */
			r--;
		}
	}
	r++;			// POR EL CAMBIO MATLAB C ????????????????
//	printf(" *************************************************** \n" );
//	printf(" ****************** SALIDO DEL BUCLE *************** \n" );
	
//printf("Conseguida matriz L (nxn) = (%i,%i)\n",L->size1,L->size2 );	
//print_gsl_matrix(L);
	
	gsl_matrix_submatrix_cpy(L,L_r, 0, 0,n,r);	//	L = L(:,1:r);	 (n x r)
//printf("Conseguida submatriz L_r (nxr) = (%i,%i) con origen (%i,%i)\n",L_r->size1,L_r->size2,0,0);	
//print_gsl_matrix(L_r);

	//  Computation of the generalized inverse of G
	M->size1 = L_r->size2;
	M->size2 = L_r->size2;	
// printf("n %i, r %i\n",n,r);	
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,L_r,L_r,0,M); // M = L'*L;	(rxr) // The function transposes the matriz L itself
//printf("Conseguida submatriz M =  L'*L	(rxr) = (%i,%i)\n",M->size1,M->size2 );	
//print_gsl_matrix(M);

	Minv->size1 = M->size2;
	Minv->size2 = M->size2;	

	gsl_inv_matrix(M,Minv);			//(rxr)
//printf("Inverse of M \n" );	
//print_gsl_matrix(Minv);

	// Calculate G = L*M*M*L';	
	Aux->size1 = Minv->size1;
	Aux->size2 = L_r->size1;
//printf("Aux =M*L' =  (%i x %i) (%i x %i) = (%i x %i)  \n",Minv->size1,Minv->size2,L_r->size2,L_r->size1,Aux->size1,Aux->size2 );	
		gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,Minv,L_r,0,Aux);	//	M*L';  (rxr)x(rxn) = (r x n)
		
	Aux2->size1 = L_r->size1;
	Aux2->size2 = Minv->size2;	
//printf("Aux2 =L*M=  (%i x %i) (%i x %i) = (%i x %i)  \n",L_r->size1,L_r->size2, Minv->size1,Minv->size2,Aux2->size1,Aux2->size2 );		
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,L_r,Minv,0,Aux2);	//  L*M;  (nxr)x(rxr) = (nxr)
//printf("Aux2 \n" );		
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Aux2,Aux,0,G);	// G = Aux2 * aux = (nxr) x (rxn) = (nxn)
//printf("G \n" );			
	if (transpose){

		gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,H,G,0,Hinv);		//  Y = H'*G   (nxm)x(mxm) = (nxm)
	}
	else{
			
		gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,G,H,0,Hinv);	// Y = G*H'   (nxn) x (n x m) = (nxm)
	}

	//------> Free memmory <-----------
	gsl_matrix_free(A);	
	gsl_vector_free(A_subvector);
	gsl_vector_free(dA);
	
	gsl_matrix_free(L);		
	gsl_matrix_free(L_submatrix	);	
	gsl_vector_free(L_subvector	);
	gsl_vector_free(L_subproduct);
	gsl_matrix_free(L_r);		

	gsl_matrix_free(Aux);
	gsl_matrix_free(Aux2);

	gsl_matrix_free(M);
	gsl_matrix_free(Minv);		
	gsl_matrix_free(G);	
	
	return Hinv;
}

// La view es como un puntero a memoria que ya existe de una matriz. Otra forma de llamar una matriz o parte de esta.
// Es como crear 2 putneros que apuntan al mismo sitio pero solo se puede liberar a traves del que no es view


// Al hacer la transformacion MATLAB C, dado que el origen de las indices en matlab es 1 y en C es 0
//	L(1:n)			= 	C(0,n-1)	-> Longitud n
//	long = (n-1)+1	= 	long = (n-1-0)+1	


/*
double get_ELM_RMSE(int *chosen_vector, int n_chosen, int Nh, int FLAGS,
				 gsl_vector ** x_train_param,gsl_vector ** x_test_param,
				 gsl_vector *y_train, gsl_vector *y_test,
				 int t_max, int t_min){
*/
	
double get_ELM_RMSE(ELM_params *p){
	int i;				 
	struct timeval time_start,time_end;		// Time variables for general use		
	float time_passed;	
	double error_Hinv;
	int errores;
	
	//-------------> PARAM VARIABLES <-----------------
	unsigned int *chosen_vector = p->chosen_vector ;
	int n_chosen = p->n_chosen;
	int Nh = p->Nh ;
	int FLAGS = p->FLAGS; 
	gsl_vector ** x_train_param = p->x_train_param  ;
	gsl_vector ** x_test_param = p->x_test_param;
	gsl_vector *y_train = p->y_train;
	gsl_vector *y_test = p->y_test;
	int t_max = p->t_max;
	int t_min = p->t_min;
	
	//-------------> DYNAMIC RESERVED VARIABLES <-----------------		 
	gsl_matrix *Xtrain;		// Input training vector matrix
	gsl_matrix *Xtest;		// Input testing vector matrix
	
	gsl_matrix *H;			// H matrix (used for training, then freed and the used for testing)
	gsl_matrix *Hinv;		// Moore-penrouse inverse of H train matrix
	
	gsl_vector *y_pred;		// Predicted result vector
	gsl_vector *beta;		// Betas
	
	gsl_matrix *W;			// Weights matrix;
	gsl_vector *b;			// Bias vector (one for neuron hidden)
	
	gsl_vector *y_test_aux;		// Copy of y_test given as input coz this funcion changes it when 
								// we desnormalize it
//	print_gsl_matrix(Xtrain);
//	print_gsl_matrix(Xtest);
// ****************************************************************************************
// ********************************  ELM ALGORITHM ****************************************
// ****************************************************************************************

	if (n_chosen == 0){		// We have to at least select one parameter.
		n_chosen = 1;
		chosen_vector[0] = 1;	// Modify the vector
	}
	y_test_aux = gsl_vector_alloc(y_test->size);

	Xtrain = get_selected_input_matrix(x_train_param, chosen_vector, n_chosen);	// Generate Xtrain
	Xtest = get_selected_input_matrix(x_test_param, chosen_vector, n_chosen);	// Get Xtest matrix
	
	if (p->n_ELMs <= 0 ){
		p->n_ELMs = 1;
	}
	p->ERMS = 0;
	
//	printf("Numero de %i \n",p->n_ELMs);
	for (i = 0; i < p->n_ELMs; i++) {
	
	gsl_vector_memcpy(y_test_aux,y_test);
//---------------------> Generate the betas <-------------------------
	if ((FLAGS & SHOW_TIME_F)&&(i == 0)){
		printf("Setting random variables and Input matrix.. ");
		errores = gettimeofday( &time_start, NULL);
	}	
	set_WeightBias(&W, &b, Nh, n_chosen);		// Set Weights and bias (Reserved by function)
	
	if ((FLAGS & SHOW_TIME_F)&&(i == 0)){
		errores = gettimeofday( &time_end, NULL);
		time_passed = get_time_passed(time_start,time_end);
		printf(" %f s\n",time_passed);
	}
	
//---------------------> Calculate train matrix H  <-------------------------
	if ((FLAGS & SHOW_TIME_F)&&(i == 0)){
		printf("Calculadno H... ");
		errores = gettimeofday( &time_start, NULL);
	}
	H = get_H_matrix(Xtrain,Nh,W,b, p->activation_f);				// Calculate H train
	
	if ((FLAGS & SHOW_TIME_F)&&(i == 0)){	
		errores = gettimeofday( &time_end, NULL);
		time_passed = get_time_passed(time_start,time_end);
		printf(" %f s\n",time_passed);	
	}
	
//---------------------> Calculate train matrix Hinv  <-------------------------
	if ((FLAGS & SHOW_TIME_F)&&(i == 0)){
		printf("Calculadno Hinv... ");
		errores = gettimeofday( &time_start, NULL);
	}
	Hinv = get_Hinv_matrix2(H);						// Calculate Hinv

	if ((FLAGS & SHOW_TIME_F)&&(i == 0)){
		errores = gettimeofday( &time_end, NULL);
		time_passed = get_time_passed(time_start,time_end);
		printf(" %f s ******************\n",time_passed);
	}

//---------------------> Checking Hinv <-------------------------
	if ((FLAGS & CHECK_HINV_F)&&(i == 0)) {
		if ((FLAGS & SHOW_TIME_F)&&(i == 0)){
			printf("Checking Hinv... ");
			errores = gettimeofday( &time_start, NULL);	
		}
		error_Hinv = check_psudoinverse(H, Hinv);				// Check pseudo inverse
		if ((FLAGS & SHOW_TIME_F)&&(i == 0)){
			errores = gettimeofday( &time_end, NULL);
			time_passed = get_time_passed(time_start,time_end);
			printf(" %f s\n",time_passed);	
		}
		printf("Error de la Hinv: %lf \n",error_Hinv);
	}
	
//---------------------> Calculating Beta solution <-------------------------
	if ((FLAGS & SHOW_TIME_F)&&(i == 0)){
		printf("Calculadno Beta... ");
		errores = gettimeofday( &time_start, NULL);
	}	
	beta = get_beta_vector(Hinv, y_train);		// Calculate Beta

	if ((FLAGS & SHOW_TIME_F)&&(i == 0)){
		errores = gettimeofday( &time_end, NULL);
		time_passed = get_time_passed(time_start,time_end);
		printf(" %f s\n",time_passed);
	}
	
//---------------------> Testing result <-------------------------
	if ((FLAGS & SHOW_TIME_F)&&(i == 0)){
		printf("Testeando... ");
		errores = gettimeofday( &time_start, NULL);
	}	
	gsl_matrix_free(H);										// Free Htrain
	H = get_H_matrix(Xtest,Nh,W,b, p->activation_f);			// Get the Htest matrix
	
	y_pred = test_ELM(H, beta);			// Get the predicted values

	desnormalize_array(y_pred->data, x_test_param[0]->size, t_min, t_max);	// Desnormalice
	desnormalize_array(y_test_aux->data, x_test_param[0]->size, t_min, t_max);
	
	p->ERMS += rmse_gsl_vector(y_pred,y_test_aux); // Calcula error
	
	if ((FLAGS & SHOW_TIME_F)&&(i == 0)){
		errores = gettimeofday( &time_end, NULL);
		time_passed = get_time_passed(time_start,time_end);
		printf(" %f s\n",time_passed);	
	}
	
//---------------------> Plotting data <-------------------------
	if ((FLAGS & PLOT_DATA_F) &&(i == 0)){
		if ((FLAGS & SHOW_TIME_F)&&(i == 0)){
			printf("Ploteando datos... ");
			errores = gettimeofday( &time_start, NULL);	
		}
		plot(y_test_aux,y_pred);							// Plot data
		print_ELM_data (p);
		if ((FLAGS & SHOW_TIME_F)&&(i == 0)){	
			errores = gettimeofday( &time_end, NULL);
			time_passed = get_time_passed(time_start,time_end);
			printf(" %f s\n",time_passed);			
		}			 
	}
			 		 
	//---------------------> Freeing memmory <-------------------------
	gsl_matrix_free(H);					// Free Hte
	gsl_vector_free(y_pred);		
	gsl_vector_free(beta);		
	gsl_matrix_free(W);			
	gsl_vector_free(b);			
	gsl_matrix_free(Hinv);						// Free Hinv
			
	}
	
	p->ERMS = p->ERMS / p->n_ELMs;
	
	gsl_vector_free(y_test_aux);
	gsl_matrix_free(Xtest);					// Free Xtest
	gsl_matrix_free(Xtrain);				// Free Xtrain
	return p->ERMS;				 
}

int ELMs_op (ELMs_data *ELMs_p){
	double ERMS_value;
	int i,k;
	unsigned int ** chosen_arrays;
	double* results;
	print_ELMs_values(ELMs_p);
	// ***************** WHITHOUT THREADS  ********************
	if (Thread_p.n_threads == 1) {
		ELM_p.FLAGS = ELMs_p->FLAGS;
		for ( k = 0; k < ELMs_p->n_Selected; k++){
			ELM_p.chosen_vector = ELMs_p->Selected[k];
			for ( i = 0; i < ELMs_p->n_ELMs; i++){
				ERMS_value = get_ELM_RMSE(&ELM_p);			// ELM !!!!		
				printf(" Resultado %lf \n", ERMS_value);
			}			
		}
	}
	// ***************** ELM THREADS ALGORITHM ********************
	else {
		chosen_arrays = (unsigned int **) malloc (sizeof(unsigned int *)*ELMs_p->n_ELMs);
		results = (double*) malloc (sizeof(double)*ELMs_p->n_ELMs);	
		for (i = 0; i < ELMs_p->n_ELMs; i++){
			chosen_arrays[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var_chosen);
		}
		
		for ( k = 0; k < ELMs_p->n_Selected; k++){
			for (i = 0; i < ELMs_p->n_ELMs; i++){
				bitvector_cpy(chosen_arrays[i], ELMs_p->Selected[k], n);
			}	
			threads_WorkOut(chosen_arrays,results, ELMs_p->n_ELMs, &Thread_p);		// Create threads and make them work !!!
			// Once the Thread input and output structures are set we call the threads.

			for (i = 0; i < ELMs_p->n_ELMs; i++){
				printf(" Resultado %lf \n", results[i]);
			}
		
			if (ELMs_p->get_distrib == 1) {
				ELM_p.FLAGS = 0;
				init_threads(&ELM_p, &Thread_p);
				get_ELM_random_distrib(ELMs_p->Selected[k],ELMs_p->n_ELMs_dis, ELMs_p->n_divs );
			}
			
		}
		for (i = 0; i < ELMs_p->n_ELMs; i++){
			free(chosen_arrays[i]);
		}
		free(chosen_arrays);
		free(results);
	}
	for (i = 0; i < ELMs_p->n_Selected; i++){
		free(ELMs_p->Selected[i]);
	}
	free(ELMs_p->Selected);	
	return 1;
}

int print_ELMs_values(ELMs_data *ELMs_p){
	int i;
	printf("********** ELMs values ***************** \n");
	printf("N ELMs: %i \n",ELMs_p->n_ELMs);		
	printf("N Selected: %i \n",ELMs_p->n_Selected);	
	for (i = 0; i < ELMs_p->n_Selected; i++) {
		print_bitvector(ELMs_p->Selected[i], n);
	}
	printf("Flags  %i \n",ELMs_p->FLAGS);		// Number of bits of the population	
	printf("Get distrib   %i \n",ELMs_p->get_distrib);		// Number of bits of the population	
	printf("********** ****************************************** \n");
	return 1;
}

int get_ELM_random_distrib(unsigned int * Selected, int n_ELMs, int div){
	int i,j;
	double max, min;
	
	double average = 0;
	double variance = 0;
	
	double * ERMS_values = (double *) malloc (sizeof(double)*n_ELMs);
	double * graph_values = (double *) malloc (sizeof(double)*div);
	
	unsigned int * order = (unsigned int *) malloc (sizeof(unsigned int)*n_ELMs);
	unsigned int ** Selected_copies = (unsigned int **) malloc (sizeof(unsigned int *)*n_ELMs);
	for (i = 0; i < n_ELMs; i++){
		Selected_copies[i] = (unsigned int *) malloc (sizeof(unsigned int)*n_var_chosen);
	}	

	//-------> Get the Fitness of every Feature <-------
	for (i = 0; i < n_ELMs; i++){
		bitvector_cpy(Selected_copies[i], Selected, n);
	}	
	
	threads_WorkOut(Selected_copies, ERMS_values, n_ELMs, &Thread_p );
	
	order_double(ERMS_values, order, n_ELMs);
	

	min = ERMS_values[n_ELMs -1];
	max = ERMS_values[0];
	
	for (i = 0; i < n_ELMs; i++){
		average += ERMS_values[i];
	}
	average /= n_ELMs;
	
	for (i = 0; i < n_ELMs; i++){
		variance += abs_d(ERMS_values[i] - average);
	}
	variance /= n_ELMs;
	variance /= average;
	printf("Variance %f \n",variance);
	//------------> Graph Stuff <-------------
	
	for (i = 0; i < div; i++){
		graph_values[i] = 0;
	}
	
	int pos;
	for (i = 0; i < n_ELMs; i++){
		pos = (ERMS_values[i] - min)/(max-min) * div;	// Normalize 0-1 and multiply by div
		if (pos == div) pos--;
		graph_values[pos]++;
	}
	
	multiple_graph mgra;
	gsl_vector **graph_ES;
	int n_curves = 1;
	
	graph_ES = (gsl_vector **)malloc((n_curves + 1)*sizeof(gsl_vector *));
	for (i = 0; i < n_curves + 1; i++){
		graph_ES[i] = gsl_vector_alloc(div);
	}
	for (i = 0; i < div; i++){
		graph_ES[0]->data[i] = min + i*(max - min)/div;
	}	
	memcpy(graph_ES[1]->data, graph_values, div*sizeof(double));
	
	mgra.y = graph_ES;
	mgra.n_curves = n_curves;
	mgra.graph_name = (char *)malloc(30*sizeof(char));
	mgra.curve_names= (char **)malloc(n_curves*sizeof(char*));
	
	mgra.X_axis_n = (char *)malloc(30*sizeof(char));
	mgra.Y_axis_n = (char *)malloc(30*sizeof(char));
	
	for (i = 0; i < n_curves; i++){
		mgra.curve_names[i] = (char *)malloc(30*sizeof(char));
	}
	sprintf(mgra.curve_names[0], "Prob. Distribution");
	
	sprintf(mgra.graph_name, "Prob. Distribution");
	
	sprintf(mgra.X_axis_n , "Ave: %f Var %f", average, variance);
	sprintf(mgra.Y_axis_n , "RMSE");

	plot_multiple_graph(&mgra);
	
	//------> Free memmory -----
	free(ERMS_values);
	free(graph_values);
	free(order);

	for (i = 0; i < n_ELMs; i++){
		free(Selected_copies[i]);
	}		
	free(Selected_copies);
	
	//--------> Free graph memmory <-----
	for (i = 0; i < n_curves + 1; i++){
			gsl_vector_free(graph_ES[i]);
	}
	free(graph_ES);
	 
	free(mgra.X_axis_n);
	free(mgra.Y_axis_n);
	free(mgra.graph_name);
	for (i = 0; i < n_curves; i++){
		free(mgra.curve_names[i]);
	}
	free(mgra.curve_names);
}

int print_ELM_values(Input_params *ELM_p){
	
	printf("********** ELM values ***************** \n");
	printf("Train_data_file %s \n",ELM_p->train_input_dir);		
	printf("Train_result_file %s \n",ELM_p->train_output_dir);		
	printf("Test_data_file	 %s \n",ELM_p->test_input_dir);		
	printf("Test_result_file %s \n",ELM_p->test_output_dir);		

	printf("Ntrain  %i \n",ELM_p->Ntrain);		// Number of bits of the population	
	printf("Ntest  %i \n",ELM_p->Ntest);		// Number of bits of the population	
	
	printf("n_input  %i \n",ELM_p->n_input);		// Number of bits of the population	
	printf("n_output %i \n",ELM_p->n_output);		// Number of bits of the population	
	
	printf("n_ELMs  %i \n",ELM_p->n_ELMs);		// Number of bits of the population	
	printf("Nh %i \n",ELM_p->Nh);		// Number of bits of the population	
	
	printf("activation_f  %i \n",ELM_p->activation_f);		// Number of bits of the population	
	printf("n_threads %i \n",ELM_p->n_threads);		// Number of bits of the population	
	
	printf("********** ****************************************** \n");	
}

