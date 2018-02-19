
#include "../ELMheader.h"

void print_gsl_vector(gsl_vector *v){
	unsigned int k=0;

	printf("(");
	for(k=0;k<v->size;k++){
		printf("%g",gsl_vector_get(v,k));
		if(k!=v->size-1) printf(" ");
	}

	printf(")\n");
}
void print_gsl_matrix(gsl_matrix *M){
	int k1=0,k2=0,n1=0,n2=0;

	n1=M->size1;
	n2=M->size2;

	printf("\n");
	printf("[");
	for(k1 = 0;k1 < n1; k1++){
	
		for(k2 = 0;k2 < n2; k2++){
			printf("%lf",gsl_matrix_get(M,k1,k2));
			if(k2!=n2-1) printf(" ");
		}
		if (k1 != n1 -1){
			printf(";\n");
		}
	}
	printf("]\n");
	printf("\n");
}

double sum_gsl_vector(gsl_vector *x){
	unsigned int k=0;
	double sum=0;

	for (k = 0;k <x->size; k++)
		sum += gsl_vector_get(x,k);
	return sum;
}
double sum_abs_gsl_vector(gsl_vector *x){
	unsigned int k = 0;
	double sum = 0;

	for(k = 0;k <x->size; k++)
		sum += abs(gsl_vector_get(x,k));
	return sum;
}

double rmse_gsl_vector(gsl_vector *x,gsl_vector *y){
	unsigned int k=0;
	double mse=0;

	for(k = 0; k < x->size; k++){
		mse += pow(gsl_vector_get(x,k)-gsl_vector_get(y,k),2);
	}
	mse /= x->size;

	return sqrt(mse);
}

void gsl_matrix_subvector_col_get(gsl_matrix *M, gsl_vector* v, int num_col, int offset, int size){
	// This funcion doesnt reserve memmory, it just copies the values
	int i;
	double aux_value;
	v->size = size;
//printf("Obtenido columna %i, (%i - %i) \n \n",num_col,offset, offset+size-1);
	for (i = 0; i < size; i++){	// Se podria hacer con memcpy
		aux_value = gsl_matrix_get(M, i + offset, num_col);
		gsl_vector_set(v,i,aux_value);
	}
	
}

int gsl_matrix_subvector_row_get(gsl_matrix *M, gsl_vector* v, int num_row, int offset, int size){
	// This funcion doesnt reserve memmory, it just copies the values
	// For this function to work, the vector size must be equal or greater to the one we will
	// extract from the matrix. 
	// This funcion will set the size of the vector to the one written. 
	// THIS DOESNT NOT AFFECT THE FREEING OF MEMMORY COZ THE ONLY THING THAT MATTERS IS THE ->DATA 
	// VECTOR ITSELF 
	int i;
	double aux_value;
	v->size = size;
	
//printf("Obtenido fila %i, (%i - %i) \n \n",num_row,offset, offset+size-1);
	
	for (i = 0; i < size; i++){	// Se podria hacer con memcpy
		aux_value = gsl_matrix_get(M, num_row, i + offset);
		gsl_vector_set(v,i,aux_value);	
	}
	return size;
}	
// Gives values to a subvector row of a matrix
int gsl_matrix_subvector_row_set(gsl_matrix *M, gsl_vector* v, int num_row, int offset, int size){
	int i;
	double aux_value;
	for (i = 0; i < size; i++){	// Se podria hacer con memcpy
		aux_value = gsl_vector_get(v,i);	
		gsl_matrix_set(M, num_row, i + offset,aux_value);
	}
	v->size = size;
	return size;	
}
// Gives values to a subvector  column of a matrix
int gsl_matrix_subvector_col_set(gsl_matrix *M, gsl_vector* v, int num_col, int offset, int size){
	int i;
	double aux_value;
	for (i = 0; i < size; i++){	// Se podria hacer con memcpy
		aux_value = gsl_vector_get(v,i);	
		gsl_matrix_set(M, i + offset,num_col,aux_value);
	}
	v->size = size;
	return size;	
}

int gsl_matrix_submatrix_cpy(gsl_matrix *M, gsl_matrix *m, int offrow, int offcol, int sizerow, int sizecol){
	// This funcion doesnt reserve memmory, it just copies the values
	int i,j;
	double aux_value;
//	printf("offrow %i, offcol %i, sizerow %i, sizecol %i \n \n",offrow,offcol,sizerow,sizecol);
	
	m->size1 = sizerow;
	m->size2 = sizecol;
	for (i = 0; i < sizerow; i++){	// Se podria hacer con memcpy
		for (j = 0; j < sizecol; j++){	// Se podria hacer con memcpy	
			aux_value = gsl_matrix_get(M, i + offrow, j + offcol);
			gsl_matrix_set(m,i,j, aux_value);	
		}
	}

	return ((m->size1) * (m->size1));
}

int gsl_inv_matrix(gsl_matrix *M,gsl_matrix *Minv ){
	// This function calculates the inverse matrix of M and places it into Minv
	// Matrices M and Minv must be same size and memmory reserved.
	unsigned int n = M->size1;
	gsl_matrix *Minvcheck;	// Matrix to check inverse of M	
	
	
	Minv->size1 = n;
	Minv->size2 = n;	
	
	Minvcheck = gsl_matrix_alloc (n,n);
	gsl_matrix_memcpy (Minv, M);
	
	gsl_linalg_cholesky_decomp (Minv);
	gsl_linalg_cholesky_invert (Minv);
	
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Minv,M,0,Minvcheck);
//print_gsl_matrix(Minvcheck);

	gsl_matrix_free(Minvcheck);
	
	return 1;
}

int gsl_vector_dotsum(gsl_vector * a, gsl_vector * b){
	unsigned int i;
//	printf("XXXX %i, %i \n", a->size, b->size);
	double aux_value, aux_value2;
	for (i = 0; i < a->size; i++){	// Se podria hacer con memcpy
		aux_value = gsl_vector_get(b,i);	
		aux_value2 = gsl_vector_get(a,i);	
		gsl_vector_set(a,i,aux_value+aux_value2);
	}
	return  a->size;		
}

int gsl_vector_dotsub(gsl_vector * a, gsl_vector * b){
	unsigned int i;
//	printf("XXXX %i, %i \n", a->size, b->size);
	double aux_value, aux_value2;
	for (i = 0; i < a->size; i++){	// Se podria hacer con memcpy
		aux_value = gsl_vector_get(a,i);	
		aux_value2 = gsl_vector_get(b,i);	
		gsl_vector_set(a,i,aux_value-aux_value2);
	}
	return  a->size;		
}
int gsl_matrix_diag(gsl_matrix * m, gsl_vector * d){
	unsigned int i;
	double aux;
	for (i = 0; i < m->size1; i++) {
		aux = gsl_matrix_get(m, i, i);
		gsl_vector_set(d,i,aux);
	}
	return d->size;
}
int gsl_vector_scale_div (gsl_vector *v, double div){
	unsigned int i;
	double aux;
	for (i = 0; i < v->size; i++) {
		aux = gsl_vector_get(v, i);
		gsl_vector_set(v,i,aux/div);
	}
	return v->size;	
}	

