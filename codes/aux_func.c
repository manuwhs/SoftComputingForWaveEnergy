// Librerias estándar de C
#include <stdio.h>
#include <stdlib.h>

// Librerias de llamadas al sistema de Unix
#include <unistd.h>
#include <fcntl.h>		
#include <sys/types.h>	
#include <sys/stat.h>	
#include "../ELMheader.h"

/* Le das el puntero a char * IN y te devuelve el caracter hexadecimal correspondiente a 
ese char y al siguiente */
unsigned int min_uint(unsigned int a, unsigned int b){
	if (a > b) 
		return b;
	return a;
}

double min_d(double a, double b){
	if (a > b) 
		return b;
	return a;
}
void print_array(double *v, int len){
	int i;
	printf("[ ");
	for (i = 0; i < len; i++){
		printf("%lf, ",v[i]);
	}
	printf("] \n");
}

void print_matrix(double **v, int len1, int len2){
	int i,j;
	for (i = 0; i < len1; i++){	// For every training vector
		printf("[ ");
		for (j = 0; j < len2; j++){	// For every parammeter of the vector.
			printf("%f, ", v[i][j]);
		}
		printf("] \n");
	}
}


int roof_int (int num, int den) {
	float f = num/den;
	int f_int = (int)f;
	if ( num % den)
		f_int++;
	return f_int;
}

double abs_d (double d){
	if (d > 0.0) 
		return d;
	else 
		return -d;
}


/* Le pasas un array y su longitud y te lo inicializa a 0 */


void zero_char (char * array, int length) {
	int i ;
	for (i = 0; i < length; i++) {
		*(array + i) = 0;
	}
}

/* ---------------------------------------------------------------------------------------------------*/

void zero_int (int * array, int length) {
	int i ;
	for (i = 0; i < length; i++) {
		*(array + i) = 0;
	}
}

/* ---------------------------------------------------------------------------------------------------*/

int str2dec(char  * cad, char digit) {
	int out = 0;
	int i, j;
	int aux;
	
	for (i = 0; i < digit; i++) { 
		aux = cad[i] - '0';
		if((digit -i ) > 0 ) {
			for (j = 1; j < digit - i; j++) {
				aux = aux * 10;
			}
		}
		out += aux;
	}
	return out;
}

/* ---------------------------------------------------------------------------------------------------*/
/* Ordena de mayor a menor un array de int */
int ordenar_int (int * lista, int * orden, int len) {
 
	int aux_big; 
	int aux_pos, aux; 
	int i,j;
	
	for (i = 0; i < len ; i++) {
		orden[i] = i;
	}

	for(i = 0 ;i < len - 1; i++) {
		aux_big = lista[i];
		aux_pos = i;

		for (j = i + 1; j < len; j++) {

			if(aux_big < lista[j]) {
				aux_big = lista [j];
				aux_pos = j;
			}
		 }
		lista[aux_pos] = lista[i];
		lista[i] = aux_big;

		aux = orden[aux_pos];
		orden[aux_pos] = orden[i];
		orden[i] = aux;
	 }
	return 1;
}
/* Ordena de mayor a menor un array de int */
int ordenar_float (float * lista, int * orden, int len) {
 
	float aux_big; 
	int aux_pos, aux; 
	int i,j;
	
	for (i = 0; i < len ; i++) {
		orden[i] = i;
	}

	for(i = 0 ;i < len - 1; i++) {
		aux_big = lista[i];
		aux_pos = i;

		for (j = i + 1; j < len; j++) {

			if(aux_big < lista[j]) {
				aux_big = lista [j];
				aux_pos = j;
			}
		 }
		lista[aux_pos] = lista[i];
		lista[i] = aux_big;

		aux = orden[aux_pos];
		orden[aux_pos] = orden[i];
		orden[i] = aux;
	 }
	return 1;
}
/* ---------------------------------------------------------------------------------------------------*/
/* Ordena de mayor a menor un array de int */
int order_double (double * lista, int * orden, int len) {
 
	double aux_big; 
	int aux_pos, aux; 
	int i,j;
	
	for (i = 0; i < len ; i++) {
		orden[i] = i;
	}

	for(i = 0 ;i < len - 1; i++) {
		aux_big = lista[i];
		aux_pos = i;

		for (j = i + 1; j < len; j++) {

			if(aux_big < lista[j]) {
				aux_big = lista [j];
				aux_pos = j;
			}
		 }
		lista[aux_pos] = lista[i];
		lista[i] = aux_big;

		aux = orden[aux_pos];
		orden[aux_pos] = orden[i];
		orden[i] = aux;
	 }
	return 1;
}
/* ---------------------------------------------------------------------------------------------------*/
/* Reordena un array tal y como le indica otro */
/* orden[i] tiene la nueva posicion en la que esta el antiguo elemento i 
 * Por ejemplo si el elemento [0] de la antigua lista esta en la posicion [43] de la nueva
 * orden[0] = 43*/
int reordenar_int (int * lista, int * orden, int len) {
 	int nueva_lista [len];
	int i;
	for(i = 0 ;i < len ; i++) {
		nueva_lista[i] = lista[orden[i]];
	 }
	 copy_vector_int(lista, nueva_lista, len);
	return 1;
}
/* ---------------------------------------------------------------------------------------------------*/
/* Reordena un array tal y como le indica otro */
int reordenar_float (float * lista, int * orden, int len) {
 	float nueva_lista [len];
	int i;
	for(i = 0 ;i < len ; i++) {
		nueva_lista[i] = lista[orden[i]];
	 }
	copy_vector_float(lista, nueva_lista, len);
	return 1;
}
/* ---------------------------------------------------------------------------------------------------*/
int desreordenar_char (char * lista, int * orden, int len) {
 	char nueva_lista [len];
	int i;
	for(i = 0 ;i < len ; i++) {
		nueva_lista[orden[i]] = lista[i];
	 }
	 copy_vector_char(lista, nueva_lista, len);
	return 1;
}

/* ---------------------------------------------------------------------------------------------------*/




int dec2str(char  * cad, int digit) {
	int i, j;
	int aux;
	
	i = 0;
	aux = digit;
	while (aux != 0 ) {
		cad[i] = aux % 10 + '0';
		aux = aux / 10;
		i++;
	}
	cad[i] = '\0';
	
	printf("%s \n", cad);
	for (j = 0; j < i/2; j++){
		aux = cad[j];
		cad[j] = cad[i - 1 -j];
		cad[i - 1 -j] = aux;
	}
	printf("%s \n", cad);
	return 1;
}

/* ---------------------------------------------------------------------------------------------------*/


int max(int a, int b) { 
	return (a > b)? a : b;
}
/* ---------------------------------------------------------------------------------------------------*/
float maxf (float * array, int num) {
	int i;
	float aux = array[0];

	for (i = 1; i < num; i++) {
		if (array[i] > aux) {
			aux = array[i];
		}
	}
	return aux;
}
/* ---------------------------------------------------------------------------------------------------*/
int findint(int * array, int len, int value) {
	int i;
	for (i = 0; i < len ; i++ ) {
		if (array[i] == value ) {
//			printf("Found on position: %i\n", i);
			return i;
		}
	}
	return -1;  /* Not found */
}

/* ---------------------------------------------------------------------------------------------------*/

int copy_vector_char(char * duplicado, char * original, int num) {
	int i;
	for (i = 0; i < num; i++) {
		duplicado[i] = original[i];
	}
	return 1;
}

/* ---------------------------------------------------------------------------------------------------*/

int copy_vector_int(int * duplicado, int * original, int num) {
	int i;
	for (i = 0; i < num; i++) {
		duplicado[i] = original[i];
	}
	return 1;
}
/* ---------------------------------------------------------------------------------------------------*/

int copy_vector_float(float * duplicado, float * original, int num) {
	int i;
	for (i = 0; i < num; i++) {
		duplicado[i] = original[i];
	}
	return 1;
}

/* ---------------------------------------------------------------------------------------------------*/


int count_char(char * array, int len, int value) {
	int i;
	int count = 0;
	for (i = 0; i < len ; i++ ) {
		if (array[i] == value ) {
//			printf("Found on position: %i\n", i);
			count++;
		}
	}
	return count;  /* Not found */
}

/* ---------------------------------------------------------------------------------------------------*/

void get_matrix_row_minmax(double **v,double *min, double *max, int len1, int len2){
	double aux_min, aux_max;
	int i,j;
	for(i = 0;i < len1; i++){
		aux_min = *(*(v+i)+0);	// v[i][0]
		aux_max = *(*(v+i)+0);	// v[i][0]
		
		for(j = 0;j < len2; j++){	
			if ( *(*(v+i)+j) > aux_max){
				aux_max = *(*(v+i)+j);
			}
			if ( *(*(v+i)+j) < aux_min){
				aux_min = *(*(v+i)+j);
			}				
		}
		min[i] = aux_min;
		max[i] = aux_max;
	}
}

void get_array_minmax(double *v,int len, double *min, double *max){
	double aux_min, aux_max;
	int i;
	aux_min = *(v+0);	// v[0]
	aux_max = *(v+0);	// v[0]
		
	for(i = 0;i < len; i++){
		if ( *(v+i) > aux_max){
			aux_max = *(v+i);
		}
		if ( *(v+i) < aux_min){
			aux_min = *(v+i);
		}				
	}
	*min = aux_min;
	*max = aux_max;	
}

void normalize_rows_matrix(double **v,double *min, double *max, int len1, int len2){
	double average, range;
	int i,j;
	for(i = 0;i < len1; i++){
		average = (max[i] + min[i])/2;
		range = (max[i] - min[i])/2;
		for(j = 0;j < len2; j++){	
				*(*(v+i)+j)= (*(*(v+i)+j)-average)/range;	
		}

	}
}

void normalize_array(double *v, int len, double min, double max){
	double average, range;
	int i;
	average = (max + min)/2;
	range = (max - min)/2;
	for(i = 0;i < len; i++){	
		*(v+i)= (*(v+i)-average)/range;	
	}
}

void desnormalize_array(double *v, int len, double min, double max){
	double average, range;
	int i;
	average = (max + min)/2;
	range = (max - min)/2;
	for(i = 0;i < len; i++){	
		*(v+i) = (*(v+i))*range + average;
	}
}


