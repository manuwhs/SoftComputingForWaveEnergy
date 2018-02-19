#ifndef AUX_FUNCHEADER_H_
#define AUX_FUNCHEADER_H_

// ------------> Aux functions

void print_array(double *v, int len);
void print_matrix(double **v, int len1, int len2);

int roof_int (int num, int den);
unsigned int min_uint(unsigned int a, unsigned int b);
int order_double (double * lista, int * orden, int len);
void zero_char (char * array, int length); 
void zero_int (int * array, int length); 

int max(int a, int b);
float maxf (float * array, int num);  /* Give it the array and will give the best */
double min_d(double a, double b);

int dec2str(char  * cad, int digit);
int str2dec(char  * cad, char digit); 

int count_char(char * array, int len, int value);
int findint(int * array, int len, int value);

int ordenar_int (int * lista, int * orden, int len);
int ordenar_float (float * lista, int * orden, int len);

int reordenar_int (int * lista, int * orden, int len);
int reordenar_float (float * lista, int * orden, int len);

int desreordenar_char (char * lista, int * orden, int len); 

int copy_vector_char(char * duplicado, char * original, int num);
int copy_vector_int(int * duplicado, int * original, int num);
int copy_vector_float(float * duplicado, float * original, int num);

double abs_d (double d);

// -------> Normalization <---------
void get_matrix_row_minmax(double **v,double *min, double *max, int len1, int len2);
void get_array_minmax(double *v,int len, double *min, double *max);

void normalize_rows_matrix(double **v,double *min, double *max, int len1, int len2);
void normalize_array(double *v, int len, double min, double max);

void desnormalize_array(double *v, int len, double min, double max);

#endif 
