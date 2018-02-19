
#ifndef GSLHEADER_H_
#define GSLHEADER_H_
// ------------> Matrix functions <---------
void gsl_matrix_subvector_col_get(gsl_matrix *M, gsl_vector* v, int num_col, int offset, int size);
int gsl_matrix_subvector_row_get(gsl_matrix *M, gsl_vector* v, int num_row, int offset, int size);
int gsl_matrix_subvector_row_set(gsl_matrix *M, gsl_vector* v, int num_row, int offset, int size);
int gsl_matrix_subvector_col_set(gsl_matrix *M, gsl_vector* v, int num_col, int offset, int size);

int gsl_matrix_submatrix_cpy(gsl_matrix *M, gsl_matrix *m, int offrow, int offcol, int sizerow, int sizecol);
int gsl_inv_matrix(gsl_matrix *M,gsl_matrix *Minv );
int gsl_vector_dotsum(gsl_vector * a, gsl_vector * b);
int gsl_vector_dotsub(gsl_vector * a, gsl_vector * b);
int gsl_matrix_diag(gsl_matrix * m, gsl_vector * d);
int gsl_vector_scale_div (gsl_vector *v, double div);

void print_gsl_vector(gsl_vector *v);
void print_gsl_matrix(gsl_matrix *M);
double sum_gsl_vector(gsl_vector *x);
double sum_abs_gsl_vector(gsl_vector *x);
double rmse_gsl_vector(gsl_vector *x,gsl_vector *y);

#endif /* ELMHEADER_H_ */
