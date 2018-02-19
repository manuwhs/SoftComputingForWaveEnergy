

#ifndef GENERAL_FUNCHEADER_H_
#define GENERAL_FUNCHEADER_H_

// -------------------> General Functions <----------------
float get_time_passed(struct timeval time_start,struct timeval time_end);
gsl_matrix* get_selected_input_matrix(gsl_vector **param,unsigned int *chosen_param, int n_chosen);
void convert_data_to_gsl(double **X, double *t, int n, int Nsamples, gsl_vector ***Xgsl, gsl_vector **tgsl);
void shuffle_train_test(double **X_tr, double *t_tr, double **X_te, double *t_te, Input_params *Input_p );

void load_data(double ***X, double **t, int n, int Nsamples, char *Xdir, char *tdir);
void load_data2(double ***X, double **t, int n, int Nsamples, char *dir);
void load_data3(double ***X, double **t, int n, int Nsamples, char *dir);
void Super_Load(double ***X_tr, double **t_tr, double ***X_te, double **t_te, Input_params *Input_p);
void Super_Load2(double ***X_tr, double **t_tr, double ***X_te, double **t_te, Input_params *Input_p);
void Super_LoadX(double ***X_tr, double **t_tr, double ***X_te, double **t_te, Input_params *Input_p);

#endif /* ELMHEADER_H_ */
