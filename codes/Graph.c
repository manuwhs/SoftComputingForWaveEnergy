/*
 * functions.c
 *
 *  Created on: Nov 3, 2013
 *      Author: alvaro
 */

#include "../ELMheader.h"



int get_ERMS_TIME_Ntrain_graph(ELM_params *p, linear_graph *l){
	int i,j,k,aux_i;						// Auxiliary variables
	double ERMS_aux;						// For getting ERMS values
	struct timeval time_start,time_end;		// Time variables for general use
	float time_aux;							// For getting time
	
	char **curve_name;			// Array where to put the curve names
	char X_axis_t[MAX_CHAR] = "Hidden Neurons";
	char Y_axis_t[MAX_CHAR] = "Time (s)";
	char g_name_t[MAX_CHAR] = "Time";
	
	char X_axis_e[MAX_CHAR] = "Hidden Neurons";
	char Y_axis_e[MAX_CHAR] = "Ntrain_ERMS";
	char g_name_e[MAX_CHAR] = "Ntrain_ERMS";
	
	int init_Nh = l->init_Xvalue;		// Initial number of neurons
	int end_Nh = l->end_Xvalue;			// Final number of neurons
	int Nh_stepsize = l->XstepSize;		// Stepsize between number of neurons
	
	int init_Ntrain = l->init_Yvalue;		// Initial number of Ntrain vectors used
	int end_Ntrain = l->end_Yvalue;			// Final number of Ntrain vectors used
	int Ntrain_stepsize;				// Stepsize between number of Ntrain vectors

	int n_curves = l->Yn_curves;	// Number of curves we will draw	
	int v_length = (int)((end_Nh - init_Nh)/Nh_stepsize) + 1;			// Number of points of those curves
	
	int n_average = l->n_average; 
	
	multiple_graph graph_p;			//Parameters of the graph
	
	// Array of vectors for the graph data, the first vector is the X axis
	gsl_vector *time_v[n_curves + 1];	// Vector where to place the time vectors
										// +1 to store de X axis values
	gsl_vector *ERMS_v[n_curves + 1];
	
	if (l->Yn_curves > 1) 	// Get Ntrain StepSize
		Ntrain_stepsize = ((end_Ntrain - init_Ntrain)/(l->Yn_curves - 1));				// Stepsize between number of Ntrain vectors
	else
		Ntrain_stepsize = 0;
	
	for (i = 0; i < n_curves + 1; i++){
		time_v[i] = gsl_vector_alloc(v_length);
		ERMS_v[i] = gsl_vector_alloc(v_length);
	}
	
	//----> Create the names of the curves <------
	curve_name = (char **) malloc (sizeof(char *)*n_curves);
	
	for (i = 0; i < n_curves; i++){
		curve_name[i] = (char *) malloc (sizeof(char)*50);
		sprintf(curve_name[i],"%i", init_Ntrain + Ntrain_stepsize*i);
	}
	
	//------------>  GET CURVES DATA <------------
//p->FLAGS = 0; 			// Eliminate special stuff
	for (j = 0; j < n_curves; j++){		// For every train vector size
	// ***************************************
	//  CHANGE THE X_param vectors ->size  !!!!!
	// *********************************
			p->x_train_param[0]->size = init_Ntrain + Ntrain_stepsize*j;
			p->y_train->size = init_Ntrain + Ntrain_stepsize*j;
			
	//*********** FOR MAKING IT WITH NUMBER OF PARAMS ************
	//p->n_chosen = init_Ntrain + Ntrain_stepsize*j;

// printf("Elementos %i \n",v_length); 
		aux_i = 0;
		for (i = init_Nh; i <= end_Nh; i+= Nh_stepsize){	// For every neuron value we have to try
			p->Nh = i;	
			ERMS_aux = 0;
			time_aux = 0;
			for(k = 0; k < n_average; k++){			// For every averaging time
				gettimeofday( &time_start, NULL);
				
				ERMS_aux += get_ELM_RMSE(p);		// ELM !!!
				
				gettimeofday( &time_end, NULL);
				time_aux += get_time_passed(time_start,time_end);
			}
			ERMS_aux = ERMS_aux/n_average;
			time_aux = time_aux /n_average;
// printf("Accediendo al %i \n",aux_i); 
			gsl_vector_set(ERMS_v[j+1],aux_i, ERMS_aux);
			gsl_vector_set(time_v[j+1],aux_i, (double) time_aux);
			aux_i++;
		}
	}
	
	// Set the X axis of the curves
	aux_i = 0;
	for (i = init_Nh; i <= end_Nh; i+= Nh_stepsize){
		gsl_vector_set(ERMS_v[0],aux_i, i);
		gsl_vector_set(time_v[0],aux_i, i);
		aux_i++;
	}

//for (i = 0;  i < n_curves +1; i++){
//	print_gsl_vector(time_v[i]);
//}
	// -------------> Plot data <------------
	graph_p.y = time_v;
	graph_p.n_curves = n_curves;
	graph_p.curve_names = curve_name;
	graph_p.X_axis_n = X_axis_t;
	graph_p.Y_axis_n = Y_axis_t;
	graph_p.graph_name = g_name_t;
	
	plot_multiple_graph (&graph_p);
	
	
	graph_p.y = ERMS_v;
	graph_p.n_curves = n_curves;
	graph_p.curve_names = curve_name;
	graph_p.X_axis_n = X_axis_e;
	graph_p.Y_axis_n = Y_axis_e;
	graph_p.graph_name = g_name_e;
	
	plot_multiple_graph (&graph_p);
	
	
	//-------> Free dynamic memmory <----------
	for (i = 0; i < n_curves + 1; i++){
		gsl_vector_free(time_v[i]);
		gsl_vector_free(ERMS_v[i]);
	}
	
	for (i = 0; i < n_curves; i++){
		free(curve_name[i]);
	}
	free(curve_name);
	return 1;
}

int get_ERMS_TIME_Nparam_graph(ELM_params *p, linear_graph *l){
	int i,j,k,aux_i;
	double ERMS_aux;		// For getting ERMS values
	struct timeval time_start,time_end;		// Time variables for general use
	float time_aux;		// For getting time

	char **curve_name;			// Array where to put the curve names
	char X_axis_t[MAX_CHAR] = "Hidden Neurons";
	char Y_axis_t[MAX_CHAR] = "Time (s)";
	char g_name_t[MAX_CHAR] = "Nparam_Time";
	
	char X_axis_e[MAX_CHAR] = "Hidden Neurons";
	char Y_axis_e[MAX_CHAR] = "ERMS";
	char g_name_e[MAX_CHAR] = "Nparam_ERMS";
	
	int init_Nh = l->init_Xvalue;		// Initial number of neurons
	int end_Nh = l->end_Xvalue;			// Final number of neurons
	int Nh_stepsize = l->XstepSize;		// Stepsize between number of neurons
	
	int init_Nparam = l->init_Yvalue;		// Initial number of Ntrain vectors used
	int end_Nparam = l->end_Yvalue;			// Final number of Ntrain vectors used
	int Nparam_stepsize;				// Stepsize between number of Ntrain vectors

	int n_curves = l->Yn_curves;	// Number of curves we will draw	
	int v_length = (int)((end_Nh - init_Nh)/Nh_stepsize) + 1;			// Number of points of those curves
	
	int n_average = l->n_average; 
	
	multiple_graph graph_p;			//Parameters of the graph
	
	// Array of vectors for the graph data, the first vector is the X axis
	gsl_vector *time_v[n_curves + 1];	// Vector where to place the time vectors
										// +1 to store de X axis values
	gsl_vector *ERMS_v[n_curves + 1];
	
	if (l->Yn_curves > 1) 	// Get Ntrain StepSize
		Nparam_stepsize = ((end_Nparam - init_Nparam)/(l->Yn_curves - 1));				// Stepsize between number of Ntrain vectors
	else
		Nparam_stepsize = 0;
	
	for (i = 0; i < n_curves + 1; i++){
		time_v[i] = gsl_vector_alloc(v_length);
		ERMS_v[i] = gsl_vector_alloc(v_length);
	}
	
	//----> Create the names of the curves <------
	curve_name = (char **) malloc (sizeof(char *)*n_curves);
	
	for (i = 0; i < n_curves; i++){
		curve_name[i] = (char *) malloc (sizeof(char)*50);
		sprintf(curve_name[i],"%i", init_Nparam + Nparam_stepsize*i);
	}
	
	//------------>  GET CURVES DATA <------------
//p->FLAGS = 0; 			// Eliminate special stuff
	for (j = 0; j < n_curves; j++){		// For every train vector size
	// ***************************************
	//  CHANGE THE NUMBER OF PARAMS !!!!!
	// *********************************

	p->n_chosen = init_Nparam + Nparam_stepsize*j;

// printf("Elementos %i \n",p->n_chosen ); 
		aux_i = 0;
		for (i = init_Nh; i <= end_Nh; i+= Nh_stepsize){	// For every neuron value we have to try
			p->Nh = i;	
			ERMS_aux = 0;
			time_aux = 0;
//			 printf("Neurona %i \n",i ); 
			for(k = 0; k < n_average; k++){			// For every averaging time
				gettimeofday( &time_start, NULL);
				
				ERMS_aux += get_ELM_RMSE(p);		// ELM !!!
				
				gettimeofday( &time_end, NULL);
				time_aux += get_time_passed(time_start,time_end);
			}
			ERMS_aux = ERMS_aux/n_average;
			time_aux = time_aux /n_average;
// printf("Accediendo al %i \n",aux_i); 
			gsl_vector_set(ERMS_v[j+1],aux_i, ERMS_aux);
			gsl_vector_set(time_v[j+1],aux_i, (double) time_aux);
			aux_i++;
		}
	}
	
	// Set the X axis of the curves
	aux_i = 0;
	for (i = init_Nh; i <= end_Nh; i+= Nh_stepsize){
		gsl_vector_set(ERMS_v[0],aux_i, i);
		gsl_vector_set(time_v[0],aux_i, i);
		aux_i++;
	}

//for (i = 0;  i < n_curves +1; i++){
//	print_gsl_vector(time_v[i]);
//}
	// -------------> Plot data <------------
	graph_p.y = time_v;
	graph_p.n_curves = n_curves;
	graph_p.curve_names = curve_name;
	graph_p.X_axis_n = X_axis_t;
	graph_p.Y_axis_n = Y_axis_t;
	graph_p.graph_name = g_name_t;
	
	plot_multiple_graph (&graph_p);
	
	
	graph_p.y = ERMS_v;
	graph_p.n_curves = n_curves;
	graph_p.curve_names = curve_name;
	graph_p.X_axis_n = X_axis_e;
	graph_p.Y_axis_n = Y_axis_e;
	graph_p.graph_name = g_name_e;
	
	plot_multiple_graph (&graph_p);
	
	
	//-------> Free dynamic memmory <----------
	for (i = 0; i < n_curves + 1; i++){
		gsl_vector_free(time_v[i]);
		gsl_vector_free(ERMS_v[i]);
	}
	

	for (i = 0; i < n_curves; i++){
		free(curve_name[i]);
	}
	free(curve_name);
	return 1;
}

int Graph_op ( linear_graph *l_graph, linear_graph *l_graph2, unsigned int * Graph_selected){

		if (1){
printf("Doing TIME graph \n");
print_lgraph_values(l_graph);
			Thread_p.ELM_p_thread[0].chosen_vector = Graph_selected;
print_bitvector(Graph_selected,n);
			Thread_p.ELM_p_thread[0].n_chosen = get_vector_weight(Graph_selected,n);	
			get_ERMS_TIME_Ntrain_graph(&Thread_p.ELM_p_thread[0], l_graph);
		}
		if (0){
printf("Doing Parameters graph \n");
			Thread_p.ELM_p_thread[0].chosen_vector = (unsigned int *) malloc (sizeof(unsigned int)*n_var_chosen);
			fill_bitvector(Thread_p.ELM_p_thread[0].chosen_vector, n);
			Thread_p.ELM_p_thread[0].n_chosen = get_vector_weight(Thread_p.ELM_p_thread[0].chosen_vector,n);
			
//			ELM_p_thread[0].x_train_param[0]->size =  ;
//			ELM_p_thread[0].y_train->size = ;
			get_ERMS_TIME_Nparam_graph(&Thread_p.ELM_p_thread[0], l_graph2);
		}
		return 1;
}

int print_lgraph_values(linear_graph *l_graph){

	printf("Eje X -> From %i to %i in steps of length %i \n",l_graph->init_Xvalue, l_graph->end_Xvalue , l_graph->XstepSize );
	printf("Eje Y -> From %i to %i in %i curves\n",l_graph->init_Yvalue, l_graph->end_Yvalue , l_graph->Yn_curves );
	return 1;
}
