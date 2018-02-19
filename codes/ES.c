
#include "../ELMheader.h"

int ES_op (ES_data * ES_p ){
	unsigned int i;
	
	unsigned int * ES_selected = ES_p->ES_selected;
	unsigned int * ES_fixed_1s = ES_p->ES_fixed_1s;
	int n_TOP = ES_p->n_TOP;
	
	unsigned int * ES_1s_pos;	// Position of de 1s in the selected vector
	
	unsigned int *ES_order;			// Used for ordering the exaustive search order.
	unsigned int **ES_chosen_vectors = ES_p->ES_bitvectors;
	
	unsigned int 	n_ES;			// Number of different bits of the Exaustive search	
	

	ES_chosen_vectors = get_ES_vectors (ES_fixed_1s, ES_selected, 1 ,n ,&n_ES);	// Get the vectors of the search
													// n_ES will be modified by the func to have the number of them.
	ES_p->n_ES = n_ES;
	
printf("n_ES %i \n", n_ES);
	
	ES_p->ES_ERMS = (double *) malloc (sizeof(double)*n_ES);
	ES_order = (unsigned int *) malloc (sizeof(unsigned int)*n_ES);	
	
//----------> Release the CRAKEN !!!... I mean threads <--------

	threads_WorkOut(ES_chosen_vectors,ES_p->ES_ERMS, n_ES, &Thread_p);
	
printf("Exahustive Search Done \n");

	plot_ES_Search_Space (ES_p);
//----------> Get the best values and print them out <------------
	Thread_p.ELM_p_thread[0].FLAGS = PLOT_DATA_F;
	order_double (ES_p->ES_ERMS,(int *) ES_order,n_ES);
	ES_1s_pos =  get_1s_postions(ES_selected,  n);
/*	
	for (i = 0; i < n_TOP; i++){
		printf("La ERMS: %lf  \n", ES_p->ES_ERMS[(n_ES-1) -i]);
	
		get_ES_selected(ES_1s_pos, Thread_p.ELM_p_thread[0].chosen_vector, ES_order[(n_ES-1) -i], n);
		or_bitvector(Thread_p.ELM_p_thread[0].chosen_vector, ES_fixed_1s, n);	
		print_bitvector(Thread_p.ELM_p_thread[0].chosen_vector, n);
		
		Thread_p.ELM_p_thread[0].n_chosen = get_vector_weight(Thread_p.ELM_p_thread[0].chosen_vector, n);
		
		get_ELM_RMSE(&(Thread_p.ELM_p_thread[0]));
	}
*/
	//---------> Free operation memmory <-------
	free(ES_p->ES_ERMS);
	free(ES_order);	
		
	free(ES_selected);
	free( ES_fixed_1s); 
	
	for (i = 0; i < n_ES; i++){
		free(ES_chosen_vectors[i]);	
	}
	free(ES_chosen_vectors);
	free(ES_1s_pos);
	
	
	
	return 1;
}

int get_ES_selected(unsigned int *ones_pos, unsigned int *v, unsigned int n, int num_bits){
	int i;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	
	int num_ones;

	for (i = 0; i < n_var_chosen; i++){
		v[i] = 0;
	}
	num_ones = get_vector_weight(&n, num_bits);	// Get the number of bits of the output.
	i = 0;
	while( i  < 32 ){
		if (n & (1 << i)){			// If we have to set a bit
			var_chosen = ones_pos[i]/(sizeof(unsigned int)*8);
			bit_chosen = ones_pos[i]%(sizeof(unsigned int)*8);
			v[var_chosen] |= 1 << bit_chosen;
		}
		i++;
	}
	return 1;
}

unsigned int ** get_ES_vectors (unsigned int *ES_fixed_1s, unsigned int *ES_selected, unsigned int min1s,unsigned int n,unsigned int *num_vectors){
	unsigned int i, aux_int;
	unsigned int n_selected_1s = get_vector_weight(ES_selected, n);
	
	unsigned int ** ES_bitvectors;
	unsigned int * ES_1s_pos;
	unsigned int * ES_go;		// Vector with the 1s in right positions.
	
	int n_var =  roof_int (n, sizeof(unsigned int)*8);
	int n_ES_go_var;

	n_ES_go_var = roof_int (n_selected_1s, sizeof(unsigned int)*8);
	ES_go =  (unsigned int *) malloc (sizeof(unsigned int )*n_ES_go_var );	
	
	//*num_vectors = math_combination(n_selected_1s,n_selected_1s - min1s);
	*num_vectors = 1 << n_selected_1s;
	ES_bitvectors =  (unsigned int **) malloc (sizeof(unsigned int *)*(*num_vectors));	
	for (i = 0; i < (*num_vectors); i++){
		ES_bitvectors[i] = (unsigned int *) malloc (sizeof(unsigned int )*n_var);
	}
	
//printf("Numero de vectors: %i \n",*num_vectors);
//	for (i = 0; i < n; i++){	// ELIMINAR DESPUES
//		ES_selected[0] |= 1 << i;
//	}
	// ES_selected[0] |= 1 << i
	
print_bitvector(ES_selected, n);
print_bitvector(ES_fixed_1s, n);

	ES_1s_pos = get_1s_postions(ES_selected, n);

for (i = 0; i < n_selected_1s; i++){
	printf("%u ",ES_1s_pos[i]);
}
printf("\n");

	//----------> Get chosen vectors of the exhaustive search <------------
	aux_int = 1 << n_selected_1s;
	for ( i = 0; i < aux_int; i++){
		get_ES_selected(ES_1s_pos, ES_bitvectors[i],bin_to_gray (ES_go[0]), n);	// Get the chosen vector

		if (get_vector_weight(ES_bitvectors[i] , n)>= min1s){	// Check it has the proper number of bits
						// ES_ERMS[i] = get_ELM_RMSE(&ELM_p_thread[0]);			// ELM !!!!x
		}else {
						//ES_ERMS[i] = RAND_MAX;
		}		
		ES_go[0]++;
		
		or_bitvector(ES_bitvectors[i] , ES_fixed_1s, n);		// ADD UP THE FIXED ONES !!!!!
		
	print_bitvector(ES_bitvectors[i], n);
	}	
	
	
	free(ES_1s_pos);
	free(ES_go);	
	
	return ES_bitvectors;
}

int plot_ES_Search_Space (ES_data * ES_p){
	int i;
	multiple_graph mgra;
	gsl_vector **graph_ES;
	
	graph_ES = (gsl_vector **)malloc(2*sizeof(gsl_vector *));
	for (i = 0; i < 2; i++){
		graph_ES[i] = gsl_vector_alloc(ES_p->n_ES);
	}
	for (i = 0; i < ES_p->n_ES; i++){
		graph_ES[0]->data[i] = i;
	}		
	memcpy(graph_ES[1]->data, ES_p->ES_ERMS, ES_p->n_ES*sizeof(double));
	
	printf("colca \n");
	mgra.y = graph_ES;
	mgra.n_curves = 1;
	mgra.graph_name = (char *)malloc(30*sizeof(char));
	mgra.curve_names= (char **)malloc(1*sizeof(char*));
	
	mgra.X_axis_n = (char *)malloc(30*sizeof(char));
	mgra.Y_axis_n = (char *)malloc(30*sizeof(char));
	
	for (i = 0; i < 1; i++){
		mgra.curve_names[i] = (char *)malloc(30*sizeof(char));
	}
	sprintf(mgra.curve_names[0], "Search Space");
	
	sprintf(mgra.graph_name, "Search Space");
	
	sprintf(mgra.X_axis_n , "Solution");
	sprintf(mgra.Y_axis_n , "Cost function");

	plot_multiple_graph(&mgra);
	
}


