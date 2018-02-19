
#ifndef ELMHEADER_H_
#define ELMHEADER_H_

#include "./headers/basic_headers.h"	
#include "./headers/gsl_headers.h"	

#include "./headers/aux_func.h"	
#include "./headers/bitvector_func.h"	
#include "./headers/ELM.h"	
#include "./headers/GA.h"	
#include "./headers/SA.h"	
#include "./headers/AIS.h"	
#include "./headers/FP.h"	
#include "./headers/general_func.h"	
#include "./headers/output_func.h"
#include "./headers/threads.h"	
#include "./headers/ES.h"
#include "./headers/gsl_func.h"
#include "./headers/process_param.h"
#include "./headers/Graph.h"

//-------------> CONSTANTS <--------------------

#define DO_ELMS				0
#define TIME_ERMS_GRAPH		1
#define GA_ALGORITHM		2
#define EXHAUSTIVE_SEARCH	3
#define SA_ALGORITHM		4
#define AIS_ALGORITHM		5
#define FP_ALGORITHM		6

#define 	MAX_CHAR		100
#define MAX_UNSIGNED    4294967295	// (1 << 32) - 1

#define PLOT_DATA  1

//---------> Main extern variables <------
extern int	n;				// Number of parameters
extern int n_var_chosen;	// Number of unsigned int variables needed for storing the parametrs
extern ELM_params ELM_p;	// ELM structure to do standard ELM


#endif /* ELMHEADER_H_ */
