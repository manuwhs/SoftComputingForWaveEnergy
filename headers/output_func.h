

#ifndef OUTPUT_FUNCHEADER_H_
#define OUTPUT_FUNCHEADER_H_
typedef struct {
		gsl_vector ** y;		// GSL vectors with the data, y[0] = X-axis
		char ** curve_names;	// Name of the curves
		int n_curves;			// Number of the curves 
		
		char *X_axis_n;			// X-axis name
		char *Y_axis_n;			// X-axis name
		char *graph_name;		// name of the graph
	}multiple_graph;
	
typedef struct {		// For getting the values, not plotting them.
	
		int init_Xvalue;	// X-Axis parameters
		int end_Xvalue;
		int XstepSize;		
		
		int init_Yvalue;	// Y-Axis parameters
		int end_Yvalue;		
		int Yn_curves;
		int n_average;			// Number if times we run the simulation and average.
		
	}linear_graph;
			
// ------------> Output functions
void plot_multiple_graph (multiple_graph *p);
void plot(gsl_vector *x, gsl_vector *y);
void print_ELM_data (ELM_params *p);
void plot2(double *x, double *y, int size);
int check_file(char * subdir, char * file);
#endif /* ELMHEADER_H_ */
