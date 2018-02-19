

#ifndef GRAPHHEADER_H_
#define GRAPHFUNCHEADER_H_

// ------------> Graph functions

int get_ERMS_TIME_Nparam_graph(ELM_params *p, linear_graph *l);
int get_ERMS_TIME_Ntrain_graph(ELM_params *p, linear_graph *l);
int Graph_op (linear_graph *l_graph, linear_graph *l_graph2, unsigned int * Graph_selected);
int print_lgraph_values(linear_graph *l_graph);
#endif /* ELMHEADER_H_ */
