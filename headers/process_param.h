
#ifndef PROCESS_PARAMHEADER_H_
#define PROCESS_PARAMELMHEADER_H_

// ------------> Get params functions
int get_ELM_params (char *config_file, Input_params *Input_p);
int get_GA_params (char *config_file, GA_data *GA_p);
int get_Graph_params (char *config_file, linear_graph *l_graph, linear_graph *l_graph2, unsigned int * Graph_selected);
int get_ELMs_params (char *config_file, ELMs_data *ELMs_p);
int get_ES_params (char *config_file, ES_data *ES_p);
int get_SA_params (char *config_file,SA_data *SA_p);
int get_AIS_params (char *config_file, AIS_data *AIS_p);
#endif /* ELMHEADER_H_ */
