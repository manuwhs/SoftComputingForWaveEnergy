/*
 * functions.c
 *
 *  Created on: Nov 3, 2013
 *      Author: alvaro
 */

#include "../ELMheader.h"

int check_file(char * subdir, char * file) {
	FILE *pf;
	char ls[MAX_CHAR] = "ls";
	char tmp_file[MAX_CHAR] = "ls_list"; 
	char command[100];
	char buffer[100];
	//---------> Check the main folder exists <------------- 
	sprintf(command,"%s ", "ls");
	sprintf(command + strlen(command)," %s > %s",subdir, tmp_file);
	system(command);	
//printf("Comando %s \n",command);
    pf = fopen(tmp_file,"r");
    if (pf == NULL){
		perror ("fopen");
		printf("Error opening %s \n", tmp_file);
	}
	while(!feof(pf)){
		fscanf(pf,"%s",buffer);	
//printf("%s \n", buffer);
		if (strcmp(buffer, file) == 0){
			fclose(pf);	
//				printf("Salida 1 \n");
			return 1;
		}
	}
	fclose(pf);	
	//remove(tmp_file);
//printf("Salida 0 \n");
	return 0;
}

int make_folders(char * folder_dir){
	
	
	
}


void plot(gsl_vector *x, gsl_vector *y){
    FILE *pf;
    unsigned int k = 0;
    int nbytes = 0;
    double RMSE;
    char buffer[150];
    char command[100];
    char result_folder[MAX_CHAR] = "Results";
    char subfolder[MAX_CHAR] = "RMSE";
    char base_folder[300];
	char mkdir[MAX_CHAR] = "mkdir ";
	char remove[MAX_CHAR] = "rm ";
	char gnuplot[MAX_CHAR] = "gnuplot ";
	char tmp_gpfile[MAX_CHAR] = "tmpor.gp";
	char regression_file[MAX_CHAR] = "regression";
	char plot_file[MAX_CHAR] = "plot.png";
	
	char dir[60];
	char RMSE_s[30];

	RMSE = rmse_gsl_vector(x,y);	
	// Locking 
	pthread_mutex_lock (&(Thread_p.output_access));
//	printf("Entro\n");
	// Get the RMSE
	sprintf(RMSE_s, "RMSE:%lf", RMSE);
	
	//---------> Make the main folder <------------- 
	if (check_file("",result_folder) == 0){
		sprintf(command,"%s ", mkdir);
		sprintf(command + strlen(command),"%s", result_folder);
		system(command);
	}

//printf("%s \n", command);

	//-------------> Make the RMSE base folder <-----------------
	sprintf(base_folder,"%s", result_folder);
		
	if (check_file(base_folder,subfolder) == 0){	
		sprintf(command,"%s", mkdir);
		sprintf(command + strlen(command),"%s/%s", base_folder,subfolder);
		system(command);	
	}
	
	//-------------> Make the own  folder <-----------------
	sprintf(base_folder,"%s", result_folder);
	sprintf(base_folder + strlen(base_folder),"/%s", subfolder);	
		
	if (check_file(base_folder,RMSE_s) == 0){	
		sprintf(command,"%s", mkdir);
		sprintf(command + strlen(command),"%s/%s/%s", result_folder,subfolder, RMSE_s);
		system(command);	
	}
	//-------------> Make the base folder <-----------------	
		// Write real and predicted results into a file to plot them 
	sprintf(base_folder,"%s/%s/%s", result_folder,subfolder, RMSE_s);	
	sprintf(dir,"%s", base_folder);
	sprintf(dir + strlen(dir),"/%s", regression_file);

    pf = fopen(dir,"w");
    if (pf == NULL){
		perror ("fopen");
		printf("Error opening for writting %s \n", dir);
	}
    for(k = 0;k < x->size; k++){
        nbytes = sprintf(buffer," %.2f %.2f\n", gsl_vector_get(x,k), gsl_vector_get(y,k));
        fwrite(buffer,sizeof(char),nbytes,pf);
    }
    fclose(pf);
    
// Create temporary file to put the GNUPLOT commands

	sprintf(dir,"%s", base_folder);
	sprintf(dir + strlen(base_folder),"/%s", tmp_gpfile);
	
//printf("%s \n",dir);
    pf = fopen(dir,"w");
      if (pf == NULL){
		perror ("fopen");
		printf("Error opening %s \n", dir);
	}
    nbytes = sprintf(buffer,"set terminal png font Helvetica 16\n");
    fwrite(buffer,sizeof(char),nbytes,pf);  
    
    sprintf(dir,"%s", base_folder);
	sprintf(dir + strlen(base_folder),"/%s", plot_file);
    nbytes = sprintf(buffer,"set output '%s'\n",dir);
    fwrite(buffer,sizeof(char),nbytes,pf);   
    
    nbytes = sprintf(buffer,"set xlabel 'Real data'\n");
    fwrite(buffer,sizeof(char),nbytes,pf);
    
    nbytes = sprintf(buffer,"set ylabel 'Predicted data'\n");
    fwrite(buffer,sizeof(char),nbytes,pf);
  
    nbytes = sprintf(buffer,"set xrange [-0.50:] \n set yrange [-0.50:]\n");
    fwrite(buffer,sizeof(char),nbytes,pf);
    
    nbytes = sprintf(buffer,"set style line 1 lc rgb '#0000A0' pt 7 ps 1 \n");	// Text style
    fwrite(buffer,sizeof(char),nbytes,pf);  
  
    nbytes = sprintf(buffer,"set style line 2 lt 1 lc rgb 000000 lw 2 \n");	// Text style
    fwrite(buffer,sizeof(char),nbytes,pf);    
    
    nbytes = sprintf(buffer,"f(x) = x \n");	// Text style
    fwrite(buffer,sizeof(char),nbytes,pf);  

//     nbytes = sprintf(buffer,"set title 'RMSE: %f' left  \n", RMSE);
 //   fwrite(buffer,sizeof(char),nbytes,pf);     
    
    sprintf(dir,"%s", base_folder);
	sprintf(dir + strlen(base_folder),"/%s", regression_file);
    nbytes = sprintf(buffer,"plot '%s' ls 1 title 'RMSE: %f' , ",dir, RMSE);
    fwrite(buffer,sizeof(char),nbytes,pf);
    

    nbytes = sprintf(buffer," f(x) title 'Ideal' with lines ls 2, '%s' using 1:2:($0+2) notitle with dots ls 1   \n ",dir); 
    
    // with labels 
    fwrite(buffer,sizeof(char),nbytes,pf);
    
    fclose(pf);
    
	// Plot the graph
    sprintf(command,"%s", gnuplot);
    sprintf(dir,"%s", base_folder);
	sprintf(dir + strlen(dir),"/%s", tmp_gpfile);
	
	sprintf(command + strlen(command),"%s",dir);
	system(command);
    
    // Remove temporary file
    sprintf(command,"%s", remove);
    
    sprintf(dir,"%s", base_folder);
	sprintf(dir + strlen(base_folder),"/%s", tmp_gpfile);
	
	sprintf(command + strlen(command),"./%s", dir);
	// system(command);

	pthread_mutex_unlock (&(Thread_p.output_access));
//	printf("------------>\n");
}

void plot_multiple_graph (multiple_graph *p) {
	// y[0] y es X-axis and the others are the y-axises
	// num_curves is the number of curves
	// curve names is the name of each of those curves.
	
	// -----> Parameters<----------	
	gsl_vector ** y = p->y;		// GSL vectors with the data, y[0] = X-axis
	char ** curve_name = p->curve_names;	// Name of the curves
	int num_curves = p->n_curves;			// Number of the curves 
		
	char *X_axis_n = p->X_axis_n;			// X-axis name
	char *Y_axis_n = p->Y_axis_n;			// X-axis name
	char *graph_name = p->graph_name;		// name of the graph	
	
    FILE *pf = NULL;
    int i ,  nbytes = 0;
    unsigned int k = 0;
    char buffer[150];
	char curves_file[MAX_CHAR] = "curves.dat";
	char temp_file[MAX_CHAR] = "temp.dat";	
	
//	printf("Got parammeters **********************\n");
	
    pf = fopen(curves_file,"w");
    if (pf == NULL){
		perror ("fopen");
		printf("Error opening %s \n", curves_file);
	}    
    // Write the vectors into the file
    for(k = 0 ;k < y[0]->size; k++) {		// For every element of the vectors
		nbytes = 0;		// Number of bytes written.
		buffer[0] = 0;	// Empty  buffer
        for(i = 0 ;i < num_curves + 1; i++) {		// For every element of the vectors
			nbytes += sprintf(buffer + strlen(buffer),"%.5f ", gsl_vector_get(y[i],k));
		}
		nbytes += sprintf(buffer + strlen(buffer),"\n");
        fwrite(buffer,sizeof(char),nbytes,pf);
    }
    fclose(pf);

    pf = fopen(temp_file,"w");
    if (pf == NULL){
		perror ("fopen");
		printf("Error opening %s \n", temp_file);
	}
    nbytes = sprintf(buffer,"set terminal png font Helvetica 16\n");	// Text style
    fwrite(buffer,sizeof(char),nbytes,pf);
    
    nbytes = sprintf(buffer,"set grid\n");							// Grid
    fwrite(buffer,sizeof(char),nbytes,pf);
    
     nbytes = sprintf(buffer,"set xrange [] \n set yrange []\n");
    fwrite(buffer,sizeof(char),nbytes,pf);   
      
    nbytes = sprintf(buffer,"set output '%s.png'\n",graph_name);		// Output file
    fwrite(buffer,sizeof(char),nbytes,pf);
    
    nbytes = sprintf(buffer,"set xlabel '%s'\n",X_axis_n);		// X label
    fwrite(buffer,sizeof(char),nbytes,pf);
    
    nbytes = sprintf(buffer,"set ylabel '%s'\n",Y_axis_n);		// Y labe;
    fwrite(buffer,sizeof(char),nbytes,pf);
    
    nbytes = sprintf(buffer,"set style line 1 lw 2 \n");	// Text style
    fwrite(buffer,sizeof(char),nbytes,pf);     
 
    nbytes = sprintf(buffer,"set style line 2 lc rgb '#0000A0' lw 2 pt 6 ps 3 \n");	// Text style
    fwrite(buffer,sizeof(char),nbytes,pf);  
       
    // Print first graph
    nbytes = sprintf(buffer,"plot '%s' using 1:%i with lines ls 2  smooth unique title '%s'",curves_file,2,curve_name[0]);
    fwrite(buffer,sizeof(char),nbytes,pf);
    

    for (i = 0; i < num_curves - 1; i++){	// Give the gnuplot order to paint the graphs
		nbytes = sprintf(buffer,",'%s' u 1:%i w l ls 1 lc %i title '%s' ",curves_file,i + 3,i,curve_name[i+1]);
		fwrite(buffer,sizeof(char),nbytes,pf);
	}

     
    fclose(pf);
    sprintf(buffer, "gnuplot %s",temp_file);
    system(buffer);
    // remove(temp_file);

}

void print_ELM_data (ELM_params *p){
 	FILE *pf;
	int nbytes = 0;
	char buffer[150];
	char command[100];
	char result_folder[MAX_CHAR] = "Results";
	char data_file[MAX_CHAR] = "data";
	char subfolder[MAX_CHAR] = "RMSE";
	char mkdir[MAX_CHAR] = "mkdir ";
	char folder_name[50];		// Name of the folder with the files
	char base_folder[300];
	char aux_bitvector[n];
	char date[60];
	char dir[60];
	int i;
	time_t fechaActual ;
	struct tm * fechaPtr;
	
	char RMSE_s[30];
	sprintf(RMSE_s, "RMSE:%lf", p->ERMS);
	// Get the date.
	fechaActual = time(0) ;		// Get the time and print it
	fechaPtr = gmtime(&fechaActual) ;
	nbytes = sprintf(date, "%i-%i-%i(%i:%i:%i)", fechaPtr->tm_mday, fechaPtr->tm_mon + 1, 
		fechaPtr->tm_year + 1900, fechaPtr->tm_hour, fechaPtr->tm_min, fechaPtr->tm_sec);	
	
	//--------> Make results folder <-------------
	if (check_file("",result_folder) == 0){
		sprintf(command,"%s ", mkdir);
		sprintf(command + strlen(command),"%s", result_folder);
		system(command);
	}
	
	//-------------> Make the RMSE base folder <-----------------
	sprintf(base_folder,"%s", result_folder);
		
	if (check_file(base_folder,subfolder) == 0){	
		sprintf(command,"%s", mkdir);
		sprintf(command + strlen(command),"%s/%s", base_folder,subfolder);
		system(command);	
	}
	
	//-------------> Make the own  folder <-----------------
	sprintf(base_folder,"%s", result_folder);
	sprintf(base_folder + strlen(base_folder),"/%s", subfolder);	
		
	if (check_file(base_folder,RMSE_s) == 0){	
		sprintf(command,"%s", mkdir);
		sprintf(command + strlen(command),"%s/%s/%s", result_folder,subfolder, RMSE_s);
		system(command);	
	}
	//-------------> Make the base folder <-----------------	
		// Write real and predicted results into a file to plot them 
	sprintf(base_folder,"%s/%s/%s", result_folder,subfolder, RMSE_s);	
	sprintf(dir,"%s", base_folder);

	//-----------> Write data into file <-----------------

	sprintf(dir,"%s/%s", base_folder,data_file);
	
	pf = fopen(dir,"w");
	 if (pf == NULL){
		perror ("fopen");
		printf("Error opening %s \n", dir);
	}

	nbytes = sprintf(buffer,"Date: %s \n", date);		// Date
	fwrite(buffer,sizeof(char),nbytes,pf);
  
 	nbytes = sprintf(buffer,"RMSE: %lf \n", p->ERMS);	// RMSE
	fwrite(buffer,sizeof(char),nbytes,pf);
	
 	nbytes = sprintf(buffer,"Selected_among: %i_parameters \n\n", n);	// Num of parameters  CGANGE %%%%%%$%^$#
	fwrite(buffer,sizeof(char),nbytes,pf);	

	get_bitvector(p->chosen_vector, n, aux_bitvector); 
	
	nbytes = sprintf(buffer,"\nSelected_vector: %s\n\n",aux_bitvector);  // Selected parameters in "bitvector" form
	fwrite(buffer,sizeof(char),nbytes,pf);	   

	nbytes = sprintf(buffer,"Hidden_neurons: %i \n", p->Nh);	// Hidden Neurons
	fwrite(buffer,sizeof(char),nbytes,pf);
	
	nbytes = sprintf(buffer,"Train_vectors: %i \n", p->y_train->size);	// Train vectors number
	fwrite(buffer,sizeof(char),nbytes,pf);
	
	nbytes = sprintf(buffer,"Test_vectors: %i \n", p->y_test->size);	// Test vectors number
	fwrite(buffer,sizeof(char),nbytes,pf);
	
    fclose(pf);
    
    //-----------> Go back to main folder <-----------------
}

void plot2(double *x, double *y, int size){
    FILE *pf;
    unsigned int k = 0;
    int nbytes = 0;
    char buffer[150];
    char command[100];
    char result_folder[MAX_CHAR] = "Results";
    char base_folder[50];
	char mkdir[MAX_CHAR] = "mkdir ";
	char remove[MAX_CHAR] = "rm ";
	char gnuplot[MAX_CHAR] = "gnuplot ";
	char tmp_gpfile[MAX_CHAR] = "tmp.gp";
	char regression_file[MAX_CHAR] = "regression";
	char plot_file[MAX_CHAR] = "plot.png";
	
	char dir[60];


		
	//---------> Make the main folder <------------- 
	sprintf(command,"%s ", mkdir);
	sprintf(command + strlen(command),"%s", result_folder);
	system(command);
	
//printf("%s \n", command);

	//-------------> Make the base folder <-----------------
	sprintf(base_folder,"%s", result_folder);
	sprintf(base_folder + strlen(base_folder),"/%s", "pene");	
	
	sprintf(command,"%s", mkdir);
	sprintf(command + strlen(command),"%s", base_folder);
	system(command);
//printf("%s \n", command);
	
	// Write real and predicted results into a file to plot them 
	sprintf(dir,"%s", base_folder);
	sprintf(dir + strlen(dir),"/%s", regression_file);
	
//printf("%s \n",dir);
    pf = fopen(dir,"w");
    if (pf == NULL){
		perror ("fopen");
	}
    for(k = 0;k < size; k++){
        nbytes = sprintf(buffer," %.2f %.2f\n", x[k], y[k]);
        fwrite(buffer,sizeof(char),nbytes,pf);
    }
    fclose(pf);
    
// Create temporary file to put the GNUPLOT commands

	sprintf(dir,"%s", base_folder);
	sprintf(dir + strlen(base_folder),"/%s", tmp_gpfile);
	
//printf("%s \n",dir);
    pf = fopen(dir,"w");
    nbytes = sprintf(buffer,"set terminal png font Helvetica 16\n");
    fwrite(buffer,sizeof(char),nbytes,pf);  
    
    sprintf(dir,"%s", base_folder);
	sprintf(dir + strlen(base_folder),"/%s", plot_file);
    nbytes = sprintf(buffer,"set output '%s'\n",dir);
    fwrite(buffer,sizeof(char),nbytes,pf);   
    
    nbytes = sprintf(buffer,"set xlabel 'Real data'\n");
    fwrite(buffer,sizeof(char),nbytes,pf);
    
    nbytes = sprintf(buffer,"set ylabel 'Predicted data'\n");
    fwrite(buffer,sizeof(char),nbytes,pf);
    
    sprintf(dir,"%s", base_folder);
	sprintf(dir + strlen(base_folder),"/%s", regression_file);
    nbytes = sprintf(buffer,"plot '%s' pt 7, ",dir);
    fwrite(buffer,sizeof(char),nbytes,pf);
    
    nbytes = sprintf(buffer," '%s' using 1:2:($0+2) with labels offset 1 notitle",dir); 
    fwrite(buffer,sizeof(char),nbytes,pf);
    
    fclose(pf);
    
	// Plot the graph
    sprintf(command,"%s", gnuplot);
    sprintf(dir,"%s", base_folder);
	sprintf(dir + strlen(base_folder),"/%s", tmp_gpfile);
	
	sprintf(command + strlen(gnuplot),"%s",dir);
	system(command);
    
    // Remove temporary file
    sprintf(command,"%s", remove);
    
    sprintf(dir,"%s", base_folder);
	sprintf(dir + strlen(base_folder),"/%s", tmp_gpfile);
	
	sprintf(command + strlen(command),"./%s", dir);
	system(command);

}
