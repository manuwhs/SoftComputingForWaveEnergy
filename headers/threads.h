
#ifndef THREADSHEADER_H_
#define THREADSHEADER_H_
#include "../ELMheader.h"
#include <pthread.h>	//Para los Hilos POSIX


typedef struct {
	int n_threads;					// Number of threads we can have simmultaneosly (max 30)
	int total_thread_jobs;			// Total number of thread jobs to do.
	int thread_job;					// Latest Job number taken so far
	
	ELM_params * ELM_p_thread;		// Structures with the data for an ELM go of a theread
	
	unsigned int **Thread_chosen_array;	// Array with the chosen vectors for every ELM thread job
	double *Thread_results;				// Array with the results of every thread.
	
	pthread_mutex_t job_access;		// Mutex for the access to the "thread_job" variable.
	pthread_mutex_t output_access;	// Mutex for the access to the output files.
	pthread_t* threads_ID;			// Array with the threads ID's 

} Thread_params;

extern Thread_params Thread_p;
int threads_WorkOut(unsigned int **Bitvectors, double *Results, int num_jobs,Thread_params *Thread_p );
void* get_ELM_RMSE_thread(void *p);	// Thread function 
int init_threads(ELM_params *ELM_p, Thread_params *Thread_p);
int set_up_threads(unsigned int num, Thread_params *Thread_p);
int destroy_threads(Thread_params *Thread_p);
#endif /* ELMHEADER_H_ */
