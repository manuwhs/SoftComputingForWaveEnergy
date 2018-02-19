#include "../ELMheader.h"
#include "../headers/threads.h"	

void* get_ELM_RMSE_thread(void *p){		// Thread function 
	 
	int doing_job;							// Job that the thread will do at every iteration
	ELM_params *p_ELM = (ELM_params *)p;  	// Give type to the pointer
	
	while (1) {
		pthread_mutex_lock (&(Thread_p.job_access)); 			// Mutex block
		if (Thread_p.thread_job < Thread_p.total_thread_jobs){	// Get any remaing job
			doing_job = Thread_p.thread_job;		// Use the job number as index for arrays
			Thread_p.thread_job++;
		} else {
			doing_job = -1;
		}
		pthread_mutex_unlock (&(Thread_p.job_access)); 	// Mutex unblock		
		
		if (doing_job == -1){
			break;										// Exit the while(1)
		}
		p_ELM->chosen_vector = Thread_p.Thread_chosen_array[doing_job];	// Get the job
		p_ELM->n_chosen = get_vector_weight(p_ELM->chosen_vector, n);	// Get the number of inputs.
		Thread_p.Thread_results[doing_job] = get_ELM_RMSE(p_ELM);		//-------->Do the ELM job<----------
	}
	pthread_exit(NULL); // Exit thread 
}

// This function does the jobs programmed in Thread_chosen_array and saves them into Thread_results
int threads_WorkOut(unsigned int **Bitvectors, double *Results, int num_jobs, Thread_params *Thread_p){
	int i;
	Thread_p->Thread_chosen_array = Bitvectors;	
	Thread_p->Thread_results = Results;				
	
	Thread_p->total_thread_jobs = num_jobs;		// Indica the number of ELM we have to do
	Thread_p->thread_job = 0;					//Initialize the jobs to 0
	
	
	for(i = 0; i < Thread_p->n_threads; i++){ 		// Create the threads
		pthread_create (&(Thread_p->threads_ID[i]), NULL, &get_ELM_RMSE_thread, &(Thread_p->ELM_p_thread[i]));	// Create thread 
	}

	for(i = 0; i < Thread_p->n_threads; i++){ 		// Wait for the threads to finish
		pthread_join (Thread_p->threads_ID[i], NULL); 
	}	
	return Thread_p->total_thread_jobs ;
}

int init_threads(ELM_params *ELM_p, Thread_params *Thread_p){
	int i;
	//------------> Initialize structures <--------------
	for (i = 0; i < Thread_p->n_threads; i++){
			Thread_p->ELM_p_thread[i].Nh = 			ELM_p->Nh;
			Thread_p->ELM_p_thread[i].n_ELMs = 			ELM_p->n_ELMs;
			Thread_p->ELM_p_thread[i].activation_f = 	ELM_p->activation_f;
			Thread_p->ELM_p_thread[i].FLAGS = 		ELM_p->FLAGS; 
				
			Thread_p->ELM_p_thread[i].x_train_param = ELM_p->x_train_param;
			Thread_p->ELM_p_thread[i].x_test_param =  ELM_p->x_test_param;
			Thread_p->ELM_p_thread[i].y_train = 		ELM_p->y_train;
			Thread_p->ELM_p_thread[i].y_test = 		ELM_p->y_test  ;
				
			Thread_p->ELM_p_thread[i].t_max = 		ELM_p->t_max ;
			Thread_p->ELM_p_thread[i].t_min = 		ELM_p->t_min;
			// "The chosen_vector" will be loaded into the structure from the "thread function"
			// The chosen vector for every ELM must be already defined in the global variable "unsigned int **Thread_chosen_array"
			// The function will use its "job_number" as an index for that structure and get it.
			// Usually we compute all the chosen vectors into another unsigned int **array and then do Thread_chosen_array=array
		}
	return 1;
}

int set_up_threads(unsigned int num, Thread_params *Thread_p){
	Thread_p->n_threads = num;
	Thread_p->threads_ID = (pthread_t *)malloc(Thread_p->n_threads*sizeof(pthread_t));	// Reserve memmory for the threads ID
	if (Thread_p->threads_ID == NULL) {
		perror("threads_ID");
		exit(-1);
	}	
		
	Thread_p->ELM_p_thread = (ELM_params *)malloc(Thread_p->n_threads*sizeof(ELM_params));	// Reserve memmory for the threads ID
	if (Thread_p->ELM_p_thread == NULL) {
		perror("ELM_p_thread");
		exit(-1);
	}			
	
	return 1;
}

int destroy_threads(Thread_params *Thread_p){
	free(Thread_p->threads_ID);		// Memmory for the threads	
	free(Thread_p->ELM_p_thread);
	return 1;
}

/* Using the threads create an appatent leak of memmory but its not
 * 
==826==    by 0x400AE1F: _dl_new_object (dl-object.c:77)
==826==    by 0x4006447: _dl_map_object_from_fd (dl-load.c:1051)
==826==    by 0x40083AF: _dl_map_object (dl-load.c:2568)
==826==    by 0x4012D6C: dl_open_worker (dl-open.c:225)
==826==    by 0x400ECCE: _dl_catch_error (dl-error.c:178)
==826==    by 0x441BE20: do_dlopen (dl-libc.c:89)
==826==    by 0x42DED4B: start_thread (pthread_create.c:308)
==826==    by 0x43E2BAD: clone (clone.S:130)
* 
==826== LEAK SUMMARY:
==826==    definitely lost: 0 bytes in 0 blocks
==826==    indirectly lost: 0 bytes in 0 blocks
==826==      possibly lost: 0 bytes in 0 blocks
==826==    still reachable: 954 bytes in 4 blocks
==826==         suppressed: 0 bytes in 0 blocks
* */







