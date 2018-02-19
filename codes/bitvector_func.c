/*
 * main.c
 *
 *  Created on: Nov 3, 2013
 *      Author: alvaro
 */

#include "../ELMheader.h"

// Generates a bit array !! MAYBE NEED TO BE REPLACED BY THE OTHER OLDER
int get_random_bitarray(unsigned int *v, int max_bits){	//Since 1 << 31 only has effect for "unsigned int"
	int i;
	int sign_bits;	// Since RAND_MAX = 2.147.483.647 = 31 LSB bits (lacking sign bit)
					// We will use the bit of this variable to set the sign bit
	int num_var = max_bits/(sizeof(unsigned int)*8);		// Number of int variables the bit array uses
	int rem_bits = max_bits%(sizeof(unsigned int)*8);	// Remaining bits	
													// We will convert it into a MASK
	unsigned int mask = 0;
	sign_bits = rand();
	
	if (rem_bits > 0 ) {
		num_var++;
		for (i = 0; i < num_var - 1; i++){
			v[i] = rand();
		}
		
		for (i = 0; i < rem_bits; i++) {
			mask |= 1 << i;
		}	
		v[num_var - 1] = 0;
		v[num_var - 1] = rand() & mask; 
			
	}
	else {
		for (i = 0; i < num_var; i++){
			v[i] = rand();
		}
		
	}
	return get_vector_weight(v, max_bits);
	
}

unsigned int hamming_distance (unsigned int *v, unsigned int *w, int max_length){
	int i;
	int hamm_d = 0;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	

	for (i = 0; i  < max_length ; i++){
		var_chosen = i/(sizeof(unsigned int)*8);
		bit_chosen = i%(sizeof(unsigned int)*8);	
		
		if ((v[var_chosen] & (1 << bit_chosen)) != (w[var_chosen] & (1 << bit_chosen))){
			hamm_d++;
		}
	}
//printf("Distancia: %i \n", hamm_d);
	return hamm_d;
}

unsigned int bitvector_setbit(unsigned int *v, int bit){
	int i;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	

	var_chosen = bit/(sizeof(unsigned int)*8);
	bit_chosen = bit%(sizeof(unsigned int)*8);	
	v[var_chosen] |= (1 << bit_chosen);

	return bit;
}

unsigned int bitvector_clearbit(unsigned int *v, int bit){
	int i;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	

	var_chosen = bit/(sizeof(unsigned int)*8);
	bit_chosen = bit%(sizeof(unsigned int)*8);	
	v[var_chosen] &= ~(1 << bit_chosen);

	return bit;
}

int bitvector_cmp(unsigned int *v, unsigned int *w){	//Since 1 << 31 only has effect for "unsigned int"
	unsigned int i;
	for (i = 0; i < n_var_chosen; i++){
		if (v[i] != w[i]){
			return -1;
		}
	}
	return 0;
}	

// Counts the number of vectors in an array of int
// Number maximum of bits 
int get_vector_weight(unsigned int *v, int length){
	int i;
	int num_ones = 0;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	
	for (i =0; i < length; i++){
		var_chosen = i/(sizeof(unsigned int)*8);
		bit_chosen = i%(sizeof(unsigned int)*8);
		if ((v[var_chosen] & (1 << bit_chosen)) > 0){
			num_ones++;
		}
	}
	return num_ones;
}

int print_bitvector(unsigned int *v, int max_length){
	int i;
	char bitarr[500];
	get_bitvector(v, max_length, bitarr);
	printf("%s \n", bitarr);
	return 1;
}

int get_bitvector(unsigned int *v, int max_length, char *bitarray){
	int i;
	int num_ones = 0;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	
	int	num_var = roof_int (max_length, sizeof(unsigned int)*8);
		
	num_ones = get_vector_weight(v, max_length);
//printf("Vector (%i 1's): ", num_ones);
	bitarray[max_length] = 0;
	for (i = 0; i  < max_length ; i++){
		var_chosen = i/(sizeof(unsigned int)*8);
		bit_chosen = i%(sizeof(unsigned int)*8);
		if (((v[var_chosen] >> bit_chosen) & 1 )> 0){
			bitarray[i] = '1';
		}
		else {
			bitarray[i] = '0';
		}
	}
	return 1;
}

int get_u_vector(unsigned int *v, int max_length, char *bitarray){
	int i,n;
	int var_chosen;			// Variables of the chosen param array 
	char bit_chosen;			// Chosen bit of that variable	
	int u_bits = sizeof(unsigned int)*8;
	int	num_var = roof_int (max_length, u_bits);	
	
	n = strlen(bitarray);
	for (i = 0; i < num_var; i++){
		v[i] = 0;
	}
	for (i = 0; i  < n ; i++){
		var_chosen = i/u_bits;
		bit_chosen = i%u_bits;
		if (bitarray[i] == '1'){
			v[var_chosen] |= (1 << bit_chosen);
		}
	}
	return 1;
}

int bitvector_cpy(unsigned int *d, unsigned int *s, int num_bits){
	int i;
	int num = roof_int (n, sizeof(unsigned int)*8);
	
	for (i = 0; i < num; i++){
		d[i] = s[i];
	}
	return 1;
}

int or_bitvector(unsigned int *d, unsigned int *s, int num_bits){
	int i;
	int num = roof_int (n, sizeof(unsigned int)*8);
	
	for (i = 0; i < num; i++){
		d[i] |= s[i];
	}
	return 1;
}

int zero_bitvector(unsigned int *v, int num_bits){
	int i;
	int num = roof_int (n, sizeof(unsigned int)*8);
	
	for (i = 0; i < num; i++){
		v[i] = 0;
	}
	return 1;
}

int fill_bitvector(unsigned int *v, int num_bits){
	int i;
	unsigned int var_chosen, bit_chosen;
	
	zero_bitvector(v, num_bits);
	for(i = 0; i < n; i++){ 
		var_chosen = i/(sizeof(unsigned int)*8);			// Variables of the chosen param array 
		bit_chosen = i%(sizeof(unsigned int)*8);			// Chosen bit of that variable		
		v[var_chosen] |= 1 << bit_chosen;	
	}


	
	return 1;
}

int mutate_bitvector(unsigned int *v, int n_bits, int n_mutations){
	int i;
	int prob;
	int var_chosen;
	int bit_chosen;	
	
	for(i = 0; i < n_mutations; i++){
			prob = (rand() % n_bits);
			var_chosen = prob/(sizeof(unsigned int)*8);
			bit_chosen = prob%(sizeof(unsigned int)*8);	
			v[var_chosen] ^= 1 << bit_chosen ;	//Mutate one bit
	}
}

int and_bitvector(unsigned int *d, unsigned int *s, int num_bits){
	int i;
	int num = roof_int (n, sizeof(unsigned int)*8);
	
	for (i = 0; i < num; i++){
		d[i] &= s[i];
	}
	return 1;
}

unsigned int * get_1s_postions(unsigned int * vector, unsigned int n){
	unsigned int i,j;
	unsigned int var_chosen, bit_chosen;
	unsigned int * position_1s = (unsigned int *) malloc (sizeof(unsigned int)*get_vector_weight(vector, n));
	//-----> Get the position of the 1s <---------
	i = 0;
	j = 0;
	while( j < n){
		var_chosen = j/(sizeof(unsigned int)*8);			// Variables of the chosen param array 
		bit_chosen = j%(sizeof(unsigned int)*8);			// Chosen bit of that variable	
		if (vector[var_chosen] & (1 << bit_chosen)){
			position_1s[i] = var_chosen*sizeof(unsigned int)*8 + bit_chosen;
			i++;
		}
		j++;
	}

	return position_1s;
}

int get_ubits(unsigned int v, int max_length, char *bitarray){
	int i;
	for (i = 0; i  < max_length ; i++){
		if (((v >> i) & 1 )> 0){
			bitarray[( max_length -1) - i] = '1';
		} else {
			bitarray[( max_length -1) - i] = '0';
		}
	}	
	bitarray[max_length] = 0;
}

unsigned int bin_to_gray (unsigned int v){
	return ((v ^ (v >> 1)));
}



