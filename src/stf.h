/* returns tinti source time function based on analytical derivation 
   Note: all allocations done in main program.  the length of stf_out, should
   be equal to the length of t. 

   inputs
   		- t : array of times to compute the source time function
		- nt : length of t
		- Tr : Rise time, equal to the duration of the yoffe portion of the function
		- Ts : Half-width of the triangle window function
		- t0 : initial time of the source function
		- slip : total slip

	outputs:
		- stf_out : pointer to array of length nt.  (allocations done in main program to prevent messy
												     memory leaks)
*/
float* tinti(float* t, int nt, float Tr, float Ts, float t0, float slip);
float tinti_scalar(float t, int nt, float Tr, float Ts, float t0, float slip);

/* returns exponential source time function at times t. */
// float* exponential(float* t, float vpeak, float t0, float slip);

/* returns draeger source time function */
// float* draeger(float* t, float xi, float tau, float n, float slip);
