#include <stdbool.h>

/* creates 1d array of values starting with "start" ending with "end" and a
   length of "length". don't forget to free buffer in main when done! */
float* linspace(float start, float end, int length);

/* creates 1d array of values starting with "start" ending with "end" and a
   spacing of dx. don't forget to free buffer in main when done! */
float* arange(float start, float end, float dx, int *n);

/* creates zero initialized array of floats.
    don't forget to free buffer in main when done! */
float* zeros_f(int n);

/* creates zero initialized array of ints.
    don't forget to free buffer in main when done! */
int* zeros_i(int n);

/* computes sum of 1d array */
float sum(float *buf, int n);

/* computes cartesian product of jagged array A and returns all combinations in B. 
	inputs
		float **A : array containing all potential values of the parameters.
		int *np : length of A[i][]
		int ndim : number of dimensions
		int *ncom (by ref) : number of possible combinations
	returns
		float **B : size B[ncom][ndim] where ncom is the total number of combinations possible
*/
float** cartesian_product( float **A, int *np, int *ncom, int ndim );

float maximum( float* buf, int n );

/* 
	compares two floating point values to a value.  typically use eps = 1 / 2^20
*/
bool compare_f( float a, float b, float eps );

/* clean function to deallocate 2d float array */
void dealloc2d_f( float** buf, int n );

/* clean function to deallocate 2d int array */
void dealloc2d_i( int** buf, int n );

/* outer product obv, these vectors must be the same length */
float **outer( float *v1, float *v2, int n );

/* rastrigin function used for optimization defaults to 2d
	inputs
		(float) x : x value
		(float) y : y value
	returns
		(float) val : rastrigin @ (x,y)
*/
float rastrigin( float x, float y );

/* ackley function used for optimization defaults to 2d
	inputs
		(float) x : x value
		(float) y : y value
	returns
		(float) val : ackley @ (x,y)
*/
float ackley( float x, float y );

/* used to find index of first occurance of value 
	inputs
		(void *) buf : buffer to search through
		(int) n : length of buf
		(char) type : any of {">" "<" "==" "=!" "<=" ">="}
	returns (int) : location of first
*/
int find_first( float* buf, int n, char* type, float val );

/* used to find index of last occurance of value satisfying condition "type"
	inputs
		(void *) buf : buffer to search through
		(int) n : length of buf
		(char) type : any of {">" "<" "==" "=!" "<=" ">="}
	returns (int) : location of first
*/
int find_last( float* buf, int n, char* type, float val );

/* returns deep slice of buffer 
	inputs
		(float*) buf : buffer you wish to slice
		(int) n : length of buf
		(int) start : starting index to slice
		(int) end : ending index to slice
		(int*) new_n : new length of slice buffer pass by reference.
	returns
		(float*) sliced_buf : new buffer of length new_n
*/
float* slice( float* buf, int n, int start, int end, int *new_n );

/* NOT IMPLEMENTED YET.  CAN SOLVE WITH FIND_FIRST
	
	trims array where start > tol and end < tol
	inputs
		(float *) buf : buffer you wish to trim
		(int) n : length of buf
		(float) tol : tol value to start and stop
		(int *) new_n : new length of trimmed array
	returns
		(float *) : pointer to buffer of trimmed array

float* trim_zeros( float *buf, int n, float tol, int *new_n );
*/
/*
	inputs
		(float*) buf : buffer to integrate
		(int) n : length of buf
		(float) dh : array step
	returns
		(float) val : integrated value
*/
float trapz( float* buf, int n, float dh );

/*
	inputs
		(float *) buf : buffer 
		(int) n : length of buf
		(float) dh : grid step
	returns
		(float*) int_buf : integrated buffer of length n
*/
float* cumtrapz( float *buf, int n, float dh );

/*
	inputs
		(float *) vector : vector
	return
		(float) norm : norm
*/
float norm( float *vector, int ndim );

/*
	inputs
		(double *) vector : vector
	return
		(double) norm : norm
*/
double norm_d( double *vector, int ndim );


/* modifies buf in-place so that norm(buf, 3) = 1.0
	inputs
		(float *) buf : buffer
		(int) ndim : length of buffer
	return
		NONE
*/
void to_unit_vector( float *buf, int ndim );

/* modifies buf in-place so that norm(buf, 3) = 1.0
	inputs
		(float *) buf : buffer
		(int) ndim : length of buffer
	return
		NONE
*/
float* unit_vector_3d( float n1, float n2, float n3 );

