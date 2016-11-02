/* creates 1d array of values starting with "start" ending with "end" and a
   length of "length". don't forget to free buffer in main when done! */
float* linspace(float start, float end, int length);

/* creates 1d array of values starting with "start" ending with "end" and a
   spacing of dx. don't forget to free buffer in main when done! */
float* arange(float start, float end, float dx, int *n);

/* computes sum of 1d array */
float sum(float *buf, int n);

/* returns value of hann window with length, n, at step i */
float hann(int i, int n);
