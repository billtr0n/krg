#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"

#define PI 3.14159265358979323846

/* helping functions to build the tinti function, r stands for range. */
float r1(float t, float Tr) {
	float p1;

	p1 = (PI*(4*t - Tr)*Tr + 8*pow(t,1.5)*sqrt(-t + Tr) + 4*Tr*sqrt(t*(-t + Tr)) - 
     2*(4*t - Tr)*Tr*atan((-2*t + Tr)/(2.*sqrt(t*(-t + Tr)))))/16.;

	return p1;
}

float r2(float t, float Tr, float Ts) {
	float p1;
	float p2;

	p1 = (Tr*sqrt(t*(-t + Tr)) + 2*sqrt(pow(t,3)*(-t + Tr)) - t*Tr*sqrt((-t + Tr + Ts)/(t - Ts)) + 
      Tr*Ts*sqrt((-t + Tr + Ts)/(t - Ts)) - 2*t*sqrt((t - Ts)*(-t + Tr + Ts)) - 
      2*Ts*sqrt((t - Ts)*(-t + Tr + Ts)) + (4*t - Tr)*Tr*asin(sqrt(t/Tr)) + 
      Tr*(-4*t + Tr)*asin(sqrt((t - Ts)/Tr)))/4.;

	p2 = (-((2*t + Tr - 6*Ts)*sqrt((t - Ts)*(-t + Tr + Ts))) + 
      Tr*(-4*t + Tr + 8*Ts)*asin(sqrt((t - Ts)/Tr)))/4.;

	return p1+p2;
}

float r3(float t, float Tr, float Ts) {
	float p1;
	float p2;

	p1 = (Tr*sqrt(t*(-t + Tr)) + 2*sqrt(pow(t,3)*(-t + Tr)) - t*Tr*sqrt((-t + Tr + Ts)/(t - Ts)) + 
      Tr*Ts*sqrt((-t + Tr + Ts)/(t - Ts)) - 2*t*sqrt((t - Ts)*(-t + Tr + Ts)) - 
      2*Ts*sqrt((t - Ts)*(-t + Tr + Ts)) + (4*t - Tr)*Tr*asin(sqrt(t/Tr)) + 
      Tr*(-4*t + Tr)*asin(sqrt((t - Ts)/Tr)))/4.;

	p2 = (2*(-2*t*sqrt((t - Ts)*(-t + Tr + Ts)) - Tr*sqrt((t - Ts)*(-t + Tr + Ts)) + 
      6*Ts*sqrt((t - Ts)*(-t + Tr + Ts)) + 2*t*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)) + 
      Tr*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)) - 4*Ts*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts))) - 
      Tr*(-4*t + Tr + 8*Ts)*atan((-2*t + Tr + 2*Ts)/(2.*sqrt((t - Ts)*(-t + Tr + Ts)))) + 
      Tr*(-4*t + Tr + 8*Ts)*atan((-2*t + Tr + 4*Ts)/(2.*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)))))/8.;

	return p1+p2;
}

float r4(float t, float Tr, float Ts) {
	float p1;
	float p2;
	
	p1 = t*Tr*acos(sqrt((t - Ts)/Tr)) + (-(sqrt((t - Ts)*(-t + Tr + Ts))*(2*t + Tr + 2*Ts)) - 
      pow(Tr,2)*acos(1.0/sqrt(Tr/(t - Ts))))/4.;

    p2 = (2*(-2*t*sqrt((t - Ts)*(-t + Tr + Ts)) - Tr*sqrt((t - Ts)*(-t + Tr + Ts)) + 
         6*Ts*sqrt((t - Ts)*(-t + Tr + Ts)) + 2*t*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)) + 
         Tr*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)) - 4*Ts*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts))) - 
      Tr*(-4*t + Tr + 8*Ts)*atan((-2*t + Tr + 2*Ts)/(2.*sqrt((t - Ts)*(-t + Tr + Ts)))) + 
      Tr*(-4*t + Tr + 8*Ts)*atan((-2*t + Tr + 4*Ts)/(2.*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)))))/8.;

    return p1+p2;
}

float r5(float t, float Tr, float Ts) {
	float p1;

	p1 = (4*(2*t + Tr - 4*Ts)*sqrt(-((t - 2*Ts)*(t - Tr - 2*Ts))) + PI*Tr*(-4*t + Tr + 8*Ts) + 
     	  2*Tr*(-4*t + Tr + 8*Ts)*atan((-2*t + Tr + 4*Ts)/(2.*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)))))/16.;

   return p1;
}

/* used for different integration range when Ts < Tr < 2 Ts */
float sr3(float t, float Tr, float Ts) {
	float p1, p2;

	p1 = t*Tr*acos(sqrt((t - Ts)/Tr)) + (-(sqrt((t - Ts)*(-t + Tr + Ts))*(2*t + Tr + 2*Ts)) - 
      pow(Tr,2)*acos(1.0/sqrt(Tr/(t - Ts))))/4.;

	p2 = (-4*(2*t + Tr - 6*Ts)*sqrt((t - Ts)*(-t + Tr + Ts)) + PI*Tr*(-4*t + Tr + 8*Ts) - 
      2*Tr*(-4*t + Tr + 8*Ts)*atan((-2*t + Tr + 2*Ts)/(2.*sqrt((t - Ts)*(-t + Tr + Ts)))))/16.;

	return p1+p2;
}

float* tinti(float* t, int nt, float Tr, float Ts, float t0, float slip) {
	int i;
	float t_shift;
	float k;
	float *stf_out;
	
	// allocate zeroed buffer
	stf_out = zeros_f(nt);

	// constant during integration
	k = 2 / (PI*Tr*pow(Ts,2));

	// printf("(in tinti): nt=%d Tr=%f Ts=%f t0=%f slip=%f\n", nt, Tr, Ts, t0, slip);

	// generate source time function
	for (i=0; i<nt; i++) {
		
		// allow shift in start of source-time function
		if (t0 > 0.0) {
			t_shift = t[i] - t0;
		}
		else {
			t_shift = t[i];
		}
		

		// condition 1 of integral: Tr > 2*Ts
		if (Tr > 2*Ts) {
			// 0 < t < ts
			if ((t_shift >= 0) && (t_shift <= Ts)) {
				stf_out[i] = k * r1(t_shift, Tr);
			}
			
			// ts < t < 2*ts
			else if ((t_shift > Ts) && (t_shift < 2*Ts)) {
				stf_out[i] = k * r2(t_shift, Tr, Ts);
			}

			// 2*ts < t < tr
			else if ((t_shift >= 2*Ts) && (t_shift < Tr)) {
				stf_out[i] = k * r3(t_shift, Tr, Ts);
			}

			// tr < t < tr + ts
			else if ((t_shift >= Tr)&&(t_shift < Tr+Ts)) {
				stf_out[i] = k * r4(t_shift, Tr, Ts);
			}

			// ts + ts < t < tr + 2*ts
			else if ((t_shift >= Tr+Ts)&&(t_shift < Tr+2*Ts)) {
				stf_out[i] = k * r5(t_shift, Tr, Ts);
			}

			else {
				stf_out[i] = 0.0;
			}
		}
		// condition 2 of integral: Ts < Tr < 2*Ts
		else if ((Tr > Ts) && (Tr <= 2*Ts)) {
			// 0 < t < ts
			if ((t_shift>=0)&&(t_shift<=Ts)) {
				stf_out[i] = k * r1(t_shift, Tr);
				// printf("(r1) %f\n", stf_out[i]);
			}

			// ts < t < tr
			else if ((t_shift > Ts) && (t_shift < Tr)) {
				stf_out[i] = k * r2(t_shift, Tr, Ts);
				// printf("(r2) %f\n", stf_out[i]);
			}

			// tr < t < 2*ts
			else if ((t_shift >= Tr) && (t_shift < 2*Ts)) {
				stf_out[i] = k * sr3(t_shift, Tr, Ts);
				// printf("(sr3) %f\n", stf_out[i]);
			}

			// 2*ts < t < ts + tr
			else if ((t_shift >= 2*Ts) && (t_shift < Ts+Tr)) {
				stf_out[i] = k * r4(t_shift, Tr, Ts);
				// printf("(r4) %f\n", stf_out[i]);
			}

			// ts + tr < t < tr + 2*ts
			else if ((t_shift >= Ts+Tr) && (t_shift<2*Ts+Tr)) {
				stf_out[i] = k * r5(t_shift, Tr, Ts);
				// printf("(r5) %f\n", stf_out[i]);
			}	

			else {
				stf_out[i] = 0.0;
				// printf("(else) %f\n", stf_out[i]);
			}
		}

		// sanity check
		else {
			fprintf(stderr, "Warning! Tr (%f) must be greater than Ts (%f)!\n", Tr, Ts);
		}

		// scale normalized slip-rate by slip.
		stf_out[i] *= slip;
	}

	
	return stf_out;
}


float tinti_scalar(float t, int nt, float Tr, float Ts, float t0, float slip) {
	int i;
	float t_shift;
	float k;
	float stf_out = 0;
  	
	// constant during integration
	k = 2 / (PI*Tr*pow(Ts,2));

	// printf("(in tinti): nt=%d Tr=%f Ts=%f t0=%f slip=%f\n", nt, Tr, Ts, t0, slip);

	// generate source time function
	for (i=0; i<nt; i++) {
		
		// allow shift in start of source-time function
		if (t0 > 0.0) {
			t_shift = t - t0;
		}
		else {
			t_shift = t;
		}
		

		// condition 1 of integral: Tr > 2*Ts
		if (Tr > 2*Ts) {
			// 0 < t < ts
			if ((t_shift >= 0) && (t_shift <= Ts)) {
				stf_out = k * r1(t_shift, Tr);
			}
			
			// ts < t < 2*ts
			else if ((t_shift > Ts) && (t_shift < 2*Ts)) {
				stf_out = k * r2(t_shift, Tr, Ts);
			}

			// 2*ts < t < tr
			else if ((t_shift >= 2*Ts) && (t_shift < Tr)) {
				stf_out = k * r3(t_shift, Tr, Ts);
			}

			// tr < t < tr + ts
			else if ((t_shift >= Tr)&&(t_shift < Tr+Ts)) {
				stf_out = k * r4(t_shift, Tr, Ts);
			}

			// ts + ts < t < tr + 2*ts
			else if ((t_shift >= Tr+Ts)&&(t_shift < Tr+2*Ts)) {
				stf_out = k * r5(t_shift, Tr, Ts);
			}

			else {
				stf_out = 0.0;
			}
		}
		// condition 2 of integral: Ts < Tr < 2*Ts
		else if ((Tr > Ts) && (Tr <= 2*Ts)) {
			// 0 < t < ts
			if ((t_shift>=0)&&(t_shift<=Ts)) {
				stf_out = k * r1(t_shift, Tr);
				// printf("(r1) %f\n", stf_out);
			}

			// ts < t < tr
			else if ((t_shift > Ts) && (t_shift < Tr)) {
				stf_out = k * r2(t_shift, Tr, Ts);
				// printf("(r2) %f\n", stf_out);
			}

			// tr < t < 2*ts
			else if ((t_shift >= Tr) && (t_shift < 2*Ts)) {
				stf_out = k * sr3(t_shift, Tr, Ts);
				// printf("(sr3) %f\n", stf_out);
			}

			// 2*ts < t < ts + tr
			else if ((t_shift >= 2*Ts) && (t_shift < Ts+Tr)) {
				stf_out = k * r4(t_shift, Tr, Ts);
				// printf("(r4) %f\n", stf_out);
			}

			// ts + tr < t < tr + 2*ts
			else if ((t_shift >= Ts+Tr) && (t_shift<2*Ts+Tr)) {
				stf_out = k * r5(t_shift, Tr, Ts);
				// printf("(r5) %f\n", stf_out);
			}	

			else {
				stf_out = 0.0;
				// printf("(else) %f\n", stf_out);
			}
		}

		// sanity check
		else {
			fprintf(stderr, "Warning! Tr (%f) must be greater than Ts (%f)!\n", Tr, Ts);
		}

		// scale normalized slip-rate by slip.
		stf_out *= slip;
	}

	
	return stf_out;
}

/*
int main() {
	float Ts = 0.1;
	float Tr = 2.0;
	float t0 = 1.0;
	float slip = 1.0;
	int i;
	float *t; 
	float *stf;
	int nt;
	float test;

	// test limiting conditions on function
	// define time vector
	t = arange(0, 10, 0.001, &nt);

	// generate tinti source-time function
	stf = tinti(t, nt, Tr, Ts, t0, slip);

	for (i=0; i<nt; i++) {
    		printf("%f\n", stf[i]);
	}

	 
	//   testing values at function boundaries.  
	//  note: values of inf are OK so long as every interval has an inequality 
	//   that provides finite values. 
	// test boundary conditions for all the integral ranges
	printf("Case 1: Tr > 2*Ts\n");
	printf("r1: %f\n", r1(0, Tr));
	printf("r1: %f\n", r1(Ts, Tr));
	printf("r2: %f\n", r2(Ts, Tr, Ts));
	printf("r2: %f\n", r2(2 * Ts, Tr, Ts));
	printf("r3: %f\n", r3(2 * Ts, Tr, Ts));
	printf("r3: %f\n", r3(Tr, Tr, Ts));
	printf("r4: %f\n", r4(Tr, Tr, Ts));
	printf("r4: %f\n", r4(Tr + Ts, Tr, Ts));
	printf("r5: %f\n", r5(Tr + Ts, Tr, Ts));
	printf("r5: %f\n", r5(Tr + 2*Ts, Tr, Ts));

	// second range only differs in the "third" integration range
	Ts = 1.1;
	Tr = 2.0;
	// generate tinti source-time function
	stf = tinti(t, nt, Tr, Ts, t0, slip);
	printf("Case 2: Ts < Tr < 2*Ts\n");
	printf("r1: %f\n", r1(0, Tr));
	printf("r1: %f\n", r1(Ts, Tr));
	printf("r2: %f\n", r2(Tr, Tr, Ts));
	printf("r2: %f\n", r2(Tr, Tr, Ts));
	printf("sr3: %f\n", sr3(Tr, Tr, Ts));
	printf("sr3: %f\n", sr3(2 * Ts, Tr, Ts));
	printf("r4: %f\n", r4(2 * Ts, Tr, Ts));
	printf("r4: %f\n", r4(Tr + Ts, Tr, Ts));
	printf("r5: %f\n", r5(Tr + Ts, Tr, Ts));
	printf("r5: %f\n", r5(Tr + 2*Ts, Tr, Ts));

	free(stf);
	free(t);

	return 0;
    }
    */
