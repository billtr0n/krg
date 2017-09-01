#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

float* linspace(float start, float end, int length) {
    float* buf;
    float dx;
    int i;
    
    dx = (end - start) / (length-1);
    buf = (float*)malloc(length*sizeof(float));
    for (i=0; i<length; i++) { 
        buf[i] = start + dx*i;
    }
    return buf;
}

float* arange(float start, float end, float dx, int *n) {
    float* buf;
    int length;
    int i;
    
    length = ceil( (end - start) / dx);
    buf = malloc(length*sizeof(float));
    for (i=0; i<length; i++) {
        buf[i] = start + dx*i;
    }
    *n = length;
    return buf;

}

float* zeros_f(int n) {
    float* buf;
    int i;

    buf = (float*)malloc(n*sizeof(float));
    for (i=0; i<n; i++) {
        buf[i] = 0.0;
    }

    return buf;
}

int* zeros_i(int n) {
    int* buf;
    int i;

    buf = (int*)malloc(n*sizeof(int));
    for (i=0; i<n; i++) {
        buf[i] = 0;
    }

    return buf;
}

float sum(float *buf, int n) {
    float total;
    int i;
    total = 0;
    for (i=0;i<n;i++) {
        total += buf[i];
    }

    return total;
}

float** cartesian_product( float **A, int *np, int *ncom, int ndim ) {
    int i, j;
    int n = 1;
    bool DEBUG = false;
    int *c;
    float **B;

    // compute ncom
    for (i=0; i<ndim; i++) {
        n *= np[i];
    }
    *ncom = n;

    // allocate count array
    c = zeros_i( ndim );

    // allocate output array
    B = (float**) malloc( sizeof(float*) * n );
    for (i=0; i<n; i++) {
        B[i] = (float*) malloc(sizeof(float)*ndim);
    }

    // prepare combinations
    if (DEBUG) {
        printf("\n");
        printf("Combinations\n");
        printf("============\n");
    }

    for (i=0; i<n; i++) {
        for (j=0; j<ndim; j++) {
            if (c[j]==np[j]) {
                c[j+1]++;
                c[j]=0;
            }
            B[i][j] = A[j][c[j]];
        }
        if (DEBUG) {
            printf("{%d, %d, %d} -> {%f, %f %f}\n", c[0], c[1], c[2], B[i][0], B[i][1], B[i][2]);
        }
        
        c[0]++;
    }
        
    // all done!
    return B;

}


bool compare_f( float a, float b, float eps ) {
    return fabs(a-b) < eps;       
}

void dealloc2d_f( float** buf, int n ) {
    int i;
    for (i=0; i<n; i++) {
        free(buf[i]);
    }
    free(buf);
}

void dealloc2d_i( int** buf, int n ) {
    int i;

    for (i=0; i<n; i++) {
        free(buf[i]);
    }
    free(buf);
}

float **outer( float *v1, float *v2, int n ) {
    int i, j;
    float **m;

    m = (float**) malloc(sizeof(float*) * n);
    for (i=0; i<n; i++) {
        m[i] = (float*) malloc(sizeof(float) * n);
    }

    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            m[i][j] = v1[i] * v2[j];
        }
    }

    return m;

}

float rastrigin( float x, float y ) {
    float f;
    float p1,p2,p3;
    int i, j;
    int d = 2;
    float pi = 3.14159265359;

    p1 = 10*d;
    p2 = x*x + y*y;
    p3 = 10 * (cos(2*pi*x) + cos(2*pi*y));
    
    f = p1 + p2 - p3;
    return f;
}

float ackley( float x, float y ) {
    float p1, p2, f;
    float pi = 3.14159265359;

    p1 = 20*exp(-0.2 * sqrt(0.5*(x*x+y*y)));
    p2 = exp(0.5*(cos(2*pi*x) + cos(2*pi*y)));
    f = -p1 - p2 + exp(1) + 20;
    return f;
}

int find_first( float* buf, int n, char* type, float val) {
    int i;

    for (i=0; i<n; i++) {
        if (strcmp(type, ">") == 0) {
            if (buf[i] > val) return i;
        }
        else if (strcmp(type, "<") == 0) {
            if (buf[i] < val) return i;
        }
        else if (strcmp(type, ">=") == 0) {
            if (buf[i] >= val) return i;
        }
        else if (strcmp(type, "==") == 0) {
            if (buf[i] == val) return i;
        }
        else if (strcmp(type, "!=") == 0) {
            if (buf[i] != val) return i;

        } else {
            printf("Error: Invalid type. Exiting program.");
            exit(-1);
        }
    }
    return -1;
}

int find_last( float* buf, int n, char* type, float val ) {
    int i;

    // to find last, just find first but looping through from end.
    for (i=n-1; i>=0; i--) {
        if (strcmp(type, ">") == 0) {
            if (buf[i] > val) return i;
        }
        else if (strcmp(type, "<") == 0) {
            if (buf[i] < val) return i;
        }
        else if (strcmp(type, ">=") == 0) {
            if (buf[i] >= val) return i;
        }
        else if (strcmp(type, "==") == 0) {
            if (buf[i] == val) return i;
        }
        else if (strcmp(type, "!=") == 0) {
            if (buf[i] != val) return i;

        } else {
            printf("Error: Invalid type. Exiting program.");
            exit(-1);
        }
    }
    return -1;
}

float* slice( float* buf, int n, int start, int end, int *new_n ) {
    int i;
    float* new;
    *new_n = n - (start + ((n-1) - end));
    printf("start=%d end=%d new_n=%d\n", start, end, *new_n);
    new = malloc(sizeof(float) * *new_n);
    for (i=0; i<*new_n; i++) {
        new[i] = buf[start+i];
    }
    return new;
}

float maximum( float* buf, int n ) {
    float left, right;
    float max;

    if (n == 1) {
        return buf[0];
    } else {
        left = maximum(buf, n/2);
        right = maximum(&buf[n/2], n-n/2);
    } 
    return (left > right) ? left : right;
}

float trapz( float* buf, int n, float dh ) {
    int i;
    float s = 0;

    for (i=1; i<n; i++) {
        s = s + (buf[i-1] + buf[i]);
    }
    s = 0.5 * dh * s;

    return s;
}

float* cumtrapz( float* buf, int n, float dh ) {
    float *s;
    int i;

    s = zeros_f(n);

    for (i=1; i<n; i++) {
        s[i] = s[i-1] + 0.5*dh*(buf[i-1] + buf[i]); 
    }

    return s;
}

float norm( float *vector, int ndim ) {
    int i;
    float s = 0;

    for (i=0; i<ndim; i++) {
        s += vector[i]*vector[i];
    }
    return sqrt(s);
}

float norm_d( double *vector, int ndim ) {
    int i;
    double s = 0;

    for (i=0; i<ndim; i++) {
        s += vector[i]*vector[i];
    }
    return s;
}

void to_unit_vector( float *buf, int ndim ) {
    int i;
    float scaling;

    scaling = 1.0 / norm(buf, 3);
    for (i=0; i<ndim; i++) {
        buf[i] *= scaling;
    }
    return;
}

float* unit_vector_3d( float n1, float n2, float n3 ) {
    int i;
    int ndim = 3;
    float s;
    float *out;

    out = malloc(sizeof(float)*ndim);
    // create array
    out[0] = n1;
    out[1] = n2;
    out[2] = n3;
    s = 1.0 / norm(out, ndim);

    // scale
    for( i=0; i<ndim; i++ ) {
        out[i] *= s;
    }
    return out;
}
/*

int main() {
    float *buf;
    float sum1;
    int i;
    buf = (float*)malloc(10*sizeof(float));

    for (i=0; i<10; i++) {
        buf[i] = i*1e8;
    }
    sum1 = sum(buf, 10);

    fprintf(stdout, "%f\n", sum1);

    return 0;
}
    
int main() {
    float start = 0;
    float end = 10;
    float dx = 1;
    int length = 11;
    float *buf1, *buf2;
    int i;
    int temp;
    
    buf1 = linspace(start, end, length);
    buf2 = arange(start, end+0.1, dx, &temp);

    
    for (i=0; i<length; i++) {
        printf("%f %f\n", buf1[i], buf2[i]);
    } 
    
    free(buf1);
    free(buf2);
    return 0;
}
*/

/*
// testing first_val
int main() {
    float *sv1;
    int i;
    int arg1,arg2;
    FILE *fin;
    int nt = 10001;

    // read test slip-rate function
    sv1 = (float*) malloc(sizeof(float) * nt);
    fin = fopen( "./testing/scratch/sv1_1300_400.bin", "rb" );
    fread( sv1, sizeof(float), nt, fin );
    fclose( fin );

    arg1 = find_first(sv1, nt, ">", 0.001);
    printf("start=%d\n", arg1);
    arg2 = find_first(&sv1[arg1+1], nt-arg1-1, "<", 0.001);
    printf("end=%d\n", arg1+arg2);
    return 0;
}*/

/* testing slice
int main() {
    float *buf;
    int i;
    int n = 10;
    int new_n;
    float* arg;

    buf = malloc(sizeof(float)*n);

    for (i=0; i<n; i++) {
        buf[i] = i*0.1;
    }
    arg = slice(buf, n, 2, 7, &new_n);
    printf("%d\n", new_n);
    return 0;
}
 */
/* testing recursive maximum 
int main() {
    float *buf;
    int i;
    int n = 10001;
    int new_n;
    float pi = 3.14159;
    float* arg;

    buf = malloc(sizeof(float)*n);

    for (i=0; i<n; i++) {
        buf[i] = cos(2*pi*i*0.1);
    }
    float max = maximum(buf, n);
    printf("max = %f\n", max);
    return 0;
}*/

/*
// testing to_unit_vector
int main() {
    float *buf;
    int i;
    int n = 100;
    int new_n;
    float pi = 3.14159;
    float arg[3];

    arg[0] = 1.0;
    arg[1] = 0.0;
    arg[2] = 1.0;

    to_unit_vector(arg, 3);

    printf("arg[0]=%f arg[1]=%f arg[2]=%f\n", arg[0], arg[1], arg[2]);
    
    return 0;
}
*/


/* testing unit_vector_3d
int main() {
    float *buf;
    int i;
    int n = 100;
    int new_n;
    float pi = 3.14159;
    float arg[3];

    arg[0] = 1.0;
    arg[1] = 0.0;
    arg[2] = 1.0;

    buf = unit_vector_3d(arg[0], arg[1], arg[2]);

    printf("buf[0]=%f buf[1]=%f buf[2]=%f\n", buf[0], buf[1], buf[2]);
    
    return 0;
} */
