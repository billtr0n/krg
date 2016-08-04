#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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
    buf = (float*)malloc(length*sizeof(float));
    for (i=0; i<length; i++) {
        buf[i] = start + dx*i;
    }
    *n = length;
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
