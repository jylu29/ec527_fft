// gcc -O1 -fopenmp fft_omp.c -lrt -lm -o fft_omp
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include <string.h>
#include <omp.h>

#define MINVAL   0.0
#define MAXVAL  1000.0
#define N 4194304
#define THREADS 16

double interval(struct timespec start, struct timespec end);
/**
 * perform fft using Cooley–Tukey algorithm
 */
void fft(const float complex *input, float complex *output, int n, double *time_array) {
    memcpy(output, input, sizeof(*input)*n);
    // Bit-reverse the input array
    for (int i = 0, j = 0; i < n; i++) {
        if (j > i) {
            float complex tmp = output[j];
            output[j] = output[i];
            output[i] = tmp;
        }
        int m = n / 2;
        while (m >= 2 && j >= m) {
            j -= m;
            m /= 2;
        }
        j += m;
    }
    struct timespec start, stop;
    // Compute the FFT using butterfly operations
    for (int s = 1; s <= log2(n); s++) {
        int m = 1 << s;
        clock_gettime(CLOCK_REALTIME, &start);
        #pragma omp parallel for
        for (int k = 0; k < n; k += m) {
            for (int j = 0; j < m / 2; j++) {
                float complex t = cexp(-2*j*M_PI*I/m);
                float complex u = output[k + j];
                t = t * output[k + j + m / 2];
                output[k + j] = u + t;
                output[k + j + m / 2] = u - t;
            }
        }
        clock_gettime(CLOCK_REALTIME, &stop);
        time_array[s-1] = interval(start, stop);
    }
}


double fRand(double fMin, double fMax) {
    double f = (double) random() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void initialize(float complex *ptr) {
    srandom(N);
    for (int i = 0; i < N; i++) {
        float complex c;
        c = (float)fRand((double) (MINVAL), (double) (MAXVAL))+(float)fRand((double) (MINVAL), (double) (MAXVAL))*I;
//        c = i+1;
        ptr[i] = c;
    }
}

double interval(struct timespec start, struct timespec end) {
    struct timespec temp;
    temp.tv_sec = end.tv_sec - start.tv_sec;
    temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    if (temp.tv_nsec < 0) {
        temp.tv_sec = temp.tv_sec - 1;
        temp.tv_nsec = temp.tv_nsec + 1000000000;
    }
    return (((double) temp.tv_sec) + ((double) temp.tv_nsec) * 1.0e-9);
}

int main() {
    float complex *x = malloc(sizeof(float complex) * N);
    float complex *result = malloc(sizeof(float complex) * N);
    // define input signal
    initialize(x);

    omp_set_num_threads(THREADS);

    struct timespec time_start, time_stop;
    double time_array[(int)log2(N)];
    // compute FFT
    clock_gettime(CLOCK_REALTIME, &time_start);
    for (int i = 0; i < 10; i++) {
        fft(x, result, N, time_array);
    }
    clock_gettime(CLOCK_REALTIME, &time_stop);
    double time = interval(time_start, time_stop);


    // print result
    printf("FFT: \n");
    for (int i = 0; i < 8; i++) {
        printf("%f + %fi\n", creal(result[i]), cimag(result[i]));
    }
    printf("\n");
    printf("Time = %8.4g s\n", time/10);

    printf("Stage Time\n");
    for (int i = 0; i < log2(N); ++i) {
        printf("%d %f\n",i, time_array[i]);
    }

    return 0;
}

