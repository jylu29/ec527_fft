// gcc -O1 fft_serial.c -lrt -lm -o fft_serial
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include <string.h>

#define MINVAL   0.0
#define MAXVAL  1000.0
#define N 4096

/**
 * perform fft using Cooley–Tukey algorithm
 */
void fft(const double complex input[], double complex output[], int n) {
    memcpy(output, input, sizeof(*input)*n);
    // Bit-reverse the input array
    for (int i = 0, j = 0; i < n; i++) {
        if (j > i) {
            double complex tmp = output[j];
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
    // Compute the FFT using butterfly operations
    for (int s = 1; s <= log2(n); s++) {
        int m = (int) pow(2, s);
        double complex wm = cexp(-2 * M_PI * I / m);
        for (int k = 0; k < n; k += m) {
            double complex w = 1.0;
            for (int j = 0; j < m / 2; j++) {
                double complex t = w * output[k + j + m / 2];
                double complex u = output[k + j];
                output[k + j] = u + t;
                output[k + j + m / 2] = u - t;
                w *= wm;
            }
        }
    }
}


double fRand(double fMin, double fMax) {
    double f = (double) random() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void initialize(double complex *ptr) {
    srandom(N);
    for (int i = 0; i < N; i++) {
        double complex c;
        c = fRand((double) (MINVAL), (double) (MAXVAL))+fRand((double) (MINVAL), (double) (MAXVAL))*I;
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
    double complex x[N];
    double complex result[N];
    // define input signal
    initialize(x);

    struct timespec time_start, time_stop;
    // compute FFT
    clock_gettime(CLOCK_REALTIME, &time_start);
    for (int i = 0; i < 10; i++) {
        fft(x, result, N);
    }
    clock_gettime(CLOCK_REALTIME, &time_stop);
    double time = interval(time_start, time_stop);


    // print result
    printf("FFT: \n");
    for (int i = 0; i < 8; i++) {
        printf("%f + %fi\n", creal(result[i]), cimag(result[i]));
    }
    printf("\n");
    printf("Time = %8.4g s\n", time);

    return 0;
}

