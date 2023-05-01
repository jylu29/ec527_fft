// nvcc -arch compute_35 fft_cuda.cu -o fft_cuda
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <cuComplex.h>
#include "cuPrintf.cu"

#define MINVAL   0.0
#define MAXVAL  1000.0
#define N 4194304
#define BLOCK_SIZE 1024

/**
 * perform fft_kernel using Cooley–Tukey algorithm
 */
__global__ void fft_kernel(cuComplex *result, int m) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    // Compute the FFT using butterfly operations
    if ((j&m)==0) {
        cuComplex t;
        float tmp = -1.0 * M_PI * j / (1.0 * m);
        t.x = cosf(tmp);
        t.y = sinf(tmp);
        t = cuCmulf(result[j + m], t);
        cuComplex u = result[j];
//    cuPrintf("t = %f + %fi, m=%d\n", t.x, t.y, m);

        result[j] = cuCaddf(u, t);
        result[j + m] = cuCsubf(u, t);
    }
}


double fRand(double fMin, double fMax) {
    double f = (double) random() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void initialize(cuComplex *ptr) {
    srandom(N);
    for (int i = 0; i < N; i++) {
        cuComplex c;
        c.y = (float) fRand((double) (MINVAL), (double) (MAXVAL));
        c.x = (float) fRand((double) (MINVAL), (double) (MAXVAL));
//        c.x = i+1;
//        c.y=0;
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

//bit-reverse order of d, store in r
__global__ void bit_rev_reorder(cuComplex *__restrict__ r, cuComplex *__restrict__ d, int s) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    r[__brev(id) >> (32 - s)] = d[id];
//    cuPrintf("id %d", id);
}

int main() {
    auto *h_x = (cuComplex *) malloc(sizeof(cuComplex) * N);
    auto *h_result = (cuComplex *) malloc(sizeof(cuComplex) * N);
    cuComplex *d_x, *d_result;
    cudaMalloc((void **) &d_x, N * sizeof(cuComplex));
    cudaMalloc((void **) &d_result, N * sizeof(cuComplex));
    // define input signal
    initialize(h_x);

    struct timespec time_start, time_stop;
    // compute FFT
    clock_gettime(CLOCK_REALTIME, &time_start);
    cudaMemcpy(d_x, h_x, N * sizeof(cuComplex), cudaMemcpyHostToDevice);
    cudaSetDevice(0);
    cudaEvent_t start;
    cudaEventCreate(&start);
    cudaEventRecord(start, 0);
    cudaPrintfInit();
    for (int i = 0; i < 10; i++) {
        int stages = (int) log2(N);
        // bit_reverse with cuda
        bit_rev_reorder<<<N / BLOCK_SIZE, BLOCK_SIZE>>>(d_result, d_x, stages);
        cudaDeviceSynchronize();
        for (int s = 0; s < stages; s++) {
            int m = 1 << s;
//            for (int k = 0; k < N; k += m) {
            fft_kernel<<<(int) (N + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE>>>(d_result, m);
                cudaDeviceSynchronize();
//            }
        }
    }
    cudaPrintfDisplay(stdout, true);
    cudaPrintfEnd();
    cudaMemcpy(h_result, d_result, N * sizeof(cuComplex), cudaMemcpyDeviceToHost);
    clock_gettime(CLOCK_REALTIME, &time_stop);
    double time = interval(time_start, time_stop);


    // print result
    printf("FFT: \n");
    for (int i = 0; i < 8; i++) {
        printf("%f + %fi\n", h_result[i].x, h_result[i].y);
    }
    printf("\n");
    printf("Time = %8.4g s\n", time);

    return 0;
}

