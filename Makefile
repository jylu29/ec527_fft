LIBS = -lm -lrt
CC = gcc
all: clean fft_serial fft_omp

fft_serial: fft_serial.c
	$(CC) -O1 fft_serial.c $(LIBS) -o fft_serial

fft_serial_precompute: fft_serial_precompute.c
	$(CC) -O1 fft_serial_precompute.c $(LIBS) -o fft_serial_precompute

fft_omp: fft_omp.c
	$(CC) -O1 -fopenmp fft_omp.c $(LIBS) -o fft_omp

fft_omp_precompute: fft_omp_precompute.c
	$(CC) -O1 -fopenmp fft_omp_precompute.c $(LIBS) -o fft_omp_precompute

fft_omp_precompute_transpose: fft_omp_precompute_transpose.c
	$(CC) -O0 -fopenmp -g fft_omp_precompute_transpose.c $(LIBS) -o fft_omp_precompute_transpose

clean:
	rm -rf fft_serial fft_omp fft_omp_precompute_transpose
