LIBS = -lm -lrt
CC = gcc
all: clean fft_serial fft_omp

fft_serial: fft_serial.c
	$(CC) -O1 -g fft_serial.c $(LIBS) -o fft_serial

fft_omp: fft_omp.c
	$(CC) -O1 -fopenmp -g fft_omp.c $(LIBS) -o fft_omp

clean:
	rm -rf fft_serial fft_omp
