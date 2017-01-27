CC = gcc
CFLAGS = -std=c99 -O3 -Wall
LIBRARY_FLAGS = -lgsl -lgslcblas -lm

main:
	$(CC) $(CFLAGS) -c mode_polynomial.c
	$(CC) $(CFLAGS) $(LIBRARY_FLAGS) -o mode_polynomial mode_polynomial.o

clean:
	rm -f mode_polynomial mode_polynomial.o
