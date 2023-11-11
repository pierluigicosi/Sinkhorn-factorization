CC=mpicc
CFLAGS=-g

all: sinkhorn compiler

sinkhorn: sinkhorn.o MyMPI.o
	$(CC) $(CFLAGS) $^ -o $@

sinkhorn.o: sinkhorn.c MyMPI.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o sinkhorn