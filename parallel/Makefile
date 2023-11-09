CC=mpicc
CFLAGS=-g

all: floyd compiler

floyd: sinkhorn.o MyMPI.o
	$(CC) $(CFLAGS) $^ -o $@

floyd.o: sinkhorn.c MyMPI.h
	$(CC) $(CFLAGS) -c $< -o $@

compiler: compiler.o
	$(CC) $(CFLAGS) $^ -o $@

compiler.o: compiler.c 
	$(CC) $(CFLAGS) -c $< -o $@

MyMPI.o: MyMPI.c MyMPI.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o floyd compiler