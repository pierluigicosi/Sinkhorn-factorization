CC=mpicc
CFLAGS=-g

all: sinkhorn

sinkhorn: sinkhorn.o ../util/MyMPI.o ../util/generator.o
	$(CC) $(CFLAGS) $^ -O3 -o $@

sinkhorn.o: sinkhorn.c ../util/MyMPI.h ../util/generator.h
	$(CC) $(CFLAGS) -c $<  -O3 -o $@

clean:
	rm *.o sinkhorn