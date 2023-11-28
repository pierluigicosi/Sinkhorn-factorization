CC=mpicc
CFLAGS=-g

all: sinkhorn

sinkhorn: sinkhorn.o ../util/MyMPI.o ../util/generator.o
	$(CC) $(CFLAGS) $^ -o $@

sinkhorn.o: sinkhorn.c ../util/MyMPI.h ../util/generator.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o sinkhorn