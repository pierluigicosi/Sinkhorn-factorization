CC = gcc
CFLAGS = -g

sinkhorn: sinkhorn.o ../util/generator.o
	$(CC) $(CFLAGS) $^ -O3 -o $@

sinkhorn.o: sinkhorn.c ../util/generator.h
	$(CC) $(CFLAGS) -O3 -c $< -o $@

clean: 
	rm *.o sinkhorn