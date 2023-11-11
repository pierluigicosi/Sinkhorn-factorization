CC = gcc
CFLAGS = -g

sinkhorn: sinkhorn.o ../util/generator.o
	$(CC) $(CFLAGS) $^ -o $@

sinkhorn.o: sinkhorn.c ../util/generator.h
	$(CC) $(CFLAGS) -c $< -o $@

clean: 
	rm *.o sinkhorn