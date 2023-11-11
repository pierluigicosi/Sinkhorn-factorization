all: generator.o MyMPI.o

generator.o : generator.c generator.h
	gcc -g -c $< -o $@

MyMPI.o: MyMPI.c MyMPI.h
	mpicc -g -c $< -o $@

clean: 
	rm generator.o MyMPI.o