#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../util/generator.h"
#include "../util/MyMPI.h"

// Funzione per calcolare la norma L2 tra due vettori
double l2_norm(double *vec1, double *vec2, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        double diff = vec1[i] - vec2[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

int main(int argc, char **argv)
{
    int rank;                   // Numero del processo attivo
    int num_procs;              // Numero totale processi
    double *matrix=NULL;        // Matrice
    int size=100;               // Dimensione matrice
    double epsilon = 1e-6;      // Soglia di convergenza
    int max_iterations = 100;   // Numero massimo di iterazioni
    double seconds,max_seconds; // Tempo
    int flag=1;                 // Flag break ciclo
    int iter=0;                 // Iterazioni ciclo sinkhorn


    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (argc!=2){
        terminate(rank,"Inserire nome file matrice");
    }

    if (num_procs>size){
        terminate(rank,"Il numero dei processi dev' essere minore della size");
    }

    // Scrittura matrice su file e lettura di essa
    if (rank==0){
        scriviNumeriSuFileTxt(argv[1],size*size);

        FILE* file = fopen(argv[1], "r");
        if (!file) {
            terminate(rank,"Errore nell'apertura del file");
        }

        // Alloca matrice
        matrix = malloc(size * size * sizeof(double));

        // Legge i valori dal file
        for (int i = 0; i < size * size; i++) {
            fscanf(file, "%lf", &matrix[i]);
        }
        // Chiude il file
        fclose(file);
    }

    MPI_Barrier (MPI_COMM_WORLD);
    seconds = - MPI_Wtime();

    /**************************************** Decomposizione matrice e vettori **************************************************/ 
    // displacement e count per scatter matrice
    int *displscatter;
    int *sendcount;
    divide_matrix(rank,num_procs,size,&sendcount,&displscatter);
    // displacement e count per gather somma righe
    int *displsgather;
    int *recvcount;
    divide_array(rank,num_procs,size,&recvcount,&displsgather);

    // Matrix row block decomposition 
    int localRows = BLOCK_SIZE(rank,num_procs,size);

    // Alloca memoria per porzione locale della matrice 
    double* localMatrix = my_malloc(rank,localRows * size * sizeof(double*));

    // Divide la matrice fra i processi
    MPI_Scatterv(matrix,sendcount,displscatter,MPI_DOUBLE,localMatrix,localRows*size,MPI_DOUBLE,0,MPI_COMM_WORLD);

    // Inizializzazione vettori di scaling
    double *D2,*ones,*row_sums,*col_sums;
    D2 = my_malloc(rank,size * sizeof(double));
    ones=malloc(size * sizeof(double));
    for (int i = 0; i < size; i++) {
        D2[i] = 1.0;
        ones[i]=1.0;
    }

    // Vettori locali
    double *local_D1 = my_malloc (rank,localRows * sizeof(double));
    double *local_D2 = my_malloc (rank,size * sizeof(double));
    double *local_col_sums = my_malloc(rank,size * sizeof(double));
    double *local_row_sums = my_malloc(rank,localRows * sizeof(double));

    if(rank==0){
        col_sums = malloc(size * sizeof(double));
        row_sums = malloc(size * sizeof(double));
    }

    /**************************************** Esecuzione Sinkhorn-Knoppin parallelo **************************************************/ 
    while(iter<max_iterations) {
        
        /********************************* Normalizza righe e colonne *********************************/
        if(iter%2==0){
            // Normalizza le righe
            for (int i = 0; i < localRows; i++) {
                local_D1[i] = 0.0;
                for (int j = 0; j < size; j++) {
                    local_D1[i] += localMatrix[i*size + j] * D2[j];
                }  
                local_D1[i] = 1.0/local_D1[i];
                
            }

        }else{
            // Normalizza le colonne
            for (int i = 0; i < size; i++) {
                local_D2[i] = 0.0;
                for (int j = 0; j < localRows; j++){
                    local_D2[i] += localMatrix[j*size+i] * local_D1[j];
                }
            }

            // Allreduce local D2 across all processes
            MPI_Allreduce(local_D2,D2,size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            for (int i = 0; i < size; i++) {
                D2[i]=1.0/D2[i];
            }
        }

        /************************************** Verifica convergenza **************************************/

        // Calcolo somma righe
        for (int i = 0; i < localRows; i++) {
            local_row_sums[i] = 0;
            for (int j = 0; j < size; j++) {
                local_row_sums[i] += local_D1[i] * localMatrix[i*size+j] * D2[j];
            }
        }

        MPI_Gatherv(local_row_sums,localRows,MPI_DOUBLE,row_sums,recvcount,displsgather,MPI_DOUBLE,0,MPI_COMM_WORLD);

        // Calcolo somma colonne
        for (int i = 0; i < size; i++) {
            local_col_sums[i] = 0;  
            for (int j = 0; j < localRows; j++) {
                local_col_sums[i] += local_D1[j] * localMatrix[j*size+i] * D2[i];
            }
        }

        MPI_Reduce(local_col_sums,col_sums,size,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

            
        if(rank==0){
            double err_row=l2_norm(row_sums,ones,size);
            double err_col=l2_norm(col_sums,ones,size);

            printf("Errore sulle righe: %f sulle colonne: %f; \n",err_row,err_col);

            if (err_row < epsilon && err_col <epsilon) {
                flag=0;
            }
        }

        MPI_Bcast(&flag,1,MPI_INT,0,MPI_COMM_WORLD);
        if(!flag) break;

        iter++;

    }


    // Stampa risultato 
    if (rank == 0) {
        printf("\n\nLa somma delle righe: ");
        for (int i = 0; i < size; i++) {
            //printf("%.3f ",row_sums[i]);
        }
        
        printf("\n\nLa somma delle colonne: ");
        for (int i = 0; i < size; i++) {
            //printf("%.3f ",col_sums[i]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    seconds += MPI_Wtime();
    MPI_Reduce (&seconds, &max_seconds, 1, MPI_DOUBLE, MPI_MAX, 0,  MPI_COMM_WORLD);
    if (!rank) printf("\nSinkhorn, matrix size %d, %d processes: %6.5f seconds",size,num_procs,max_seconds);

    // Libero memoria
    if(rank==0){
        free(matrix);
        free(col_sums);
        free(row_sums);
    }
    free(D2);
    free(local_D1);
    free(local_D2);
    free(localMatrix);
    free(local_row_sums);
    free(local_col_sums);
    free(ones);

    MPI_Finalize();
    return 0;
}