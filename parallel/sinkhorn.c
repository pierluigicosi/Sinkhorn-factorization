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

// Funzione per calcolare la somma delle righe di una matrice
void calcolaSommaRighe(double **matrix, double *sommaRighe,int size) {
    for (int i = 0; i < size; i++) {
        sommaRighe[i] = 0;  // Inizializza la somma per ogni riga
        for (int j = 0; j < size; j++) {
            sommaRighe[i] += matrix[i][j];
        }
    }
}

// Funzione per calcolare la somma delle colonne di una matrice
void calcolaSommaColonne(double **matrix, double *sommaColonne,int size) {
    for (int j = 0; j < size; j++) {
        sommaColonne[j] = 0;  // Inizializza la somma per ogni colonna
        for (int i = 0; i < size; i++) {
            sommaColonne[j] += matrix[i][j];
        }
    }
}

int main(int argc, char **argv)
{
    int rank;                   // Numero del processo attivo
    double *matrix=NULL;        // Matrice
    int num_procs;              // Numero totale processi
    int size=5;                 // Dimensione matrice
    double epsilon = 1e-6;      // Soglia di convergenza
    int max_iterations = 200;     // Numero massimo di iterazioni
    double seconds,max_seconds; // Tempo


    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (argc!=2){
        perror("Errore: inserire file matrice");
        exit(EXIT_FAILURE);
    }

    if (num_procs>size){
        perror("Errore: troppi processi");
        exit(EXIT_FAILURE);
    }

    // inizializzazione matrice 
    if (rank==0){
        scriviNumeriSuFileTxt(argv[1],size*size);  // scrive matrice random su file binario

        FILE* file = fopen(argv[1], "r");
        if (!file) {
            perror("Errore durante l'apertura del file");
            exit(EXIT_FAILURE);
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


    // Row block decomposition
    int localRows = BLOCK_SIZE(rank,num_procs,size);

    //printf("rank %d block size %d",rank,localRows);

    // Allocate memory for the local portion of the matrix
    double* localMatrix = my_malloc(rank,localRows * size * sizeof(double*));


    // Scatter the matrix among processes using checkerboard decomposition
    int *displscatter;
    int *sendcount;
    scatter_arrays(rank,num_procs,size,&sendcount,&displscatter);
    MPI_Scatterv(matrix,sendcount,displscatter,MPI_DOUBLE,localMatrix,localRows*size,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //printf("rank %d %lf ",rank,localMatrix[6]);

    // inizializzazione vettori di scaling
    double *D,*D_prime,*ones;
    D = my_malloc(rank,size * sizeof(double));
    D_prime = my_malloc(rank,size * sizeof(double));
    ones=malloc(size * sizeof(double));
    for (int i = 0; i < size; i++) {
        D[i] = 1.0;
        D_prime[i] = 1.0;
        ones[i]=1.0;
    }

    // esecuzione sinkhorn knopp
    MPI_Barrier (MPI_COMM_WORLD);
    seconds = - MPI_Wtime();

    int *displs;
    int *recvcount;
    create_mixed_xfer_arrays(rank,num_procs,size,&recvcount,&displs);

    double *rowsums = my_malloc (rank,size * sizeof(double));
    double *colsums = my_malloc (rank,size * sizeof(double));
    double *local_D = my_malloc (rank,localRows * sizeof(double));
    double *local_colsums = my_malloc (rank,size * sizeof(double));

    int iter=0;
    double *col_sums = malloc(size * sizeof(double));
    double *row_sums = malloc(size * sizeof(double));
    int flag=1;
    
    while(iter<max_iterations) {
        
        /**************************** Normalizza righe e colonne in parallelo *******************************/
        if(iter%2==0){
            // Normalizza le righe
            for (int i = 0; i < localRows; i++) {
                rowsums[i] = 0.0;
                for (int j = 0; j < size; j++) {
                    rowsums[i] += localMatrix[i*size + j] * D_prime[j];
                    //printf("%lf ",rowsums[i]);
                }  
                local_D[i] = 1.0/rowsums[i];
                
            }
            
            // Allreduce local D across all processes
            MPI_Gatherv(local_D,localRows,MPI_DOUBLE,D,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);

        }else{
            // Normalizza le colonne
            for (int i = 0; i < size; i++) {
                local_colsums[i] = 0.0;
                for (int j = 0; j < localRows; j++){
                    local_colsums[i] += localMatrix[j*size+i] * local_D[j];
                }
            }

            // Allreduce local D_prime across all processes
            MPI_Allreduce(local_colsums,colsums,size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            for (int i = 0; i < size; i++) {
                D_prime[i]=1.0/colsums[i];
            }
        }

        /******************** Verifica convergenza **********************/
        if (rank==0){
            double** approx_matrix = malloc(size * sizeof(double*));
            for (int i = 0; i < size; i++) {
                approx_matrix[i] = malloc(size * sizeof(double));
            }

            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    approx_matrix[i][j] = D[i] * matrix[i*size+j] * D_prime[j];
                }
            }

            //fare gather

            calcolaSommaColonne(approx_matrix,col_sums,size);
            calcolaSommaRighe(approx_matrix,row_sums,size);

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


    // TEST RISULTATO
    if (rank == 0) {
        printf("\n\nLa somma delle righe: ");
        for (int i = 0; i < size; i++) {
            printf("%.3f ",row_sums[i]);
        }
        
        printf("\n\nLa somma delle colonne: ");
        for (int i = 0; i < size; i++) {
            printf("%.3f ",col_sums[i]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    seconds += MPI_Wtime();
    MPI_Reduce (&seconds, &max_seconds, 1, MPI_DOUBLE, MPI_MAX, 0,  MPI_COMM_WORLD);
    if (!rank) printf("\nSinkhorn, matrix size %d, %d processes: %6.5f seconds",size,num_procs,max_seconds);

    // libero memoria
    if(rank==0){
        free(matrix);
    }
    free(rowsums);
    free(colsums);
    free(local_D);
    free(local_colsums);
    free(D);
    free(D_prime);
    free(localMatrix);

    MPI_Finalize();
    return 0;
}