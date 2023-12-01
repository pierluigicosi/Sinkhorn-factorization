#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../util/generator.h"


// Funzione per calcolare la norma L2 tra due vettori
double l2_norm(double *vec1, double *vec2, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        double diff = vec1[i] - vec2[i];
        sum += diff * diff;
    }
    return sqrt(sum);
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


// Funzione per eseguire l'algoritmo di Sinkhorn-Knopp
void sinkhorn_knopp(double **matrix, int size, double epsilon, int max_iterations) {

    // Allocazione di memoria per i vettori di scaling D1 e D1'
    double *D1 = malloc(size * sizeof(double));
    double *D2 = malloc(size * sizeof(double));
    double** approx_matrix = malloc(size * sizeof(double*));
    for (int i = 0; i < size; i++) {
        approx_matrix[i] = malloc(size * sizeof(double));
    }
    // Allocazione vettori per somma colonne e righe ad ogni iterazione
    double *col_sums = malloc(size * sizeof(double));
    double *row_sums = malloc(size * sizeof(double));
    double *ones = malloc(size * sizeof(double));

    // Inizializzazione dei vettori D1 e D1' con valori tutti uguali a 1
    for (int i = 0; i < size; i++) {
        D1[i] = 1.0;
        D2[i] = 1.0;
        ones[i]=1.0;
    }

    // Variabile per conteggiare il numero di iterazioni
    int iteration = 0;

    // Ciclo di iterazione dell'algoritmo
    while (iteration < max_iterations) {

        if (iteration%2==0){
            // Normalizza le righe
            for (int i = 0; i < size; i++) {
                double row_sum = 0.0;
                for (int j = 0; j < size; j++) {
                    row_sum += matrix[i][j] * D2[j];
                }           
                D1[i] = 1.0 / row_sum;
            }
        }
        else {
            // Normalizza le colonne
            for (int j = 0; j < size; j++) {
                double col_sum = 0.0;
                for (int i = 0; i < size; i++) {
                    col_sum += matrix[i][j] * D1[i];
                }
                D2[j] = 1.0 / col_sum;
            }
        } 

        /******************** Verifica convergenza **********************/
        for (int i = 0; i < size; i++) {
            row_sums[i]=0.0;
            for (int j = 0; j < size; j++) {
                approx_matrix[i][j] = D1[i] * matrix[i][j] * D2[j];
                row_sums[i]+=approx_matrix[i][j];
            }
        }

        calcolaSommaColonne(approx_matrix,col_sums,size);

        double err_row=l2_norm(row_sums,ones,size);
        double err_col=l2_norm(col_sums,ones,size);

        printf("Errore sulle righe: %f sulle colonne: %f; \n",err_row,err_col);

        if (err_row < epsilon && err_col <epsilon) {
            break;
        }

        iteration++;
    }

    /************** Controlla se è doubly stochastic ***************/
    printf("\n\nLa somma delle righe: ");
    for (int i = 0; i < size; i++) {
        //printf("%.3f ",row_sums[i]);
    }
    
    printf("\n\nLa somma delle colonne: ");
    for (int i = 0; i < size; i++) {
         //printf("%.3f ",col_sums[i]);
    }
    

    // Deallocazione della memoria
    free(D1);
    free(D2);
    for (int i = 0; i < size; i++) {
        free(approx_matrix[i]);
    }
    free(approx_matrix);
    free(row_sums);
    free(col_sums);
    free(ones);
}


int main(int argc, char *argv[]) {
    int size = 100;               // Dimensione della matrice quadrata
    double epsilon = 1e-6;      // Soglia di convergenza
    int max_iterations = 100;   // Numero massimo di iterazioni

    scriviNumeriSuFileTxt(argv[1],size*size);  // Scrive matrice random su file

    FILE* file = fopen(argv[1], "r");
    if (!file) {
        perror("Errore durante l'apertura del file");
        exit(EXIT_FAILURE);
    }

    // Alloca matrice
    double** matrix = malloc(size * sizeof(double*));
    for (int i = 0; i < size; i++) {
        matrix[i] = malloc(size * sizeof(double));
    }

    // Legge i valori dal file
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            fscanf(file, "%lf", &matrix[i][j]);
        }
    }

    // Chiude il file
    fclose(file);

    clock_t start = clock();

    // Esecuzione dell'algoritmo di Sinkhorn-Knopp
    sinkhorn_knopp(matrix, size, epsilon, max_iterations);

    clock_t end = clock();

    printf("\nSinkhorn, matrix size %d: %6.5f seconds", size,(double)(end-start)/CLOCKS_PER_SEC);

    // Dealloca la matrice
    for (int i = 0; i < size; i++) {
        free(matrix[i]);
    }
    free(matrix);

    return 0;
}
