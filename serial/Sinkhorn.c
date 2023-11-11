#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../util/generator.h"


double* reverse_array(const double arr[], int size) {
    double* reversed_array = malloc(size * sizeof(double));

    if (reversed_array == NULL) {
        printf("Errore: memoria non disponibile.\n");
        exit(1);
    }

    int end = size - 1;
    for (int i = 0; i < size; i++) {
        reversed_array[i] = arr[end - i];
    }

    return reversed_array;
}

// Funzione per calcolare la norma L2 tra due vettori
double l2_norm(double *vec1, double *vec2, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        double diff = vec1[i] - vec2[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

// Funzione per eseguire l'algoritmo di Sinkhorn-Knopp
void sinkhorn_knopp(double **matrix, int size, double epsilon, int max_iterations) {
    // Allocazione di memoria per i vettori di scaling D e D'
    double *D = (double *)malloc(size * sizeof(double));
    double *D_prime = (double *)malloc(size * sizeof(double));
    double** approx_matrix = malloc(size * sizeof(double*));
    for (int i = 0; i < size; i++) {
        approx_matrix[i] = malloc(size * sizeof(double));
    }

    // Inizializzazione dei vettori D e D' con valori tutti uguali a 1
    for (int i = 0; i < size; i++) {
        D[i] = 1.0;
        D_prime[i] = 1.0;
    }

    // Variabile per conteggiare il numero di iterazioni
    int iteration = 0;

    // Ciclo di iterazione dell'algoritmo
    while (iteration < max_iterations) {

        // Scaling delle righe di D
        for (int i = 0; i < size; i++) {
            double row_sum = 0.0;
            for (int j = 0; j < size; j++) {
                row_sum += matrix[i][j] * D_prime[j];
            }           
            D[i] = 1.0 / row_sum;
        }
        
        // Scaling delle colonne di D'
        for (int j = 0; j < size; j++) {
            double col_sum = 0.0;
            for (int i = 0; i < size; i++) {
                col_sum += matrix[i][j] * D[i];
            }
            D_prime[j] = 1.0 / col_sum;
        }


        // Verifica della convergenza
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                approx_matrix[i][j] = D[i] * matrix[i][j] * D_prime[j];
            }
        }
        double norm = l2_norm(D,D_prime,size);

        printf("%f ",norm);

        if (norm < epsilon) {
            break;
        }

        iteration++;
    }

    // Stampa delle matrici diagonali D e D'
    printf("\nMatrice D:\n");
    for (int i = 0; i < size; i++) {
        printf("%lf ", D[i]);
    }

    printf("\n\nMatrice D':\n");
    for (int i = 0; i < size; i++) {
        printf("%lf ", D_prime[i]);
    }

    double *inv_D = reverse_array(D,size);
    double *inv_D_prime = reverse_array(D_prime,size);

    // Controllo se è doubly stochastic

    // Somma delle righe
    printf("\n\nSomma delle righe:\n");
    for (int i = 0; i < size; i++) {
        float sum=0;
        for (int j = 0; j < size; j++) {
            sum+=approx_matrix[i][j];
            //printf("%f ",approx_matrix[i][j]);
        }
        printf("%f ",sum);
    }

    // Somma delle colonne
    printf("\n\nSomma delle colonne:\n");
    for (int i = 0; i < size; i++) {
        float sum=0;
        for (int j = 0; j < size; j++) {
            sum+=approx_matrix[j][i];
        }
        printf("%f ",sum);
    }

    // Deallocazione della memoria
    free(D);
    free(D_prime);
    free(approx_matrix);
}


int main(int argc, char *argv[]) {
    int size = 100;            // Dimensione della matrice quadrata
    double epsilon = 1e-6;     // Soglia di convergenza
    int max_iterations = 50;   // Numero massimo di iterazioni

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

    printf("\nTempo di esecuzione: %f", (double)(end-start)/CLOCKS_PER_SEC);

    // Dealloca la matrice
    for (int i = 0; i < size; i++) {
        free(matrix[i]);
    }
    free(matrix);

    return 0;
}
