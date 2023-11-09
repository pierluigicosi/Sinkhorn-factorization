#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


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
void sinkhorn_knopp(double *matrix, int size, double epsilon, int max_iterations) {
    // Allocazione di memoria per i vettori di scaling D e D'
    double *D = (double *)malloc(size * sizeof(double));
    double *D_prime = (double *)malloc(size * sizeof(double));
    double *approx_matrix = (double *)malloc(size * size * sizeof(double));

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
                row_sum += matrix[i * size + j] * D_prime[j];
            }           
            D[i] = 1.0 / row_sum;
        }
        
        // Scaling delle colonne di D'
        for (int j = 0; j < size; j++) {
            double col_sum = 0.0;
            for (int i = 0; i < size; i++) {
                col_sum += matrix[i * size + j] * D[i];
            }
            D_prime[j] = 1.0 / col_sum;
        }


        // Verifica della convergenza
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                approx_matrix[i * size + j] = D[i] * matrix[i * size + j] * D_prime[j];
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

    printf("\nMatrice D':\n");
    for (int i = 0; i < size; i++) {
        printf("%lf ", D_prime[i]);
    }

    double *inv_D = reverse_array(D,size);
    double *inv_D_prime = reverse_array(D_prime,size);

    // Controllo se è doubly stochastic

    //somma delle righe
    printf("\nSomma delle righe:\n");
    for (int i = 0; i < size; i++) {
        float sum=0;
        for (int j = 0; j < size; j++) {
            sum+=approx_matrix[i*size+j];
            printf("%f ",approx_matrix[i * size + j]);
        }
        printf("Somma: %f\n",sum);
    }

    //somma delle colonne
    printf("\nSomma delle colonne:\n");
    for (int i = 0; i < size; i++) {
        float sum=0;
        for (int j = 0; j < size; j++) {
            sum+=approx_matrix[j*size+i];
        }
        printf("%f ",sum);
    }

    // Deallocazione della memoria
    free(D);
    free(D_prime);
    free(approx_matrix);
}


int main() {
    int size = 10; // Dimensione della matrice quadrata
    double matrix[10][10] = { // Inserire la matrice con elementi positivi
        {1, 2, 1, 2, 0.2,1, 2, 1, 2, 0.2},
        {5, 6, 0.4, 0.3, 0.1,5, 6, 0.4, 0.3, 0.1},
        {5, 6, 0.4, 0.3, 0.1,5, 6, 0.4, 0.3, 0.1},
        {5, 6, 0.4, 0.3, 0.1,5, 6, 0.4, 0.3, 0.1},
        {2.0, 7, 10, 0.1, 1.0,2.0, 7, 10, 0.1, 1.0},
        {1, 2, 1, 2, 0.2,1, 2, 1, 2, 0.2},
        {5, 6, 0.4, 0.3, 0.1,5, 6, 0.4, 0.3, 0.1},
        {5, 6, 0.4, 0.3, 0.1,5, 6, 0.4, 6, 0.1},
        {5, 6, 0.4, 0.3, 0.1,5, 6, 0.4, 0.3, 0.1},
        {2.0, 7, 10, 0.1, 1.0,2.0, 7, 10, 0.1, 1.0}
    };

    double epsilon = 1e-6; // Soglia di convergenza
    int max_iterations = 50; // Numero massimo di iterazioni

    clock_t start = clock();

    // Esecuzione dell'algoritmo di Sinkhorn-Knopp
    sinkhorn_knopp(&matrix[0][0], size, epsilon, max_iterations);

    clock_t end = clock();

    printf("\nTempo di esecuzione: %f", (double)(end-start)/CLOCKS_PER_SEC);

    return 0;
}
