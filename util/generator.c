#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Funzione per scrivere numeri razionali su un file
void scriviNumeriSuFileTxt(const char* filetxt, int numeroNumeri) {
    FILE* file = fopen(filetxt, "w");
    if (!file) {
        perror("Errore durante l'apertura del file");
        exit(EXIT_FAILURE);
    }

    srand(time(NULL));

    // Scrivi i numeri razionali nel file di testo e in quello binario
    for (int i = 0; i < numeroNumeri; i++) {
        int numero = rand() % 10001;
        fprintf(file, "%d ", numero);
    }

    // Chiudi il file
    fclose(file);
}

// Funzione per scrivere numeri razionali su un file
void scriviNumeriSuFileBin(const char* filebin, int dim) {
    FILE* file = fopen(filebin, "w");
    if (!file) {
        perror("Errore durante l'apertura del file");
        exit(EXIT_FAILURE);
    }

    fwrite(&dim, sizeof(int), 1, file);
    fwrite(&dim, sizeof(int), 1, file);

    srand(time(NULL));

    // Scrivi i numeri razionali nel file di testo e in quello binario
    for (int i = 0; i < dim*dim; i++) {
        int numero = rand() % 10001;
        fwrite(&numero, sizeof(int), 1, file);
    }

    // Chiudi il file
    fclose(file);
}
