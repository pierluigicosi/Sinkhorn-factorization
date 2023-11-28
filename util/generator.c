#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Definizione della struct per il numero razionale
typedef struct {
    int numeratore;
    int denominatore;
} NumeroRazionale;

// Funzione per generare un numero razionale casuale tra 0 e 1
NumeroRazionale generaNumeroRazionaleCasuale() {
    NumeroRazionale numero;

    // Genera numeratore e denominatore casuali
    numero.numeratore = rand() % 1001;  // Numeratore tra 1 e 100
    numero.denominatore = rand() % 1000 +1;  // Denominatore tra 1 e 100

    return numero;
}

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
        NumeroRazionale numeroRazionale = generaNumeroRazionaleCasuale();
        double numero = (double)numeroRazionale.numeratore/numeroRazionale.denominatore;
        fprintf(file, "%lf ", numero);
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
        NumeroRazionale numeroRazionale = generaNumeroRazionaleCasuale();
        double numero = (double)numeroRazionale.numeratore/numeroRazionale.denominatore;
        fwrite(&numero, sizeof(int), 1, file);
    }

    // Chiudi il file
    fclose(file);
}
