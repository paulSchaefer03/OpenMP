#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 1000000000

int main() {
    

    double *A = (double * )malloc(N * sizeof(double));

    if(A == NULL) {
        fprintf(stderr, "Speicher allocation fehlgeschalgen\n");
        return 1;
    }

    // Init Array
    double start = (double)clock() / CLOCKS_PER_SEC;
    // Initialisierung des Arrays mit Zufallszahlen
    srand(start);
    for (int i = 0; i < N; i++) {
        A[i] = ((double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
    }
    double end = (double)clock() / CLOCKS_PER_SEC;
    printf("Zeit zur Initialisierung des Arrays: %f Sekunden\n", end - start);


    // Berechen der Summe
    start = (double)clock() / CLOCKS_PER_SEC;
    double erg = 0.0;
    for (int i = 0; i < N; i++) {
        erg += A[i];
    }
    end = (double)clock() / CLOCKS_PER_SEC;
    printf("Summe der Elemente im Array: %f in %f Sekunden\n", erg, end - start);
    
    free(A);
    A = NULL; 

    return 0;
}