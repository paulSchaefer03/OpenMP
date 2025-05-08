#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 1000000000

int main() {
    

    double *A = (double * )malloc(N * sizeof(double));

    if(A == NULL) {
        fprintf(stderr, "Speicher allocation fehlgeschalgen\n");
        return 1;
    }

    // Init Array
    double start = omp_get_wtime();
    // Initialisierung des Arrays mit Zufallszahlen
    srand(start);
    for (int i = 0; i < N; i++) {
        A[i] = ((double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
    }
    double end = omp_get_wtime();
    printf("Zeit zur Initialisierung des Arrays: %f Sekunden\n", end - start);


    // Berechen der Summe
    start = omp_get_wtime();
    double erg = 0.0;
    for (int i = 0; i < N; i++) {
        erg += A[i];
    }
    end = omp_get_wtime();
    printf("Summe der Elemente im Array: %f in %f Sekunden\n", erg, end - start);
    
    free(A);
    A = NULL; 

    return 0;
}