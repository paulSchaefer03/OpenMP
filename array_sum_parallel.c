#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 1000000000
#define NUM_THREADS 12
#define ANZ_RUNDEN 1000

int main() {
    double global_zeit_reduce = 0.0;
    double global_zeit_critical = 0.0;


    for (int i = 0; i < ANZ_RUNDEN; i++) {
        double *A = (double * )malloc(N * sizeof(double));

        if(A == NULL) {
            fprintf(stderr, "Speicher allocation fehlgeschalgen\n");
            return 1;
        }
    
        double start = 0.0;
        double ende = 0.0;
        double erg = 0.0;
        double erg2 = 0.0;
    
    
        // Init  Arrays mit Zufallszahlen
        #pragma omp parallel num_threads(NUM_THREADS)
        {
            #pragma omp single
            {
                start = omp_get_wtime();
                srand(start);
            }
    
            int init = (int)start + 56 + omp_get_thread_num();
            #pragma omp for
            for (long long i = 0; i < N; i++) {
                double r = (double)rand_r(&init) / RAND_MAX;
                int s = rand_r(&init) % 2 == 0 ? 1 : -1;
                A[i] = r * s;
                //A[i] = i;
            }
    
            #pragma omp single
            {
                ende = omp_get_wtime();
                //printf("Zeit zur Initialisierung des Arrays: %f Sekunden\n", ende - start);
                start = omp_get_wtime();
            }
    
    
            // Berechen der Summe
            double partial_erg = 0.0;
            #pragma omp for
            for (long long i = 0; i < N; i++) {
                partial_erg += A[i];
            }
    
            #pragma omp critical
            {
                erg += partial_erg;
                //printf("erg: %f, partial_erg: %f (Thread %d)\n",erg, partial_erg, omp_get_thread_num());
            }
            
            #pragma omp barrier
            #pragma omp single
            {
                ende = omp_get_wtime();
                //printf("Summe der Elemente im Array(critical): %f in %f Sekunden (Thread %d)\n", erg, ende - start, omp_get_thread_num());
                global_zeit_critical += (ende - start);
            }
    
        }
    
        start = omp_get_wtime();
        erg2 = 0.0;
        #pragma omp parallel for reduction(+:erg2) num_threads(NUM_THREADS)
        for (long long i = 0; i < N; i++) {
                erg2 += A[i];
        }
        ende = omp_get_wtime();
        //printf("Summe der Elemente im Array(Reduce): %f in %f Sekunden (Thread %d)\n", erg2, ende - start, omp_get_thread_num());
        global_zeit_reduce += (ende - start);

        free(A);
        A = NULL; 
    }

    printf("Summe der Zeiten (critical) : %f Sekunden\n", global_zeit_critical / ANZ_RUNDEN);
    printf("Summe der Zeiten (reduction): %f Sekunden\n", global_zeit_reduce / ANZ_RUNDEN);


    return 0;
}