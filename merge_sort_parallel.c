#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdbool.h>
#include <math.h>

#define N 100000000
#define RUNDEN 100
#define NUM_THREADS 12


#define TYP double
#define FORMAT "%f"
#define IST_INT 0
/* #define TYP int
#define FORMAT "%d" 
#define IST_INT 1 */

//Klassische Merge Sort, orientiert an https://www.geeksforgeeks.org/merge-sort/ Rekursiv
void rekTeilenUndSortiern(TYP *A, int start, int ende) {
    if (start < ende) {
        int mitte = (int)round((start + ende) / 2.0);
        rekTeilenUndSortiern(A, start, mitte - 1);
        rekTeilenUndSortiern(A, mitte, ende);

        int groeßeLinks = mitte - start;
        int groeßeRechts = ende - mitte + 1;

        TYP *ALinks = (TYP *)malloc(groeßeLinks * sizeof(TYP));
        TYP *ARechts = (TYP *)malloc(groeßeRechts * sizeof(TYP));

        if (ALinks == NULL || ARechts == NULL) {
            fprintf(stderr, "Speicher allocation fehlgeschlagen\n");
            exit(1);
        }

        //Temp Arrays
        for (int i = 0; i < groeßeLinks; i++) ALinks[i] = A[start + i];
        for (int i = 0; i < groeßeRechts; i++) ARechts[i] = A[mitte + i];

        int i = 0, j = 0;

        //Sortieren
        for (; i < groeßeLinks && j < groeßeRechts; start++) {
            if (ALinks[i] <= ARechts[j]) {
                A[start] = ALinks[i++];
            } else {
                A[start] = ARechts[j++];
            }
        }

        //Rest Links
        for (; i < groeßeLinks; start++) A[start] = ALinks[i++];

        //Rest Rechts
        for (; j < groeßeRechts; start++) A[start] = ARechts[j++];

        free(ALinks);
        free(ARechts);
    }
}

//Iterative Merge Sort, orientiert an https://www.geeksforgeeks.org/iterative-merge-sort/ (aber besser)
void iterTeilenUndSortiern(TYP *A) {
    TYP *B = (TYP *)malloc(N * sizeof(TYP));
    if (B == NULL) {
        fprintf(stderr, "Speicherallocation fehlgeschlagen\n");
        exit(1);
    }

    TYP *buffer[2] = { A, B };
    int curr = 0;  // 0 = A ist Input, 1 = B ist Output

    #pragma omp parallel num_threads(NUM_THREADS)
    {
        for (int width = 1; width < N; width *= 2) {
            TYP *localOrg = buffer[curr];
            TYP *localZws = buffer[1 - curr];

            if (width == 1) {
                #pragma omp for schedule(guided, 1) nowait
                for (int i = 0; i < N / 2; i++) {
                    int idx = 2 * i;
                    if (localOrg[idx] <= localOrg[idx + 1]) {
                        localZws[idx] = localOrg[idx];
                        localZws[idx + 1] = localOrg[idx + 1];
                    } else {
                        localZws[idx] = localOrg[idx + 1];
                        localZws[idx + 1] = localOrg[idx];
                    }
                }
                #pragma omp single nowait
                if (N % 2 == 1)
                    localZws[N - 1] = localOrg[N - 1];
            } else {
                #pragma omp for schedule(guided, 1) nowait
                for (int i = 0; i < N; i += 2 * width) {
                    int start = i;
                    int mitte = i + width;
                    int ende = (i + 2 * width < N) ? i + 2 * width : N;

                    if (mitte >= N) {
                        for (int k = start; k < ende; k++)
                            localZws[k] = localOrg[k];
                        continue;
                    }

                    int links = start, rechts = mitte;
                    for (; links < mitte && rechts < ende; start++) {
                        if (localOrg[links] <= localOrg[rechts]) {
                            localZws[start] = localOrg[links++];
                        } else {
                            localZws[start] = localOrg[rechts++];
                        }
                    }

                    for (; links < mitte; start++) localZws[start] = localOrg[links++];
                    for (; rechts < ende; start++) localZws[start] = localOrg[rechts++];
                }
            }

            // Umschalten des Puffers nach abgeschlossener Stufe
            #pragma omp single
            {
                curr = 1 - curr;
            }

            #pragma omp barrier
        }

        // Finaler Zustand in A zurückschreiben
        if (curr != 0) {
            #pragma omp for schedule(static)
            for (int i = 0; i < N; i++) {
                A[i] = buffer[curr][i];
            }
        }
    }

    free(B);
}



void befullenArray(TYP *A, int init_basis) {
    
    #pragma omp for
    for (int i = 0; i < N; i++) {
        unsigned int init = init_basis + omp_get_thread_num() + i;
#if IST_INT
        A[i] = (rand_r(&init) % 20000001) - 10000000;
#else
        A[i] = ((rand_r(&init) % 20000001) - 10000000) + ((double)rand_r(&init) / RAND_MAX) * (rand_r(&init) % 2 == 0 ? 1 : -1);
#endif
        //A[i] = N - i;
    }

}

void ausgabeArray(TYP *A) {
    printf("Array:\n");
    for (int i = 0; i < N; i++) {
        printf(FORMAT " \n", A[i]);
    }
    printf("\n");
}

bool checkSort(TYP *A) {
    for (int i = 0; i < N - 1; i++) {
        if (A[i] > A[i + 1]) {
            return false;
        }
    }
    return true;
}

int main() {
    TYP *A = (TYP * )malloc(N * sizeof(TYP));
    double rekGlobalTime, iterGlobalTime, start, ende;

    if(A == NULL) {
        fprintf(stderr, "Speicher allocation fehlgeschalgen\n");
        return 1;
    }

    // Init Array
    for(int i = 0; i < RUNDEN; i++) {
/*         befullenArray(A);
        //ausgabeArray(A);
        start = omp_get_wtime();
        rekTeilenUndSortiern(A, 0, (N-1));
        ende = omp_get_wtime();
        //ausgabeArray(A);
        rekGlobalTime += (ende - start);
        //printf("[REKURSIV] Zeit zur Sortierung des Arrays: %f Sekunden\n", ende - start);
        if(!checkSort(A)){
            fprintf(stderr, "Array ist nicht sortiert\n");
            ausgabeArray(A);
            exit(1);
        } */

        #pragma omp parallel num_threads(NUM_THREADS)
        {

            #pragma omp single
            {
                start = omp_get_wtime();
            }
            int init_basis = (int)start + 56;
            befullenArray(A, init_basis);
            #pragma omp barrier
            #pragma omp single
            {
                ende = omp_get_wtime();
                //printf("Zeit zur Initialisierung des Arrays: %f Sekunden\n", ende - start);
                //ausgabeArray(A);
            }
        }    
                
        start = omp_get_wtime();

        iterTeilenUndSortiern(A);

        ende = omp_get_wtime();
        iterGlobalTime += (ende - start);
        //printf("Zeit zur Sortierung des Arrays: %f Sekunden\n", ende - start);
        //ausgabeArray(A);
        if(!checkSort(A)){
            fprintf(stderr, "Array ist nicht sortiert\n");
            ausgabeArray(A);
            exit(1);
        }
            
    }
    free(A);

    printf("Durchschnittliche Zeit für rekursive Sortierung mit %d Elementen in %d Runden mit %s Array: %f Sekunden\n", N, RUNDEN, IST_INT ? "int" : "double", rekGlobalTime / RUNDEN);
    printf("Durchschnittliche Zeit für iterative Sortierung mit %d Elementen in %d Runden mit %s Array: %f Sekunden\n", N, RUNDEN, IST_INT ? "int" : "double", iterGlobalTime / RUNDEN);

    return 0;
}
