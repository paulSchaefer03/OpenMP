#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdbool.h>
#include <math.h>

#define N 1000000000
#define RUNDEN 1
#define NUM_THREADS 12


//#define MULTITHREAD_GRENZE 1000000 //Cluster
#define MULTITHREAD_GRENZE 10000   //Lokal

#define MERGE_GRENZE 8192 //Lokal
//#define MERGE_GRENZE 65536 //Cluster

/* #define TYP double
#define FORMAT "%f"
#define IST_INT 0 */
#define TYP int
#define FORMAT "%d" 
#define IST_INT 1

//Klassische Merge Sort, orientiert an https://www.geeksforgeeks.org/merge-sort/ Rekursiv
void rekTeilenUndSortiern(TYP *A, long start, long ende) {
    if (start < ende) {
        long mitte = (long)round((start + ende) / 2.0);
        if ((ende - start) > MERGE_GRENZE) {
            #pragma omp task shared(A) firstprivate(start, mitte, ende)
            rekTeilenUndSortiern(A, start, mitte - 1);
            #pragma omp task shared(A) firstprivate(mitte, ende)
            rekTeilenUndSortiern(A, mitte, ende);
            #pragma omp taskwait
        } else {
            rekTeilenUndSortiern(A, start, mitte - 1);
            rekTeilenUndSortiern(A, mitte, ende);
        }
            
        long groeßeLinks = mitte - start;
        long groeßeRechts = ende - mitte + 1;

        TYP *ALinks = (TYP *)malloc(groeßeLinks * sizeof(TYP));
        TYP *ARechts = (TYP *)malloc(groeßeRechts * sizeof(TYP));

        if (ALinks == NULL || ARechts == NULL) {
            fprintf(stderr, "Speicher allocation fehlgeschlagen\n");
            exit(1);
        }

        //Temp Arrays
        for (long i = 0; i < groeßeLinks; i++) ALinks[i] = A[start + i];
        for (long i = 0; i < groeßeRechts; i++) ARechts[i] = A[mitte + i];

        long i = 0, j = 0;

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

    #pragma omp parallel num_threads(NUM_THREADS) if(N >= MULTITHREAD_GRENZE)
    {
        for (long width = 1; width < N; width *= 2) {
            TYP *localOrg = buffer[curr];
            TYP *localZws = buffer[1 - curr];

            if (width == 1) {
                #pragma omp for schedule(guided, 1) nowait
                for (long i = 0; i < N / 2; i++) {
                    long idx = 2 * i;
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
                for (long i = 0; i < N; i += 2 * width) {
                    long start = i;
                    long mitte = i + width;
                    long ende = (i + 2 * width < N) ? i + 2 * width : N;

                    if (mitte >= N) {
                        for (long k = start; k < ende; k++)
                            localZws[k] = localOrg[k];
                        continue;
                    }

                    long links = start, rechts = mitte;
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
            for (long i = 0; i < N; i++) {
                A[i] = buffer[curr][i];
            }
        }
    }

    free(B);
}



void befullenArray(TYP *A, int init_basis) {
    
    #pragma omp for
    for (long i = 0; i < N; i++) {
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
    for (long i = 0; i < N; i++) {
        printf(FORMAT " \n", A[i]);
    }
    printf("\n");
}

bool checkSort(TYP *A) {
    for (long i = 0; i < N - 1; i++) {
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

        #pragma omp parallel num_threads(NUM_THREADS) if(N >= MULTITHREAD_GRENZE)
        {

            #pragma omp single
            start = omp_get_wtime();
            
            int init_basis = (int)start + 56;
            befullenArray(A, init_basis);
            #pragma omp barrier
            #pragma omp single
            {
                ende = omp_get_wtime();
                //printf("Zeit zur Initialisierung 1 des Arrays: %f Sekunden\n", ende - start);
                //ausgabeArray(A);
            }

            #pragma omp single
            {
                start = omp_get_wtime();
                rekTeilenUndSortiern(A, 0, (N-1));
                //ausgabeArray(A);
                ende = omp_get_wtime();
                rekGlobalTime += (ende - start);
                if(!checkSort(A)){
                    fprintf(stderr, "Array ist nicht sortiert\n");
                    //ausgabeArray(A);
                    exit(1);
                }
                start = omp_get_wtime();
            }

            init_basis = (int)start + 56;
            befullenArray(A, init_basis);
            #pragma omp barrier
            #pragma omp single
            {
                ende = omp_get_wtime();
                //printf("Zeit zur Initialisierung 2 des Arrays: %f Sekunden\n", ende - start);
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
            //ausgabeArray(A, N);
            exit(1);
        }
            
    }
    free(A);

    printf("Durchschnittliche Zeit für rekursive Sortierung mit %ld Elementen in %d Runden mit %s Array und %d Threads: %f Sekunden\n", (long)N, RUNDEN, IST_INT ? "int" : "double", NUM_THREADS, rekGlobalTime / RUNDEN);
    printf("Durchschnittliche Zeit für iterative Sortierung mit %ld Elementen in %d Runden mit %s Array und %d Threads: %f Sekunden\n", (long)N, RUNDEN, IST_INT ? "int" : "double", NUM_THREADS, iterGlobalTime / RUNDEN);

    return 0;
}
