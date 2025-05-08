#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdbool.h>
#include <math.h>

#define N 1000000
#define ROUNDEN 1000

//Klassische Merge Sort(orientiert an https://www.geeksforgeeks.org/iterative-merge-sort/) Rekursiv
void rekTeilenUndSortiern(double *A, int start, int ende) {
    if (start < ende) {
        int mitte = (int)round((start + ende) / 2.0);
        rekTeilenUndSortiern(A, start, mitte - 1);
        rekTeilenUndSortiern(A, mitte, ende);

        int groeßeLinks = mitte - start;
        int groeßeRechts = ende - mitte + 1;

        double *ALinks = (double *)malloc(groeßeLinks * sizeof(double));
        double *ARechts = (double *)malloc(groeßeRechts * sizeof(double));

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

//Iterative Merge Sort (orientiert an https://www.geeksforgeeks.org/iterative-merge-sort/ aber besser)
void iterTeilenUndSortiern(double * A) {
    double *B = (double *)malloc(N * sizeof(double));
    if (B == NULL) {
        fprintf(stderr, "Speicherallocation fehlgeschlagen\n");
        exit(1);
    }

    double *orginal = A;
    double *zws = B;
    for (int width = 1; width < N; width *= 2) {
        if (width == 1) {
            for (int i = 0; i + 1 < N; i += 2) {
                if (orginal[i] <= orginal[i + 1]) {
                    zws[i] = orginal[i];
                    zws[i + 1] = orginal[i + 1];
                } else {
                    zws[i] = orginal[i + 1];
                    zws[i + 1] = orginal[i];
                }
            }
            if (N % 2 == 1) zws[N - 1] = orginal[N - 1];
            double *temp = orginal; orginal = zws; zws = temp;
            continue;
        }

        for (int i = 0; i < N; i += 2 * width) {
            int start = i;
            int mitte = i + width;
            int ende = (i + 2 * width < N) ? i + 2 * width : N;

            if (mitte >= N) {
                for (int k = start; k < ende; k++) zws[k] = orginal[k];
                continue;
            }

            int links = start, rechts = mitte;

            for (; links < mitte && rechts < ende; start++) {
                if (orginal[links] <= orginal[rechts]) {
                    zws[start] = orginal[links++];
                } else {
                    zws[start] = orginal[rechts++];
                }
            }

            //Rest Links

            for (; links < mitte; start++) zws[start] = orginal[links++];
            for (; rechts < ende; start++) zws[start] = orginal[rechts++];

        }

        // Tauschen
        double *temp = orginal;
        orginal = zws;
        zws = temp;
    }

    // 
    if (orginal != A) {
        for (int i = 0; i < N; i++)
            A[i] = orginal[i];
    }

    free(B);
}


void befullenArray(double *A) {
    double start = omp_get_wtime();
    // Initialisierung des Arrays mit Zufallszahlen
    srand(start);
    for (int i = 0; i < N; i++) {
        A[i] = ((double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
    }

    double ende = omp_get_wtime();
    printf("Zeit zur Initialisierung des Arrays: %f Sekunden\n", ende - start);
}

void ausgabeArray(double *A) {
    printf("Array:\n");
    for (int i = 0; i < N; i++) {
        printf("%f ", A[i]);
    }
    printf("\n");
}

bool checkSort(double *A) {
    for (int i = 0; i < N - 1; i++) {
        if (A[i] > A[i + 1]) {
            return false;
        }
    }
    return true;
}

int main() {
    double *A = (double * )malloc(N * sizeof(double));
    double rekGlobalTime, iterGlobalTime, start, ende;

    if(A == NULL) {
        fprintf(stderr, "Speicher allocation fehlgeschalgen\n");
        return 1;
    }

    // Init Array
    for(int i = 0; i < ROUNDEN; i++) {
        befullenArray(A);
        start = omp_get_wtime();
        rekTeilenUndSortiern(A, 0, (N-1));
        ende = omp_get_wtime();
        rekGlobalTime += (ende - start);
        printf("[REKURSIV] Zeit zur Sortierung des Arrays: %f Sekunden\n", ende - start);
        if(!checkSort(A)){
            fprintf(stderr, "Array ist nicht sortiert\n");
            ausgabeArray(A);
            exit(1);
        }

        befullenArray(A);
        start = omp_get_wtime();
        iterTeilenUndSortiern(A);
        ende = omp_get_wtime();
        iterGlobalTime += (ende - start);
        printf("[ITERATIV] Zeit zur Sortierung des Arrays: %f Sekunden\n", ende - start);
        if(!checkSort(A)){
            fprintf(stderr, "Array ist nicht sortiert\n");
            ausgabeArray(A);
            exit(1);
        }
    }
    free(A);

    printf("Durchschnittliche Zeit für rekursive Sortierung: %f Sekunden\n", rekGlobalTime / ROUNDEN);
    printf("Durchschnittliche Zeit für iterative Sortierung: %f Sekunden\n", iterGlobalTime / ROUNDEN);

    return 0;
}

//1000 Runden, 1000000 Elemente: (MIT gcc -O3 -march=native -funroll-loops)
//Durchschnittliche Zeit für rekursive Sortierung: 0.094917 Sekunden
//Durchschnittliche Zeit für iterative Sortierung: 0.072731 Sekunden

//1000 Runden, 10000000 Elemente: (OHNE gcc -O3 -march=native -funroll-loops)
//Durchschnittliche Zeit für rekursive Sortierung: 0.142225 Sekunden
//Durchschnittliche Zeit für iterative Sortierung: 0.101896 Sekunden