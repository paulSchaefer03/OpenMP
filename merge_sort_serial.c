#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

#define N 100000000
#define RUNDEN 1


/* #define TYP double
#define FORMAT "%f"
#define IST_INT 0 */
#define TYP int
#define FORMAT "%d" 
#define IST_INT 1

//Klassische Merge Sort(orientiert an https://www.geeksforgeeks.org/iterative-merge-sort/) Rekursiv
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

//Iterative Merge Sort (orientiert an https://www.geeksforgeeks.org/iterative-merge-sort/ aber besser)
void iterTeilenUndSortiern(TYP * A) {
    TYP *B = (TYP *)malloc(N * sizeof(TYP));
    if (B == NULL) {
        fprintf(stderr, "Speicherallocation fehlgeschlagen\n");
        exit(1);
    }

    TYP *orginal = A;
    TYP *zws = B;
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
            TYP *temp = orginal; orginal = zws; zws = temp;
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
        TYP *temp = orginal;
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


void befullenArray(TYP *A) {
    double start = (double)clock() / CLOCKS_PER_SEC;
    // Initialisierung des Arrays mit Zufallszahlen
    srand(start);
    for (int i = 0; i < N; i++) {

#if IST_INT
        A[i] = (rand() % 20000001) - 10000000;
#else
        A[i] = ((rand() % 20000001) - 10000000) + ((double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
#endif
        //A[i] = N - i;
    }

    double ende = (double)clock() / CLOCKS_PER_SEC;
    //printf("Zeit zur Initialisierung des Arrays: %f Sekunden\n", ende - start);
}

void ausgabeArray(TYP *A) {
    printf("Array:\n");
    for (int i = 0; i < N; i++) {
        printf(FORMAT " ", A[i]);
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
        befullenArray(A);
        //ausgabeArray(A);
        start = (double)clock() / CLOCKS_PER_SEC;
        rekTeilenUndSortiern(A, 0, (N-1));
        ende = (double)clock() / CLOCKS_PER_SEC;
        //ausgabeArray(A);
        rekGlobalTime += (ende - start);
        //printf("[REKURSIV] Zeit zur Sortierung des Arrays: %f Sekunden\n", ende - start);
        if(!checkSort(A)){
            fprintf(stderr, "Array ist nicht sortiert\n");
            ausgabeArray(A);
            exit(1);
        }

        befullenArray(A);
        //ausgabeArray(A);
        start = (double)clock() / CLOCKS_PER_SEC;
        iterTeilenUndSortiern(A);
        ende = (double)clock() / CLOCKS_PER_SEC;
        //ausgabeArray(A);
        iterGlobalTime += (ende - start);
        //printf("[ITERATIV] Zeit zur Sortierung des Arrays: %f Sekunden\n", ende - start);
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

//1000 Runden, 1000000 Elemente: (MIT gcc -O3 -march=native -funroll-loops)
//Durchschnittliche Zeit für rekursive Sortierung: 0.094917 Sekunden
//Durchschnittliche Zeit für iterative Sortierung: 0.072731 Sekunden

//1000 Runden, 10000000 Elemente: (OHNE gcc -O3 -march=native -funroll-loops)
//Durchschnittliche Zeit für rekursive Sortierung: 0.142225 Sekunden
//Durchschnittliche Zeit für iterative Sortierung: 0.101896 Sekunden