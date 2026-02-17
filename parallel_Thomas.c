#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>
#include <time.h>

#define BLOCK_LOW(id,p,n)    ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n)   (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n)   (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(id,p,n)  (((p)*((id)+1)-1)/(n))

/*
 * --------------------------------- thomas_algorithm ------------------------------------//
 * Implementazione dell'algpritmo di Thomas standard; in input richiede:
 * - N: intero che indica la dimensione dei vettori A, B, C, D e X
 * - A: puntatore all'array che contiene gli elementi sotto la diagonale principale
 * - B: puntatore all'array che contiene gli elementi della diagonale principale
 * - C: puntatore all'array che contiene gli elementi sopra la diagonale principale
 * - D: puntatore all'array che contiene i termini noti del sistema
 * - X: puntatore all'array dove verranno inserite le soluzioni del sistema
 */
void thomas_algorithm(int N, double *A, double *B, double *C, double *D, double *X) {

    double w;

    // FORWARD: L'obiettivo è quello di eliminare gli elementi A[i] cioè gli elementi al di sotto
    // della diagonale principale.
    for (int i = 1; i < N; i++) {
        w = A[i] / B[i - 1];

        B[i] = B[i] - w * C[i - 1];
        D[i] = D[i] - w * D[i - 1];
    }

    // BACKWARD: Partendo dall'ultima incognita e risalendo alla prima vengono ricavate
    // le soluzioni del sistema
    X[N - 1] = D[N - 1] / B[N - 1];

    for (int i = N - 2; i >=0; i--) {
        X[i] = (D[i] - C[i] * X[i + 1]) / B[i];
    }
}

/*
 * --------------------------------- thomas_v1 ------------------------------------//
 * Implementazione dell'algoritmo di Thomas modificato; in input richiede:
 * - N: intero che indica la dimensione dei vettori A, B, C e D
 * - A: puntatore all'array che contiene gli elementi sotto la diagonale principale
 * - B: puntatore all'array che contiene gli elementi della diagonale principale
 * - C: puntatore all'array che contiene gli elementi sopra la diagonale principale
 * - D: puntatore all'array che contiene i termini noti del sistema
 */
void thomas_v1(int N, double *A, double *B, double *C, double *D) {

    double w;

    // Vengono normalizzate le prime 2 righe per B[i]
    for (int i = 0; i < N; i++) {
        if (i == 0 || i == 1) {
            A[i] = A[i]/B[i];
            C[i] = C[i]/B[i];
            D[i] = D[i]/B[i];
            B[i] = B[i]/B[i];
        }
        else {
            /*
             * Devo eliminare elemento sotto la diagonale principale. Partendo dalla 2° riga (e non dalla 1°)
             * come in Thomas standard, si fa in modo di riscrivere ogni riga interna in funzione della variabile
             * di bordo sinistra.
             * La formula è la stessa usata nel passo di Forward in Thomas standard:
             *
             *                Ri = Ri - w * R_{i-1}   con w = (Ai)/(B_{i-1})
             *
            */

            // Oltre all'aggiornamento della riga, viene effettuata anche la normalizzazione per B[i]
            w = A[i] / B[i - 1];
            B[i] = B[i] - w * C[i - 1];
            D[i] = (D[i] - w * D[i - 1]) / B[i];
            A[i] = (- w * A[i - 1]) / B[i];

            C[i] = C[i] / B[i];
            B[i] = 1.0;

        }
    }

    /*
     * Devo eliminare elemento sopra la diagonale principale. Partendo dalla terzultima riga e applicando
     * la formula descritto in seguito si eliminano gli elementi al di sopra della diagonale principale e
     * si fa in modo che tutte le equazioni interne dipendano solo dalla variabile di bordo destra.
     * In questo modo alla fine le righe interne avranno la dipendenza SOLO dallle variabili di bordo
     * e non dalle incognite vicine.
     *
     * La formula è la seguente:
     *
     *                Ri = Ri - w * R_{i+1}   con w = (Ci)/(B_{i+1})
     *
     * La 1° riga viene trattata diversamente in quanto  se viene eseguita la stessa formula
     * si ha che l'elemento B[0] non è più pari ad 1 in quanto subisce l'aggiornamento dell'A[1] che gli sta sotto.
     * Bisogna trattare la 1° riga in modo differente.
    */
    for (int i = N - 3; i >= 0; i--) {
        if (i == 0) {
            w = C[i] / B[i + 1];
            B[i] = B[i] - w * A[i + 1];
            D[i] = (D[i] - w * D[i + 1]) / B[i];
            A[i] = A[i] / B[i];
            C[i] = - (w * C[i + 1]) / B[i];
            B[i] = 1.0;
        }
        else {
            w = C[i] / B[i + 1];
            D[i] = D[i] - w * D[i + 1];
            A[i] = A[i] - w * A[i + 1];
            C[i] = - w * C[i + 1];
        }
    }
}

/*
 * --------------------------------- distribute_input ------------------------------------//
 * Funzione che distribuisce l'input tra i vari processi; in input prende:
 * - *file_name: Nome del file da cui leggere l'input.
 * - *n: Puntatore ad un intero in cui salvare la dimensione del vettore da distribuire tra i vari processi
 * - **v: Puntatore a puntatore da riempire la parte del vettore di input che spetta ad ogni processo
 * - **x: Puntatote a puntatore da allocare per ogni processo, conterrà le soluzioni del sistema (ogni processo lavora su un pezzo differente)
 */
void distribute_input(char *file_name, int *n, double **v, double **x) {
    int id;
    int p;
    MPI_Status status ;

    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    FILE *fp = fopen(file_name, "r");
    if (fp == NULL) {
        if (id == 0) {
            printf("Non è stato possibile aprire il file %s\n", file_name);
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    else {
        /* E' l'ultimo processo a leggere di volta in volta il pezzo dell'input da dover distribuire, questo perchè per
         * come sono state definite le MACRO all'ultimo processo spettere sempre parte intera superiore di n/p e quindi
         * avrà sempre la memoria necessaria per memorizzare il pezzo da inviare.
         */
        if (id == p - 1) {
            if (fscanf(fp, "%d", n) != 1) {
                printf("Errore: lettura dimensione input fallita\n");
                fflush(stdout);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        // Viene effettuata una broadcast per inviare "n" ad ogni processo.
        MPI_Bcast (n, 1, MPI_INT, p-1, MPI_COMM_WORLD);

        // Ogni processo alloca il vettore per salvare i dati che gli verranno spediti dal processo p-1
        int block_size = BLOCK_SIZE(id, p, *n);
        *v = malloc(block_size * sizeof(double));

        if (*v == NULL) {
            printf("Vettore V non allocato correttamente\n");
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Alloco anche il vettore x che mi serve per memorizzare le soluzioni del sistema successivamente
        *x = malloc(block_size * sizeof(double));

        if (*x == NULL) {
            printf("Vettore X non allocato correttamente\n");
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Il processo P-1 legge l'input e lo distribuisce al processo corrispondente.
        if (id == p - 1) {
            for (int i = 0; i < p-1; i++) {
                block_size = BLOCK_SIZE(i, p, *n);
                for (int j = 0; j < block_size; j++) {
                    if (fscanf(fp, "%lf", &(*v)[j]) != 1) {
                        printf("Errore, lettura elementi da file fallita!! \n");
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
                }
                MPI_Send (*v, block_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }

            // In questa sezione di codice legge i dati che gli spettano
            block_size = BLOCK_SIZE(id, p, *n);
            for (int j = 0; j < block_size; j++) {
                if (fscanf(fp, "%lf", &(*v)[j]) != 1) {
                    printf("Errore, lettura elementi da file fallita!! \n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            fclose(fp);
        }
        else {
            MPI_Recv (*v, block_size, MPI_DOUBLE, p-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }
    }
}



/*
 * --------------------------------- check_parallel_thomas ------------------------------------//
 * Funzione utilizzata per testare la correttezza della versione parallela dell'algoritmo di Thomas.
 * Banalmente viene eseguito, sullo stesso input su cui è stato eseguito Thomas Parallelo, l'algoritmo di
 * Thomas sequenziale e si confrontano i risultati ottenuti. (Vengono confrontati i risultati ottenuti da un processo
 * con quelli del blocco corrispondente in sequenziale)
 * In input prende:
 * - *file_nameA: Nome del file da cui leggere il vettore A.
 * - *file_nameB: Nome del file da cui leggere il vettore B.
 * - *file_nameC: Nome del file da cui leggere il vettore C.
 * - *file_nameD: Nome del file da cui leggere il vettore D.
 * - *X_parallel: Puntatore al vettore contente le soluzioni trovate con la versione parallela dell'algoritmo di Thomas
 * - block_size: Dimensione del blocco di dati del processo di cui si vuole verificare la soluzione
 * - start_i: Indica l'indice di partenza nel vettore sequenziale X per confrontarlo con X_parallel
 */

int check_parallel_thomas(char *file_nameA, char *file_nameB, char *file_nameC, char *file_nameD, double *X_parallel, int block_size, int start_i) {

    int n;

    FILE *fpA = fopen(file_nameA, "r");
    FILE *fpB = fopen(file_nameB, "r");
    FILE *fpC = fopen(file_nameC, "r");
    FILE *fpD = fopen(file_nameD, "r");

    if (fpA == NULL || fpB == NULL || fpC == NULL || fpD == NULL) {
        printf("Non è stato possibile aprire i file \n");
        return 1;
    }

    fscanf(fpA, "%d", &n);
    fscanf(fpB, "%d", &n);
    fscanf(fpC, "%d", &n);
    fscanf(fpD, "%d", &n);


    double *A = malloc(n * sizeof(double));
    double *B = malloc(n * sizeof(double));
    double *C = malloc(n * sizeof(double));
    double *D = malloc(n * sizeof(double));
    double *X = malloc(n * sizeof(double));

    if (A == NULL || B == NULL || C == NULL || D == NULL || X == NULL) {
        printf("Vettori non allocati correttamente\n");
        return 1;
    }

    for (int i = 0; i < n; i++) {
        if (fscanf(fpA, "%lf", &A[i]) != 1) {
            printf("Errore, lettura vettore A fallita\n");
        }
        if (fscanf(fpB, "%lf", &B[i]) != 1) {
            printf("Errore, lettura vettore B fallita\n");
        }
        if (fscanf(fpC, "%lf", &C[i]) != 1) {
            printf("Errore, lettura vettore C fallita\n");
        }
        if (fscanf(fpD, "%lf", &D[i]) != 1) {
            printf("Errore, lettura vettore D fallita\n");
        }
    }

    thomas_algorithm(n, A, B, C, D, X);


    for (int i = 0; i < block_size; i++) {
        if (fabs(X[start_i + i] - X_parallel[i]) > 0.001) {
            printf("Caso sequenziale e parallelo non coincidono \n");

            free(A);
            free(B);
            free(C);
            free(D);
            free(X);

            fclose(fpA);
            fclose(fpB);
            fclose(fpC);
            fclose(fpD);

            return 1;
        }
    }

    free(A);
    free(B);
    free(C);
    free(D);
    free(X);

    fclose(fpA);
    fclose(fpB);
    fclose(fpC);
    fclose(fpD);

    return 0;
}



int main(int argc, char *argv[]) {

    int n;
    int id;
    int p;
    double *A;
    double *B;
    double *C;
    double *D;
    double *X;
    int size_reduce_system;
    int block_size;
    double elapsed_time;


    // Inizializzo MPI
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    // La dimensionde del sistema ridotto è 2P in quanto prende solo 2 equazioni per processo
    size_reduce_system = 2 * p;

    // Array utilizzato per memorizzare intera 1° riga e ultima riga di ogni processo (Solo A, C e D. B non lo invio, so già che è 1)
    double global_send_raw[size_reduce_system * 3];

    /* -----------------------------  Passo 1  ----------------------------------------//
     * Il passo 1 consiste nel distribuire l'input tra i vari processi.
     * Il processo p-1 legge di volta in volta il blocco di dati destinati al processo i-esimo
     * e glielo invia con una comunicazione point to point.
     * Diamo la possibilità di leggere il blocco al processo P-1 in quanto per come abbiamo definito le
     * MACRO avrà sempre parte intera superiore di n/p e quindi ha sempre la memoria necessaria per
     * memorizzare il blocco di dati che legge dal file.
     */


    distribute_input(argv[1], &n, &A, &X);
    distribute_input(argv[2], &n, &B, &X);
    distribute_input(argv[3], &n, &C, &X);
    distribute_input(argv[4], &n, &D, &X);

    block_size = BLOCK_SIZE(id, p, n);

    MPI_Barrier(MPI_COMM_WORLD);

    // L'applicazione parallela parte e iniziamo a monitorare il tempo di esecuzione
    elapsed_time = -MPI_Wtime();

    /* -----------------------------  Passo 2  ----------------------------------------//
     * Ogni processo applica thomas_v1 al suo blocco di dati.
     */

    thomas_v1(block_size, A, B, C, D);

    /* -----------------------------------  Passo 3  -------------------------------------------------//
     * Ogni processo manda la sua prima e ultima riga al processo con rango 0 (usando una MPI_Gather)
     */

    block_size = BLOCK_SIZE(id, p, n);
    double send_raw[6];

    // Metto prima e ultima riga in un unico vettore in modo da effettuare una singola Gather.
    send_raw[0] = A[0];
    send_raw[1] = A[block_size - 1];
    send_raw[2] = C[0];
    send_raw[3] = C[block_size - 1];
    send_raw[4] = D[0];
    send_raw[5] = D[block_size - 1];

    MPI_Gather(send_raw, 6, MPI_DOUBLE, global_send_raw, 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* -----------------------------------  Passo 4  -------------------------------------------------//
     * Il processo con rango 0 crea un sistema di dimensione 2P contenente tutte le righe ricevute
     * dagli altri processi (2 righe ricevute per processo) e risolve il nuovo sistema con thomas sequenziale.
     * Le soluzioni ottenute sono i valori delle variabili di bordo dei sistemi locali nei processi.
     */

    double *reduce_X = (double *)malloc(size_reduce_system * sizeof(double));
    if (id == 0) {

      // Array che contengono i valori di A, B, C e D del sistema ridotto.
      double globalA[size_reduce_system];
      double globalB[size_reduce_system];
      double globalC[size_reduce_system];
      double globalD[size_reduce_system];

      // Scompongo le righe ricevute nelle componenti A, B, C e D.
      int count = 0;
      for (int i = 0; i < size_reduce_system; i += 2) {
        globalA[i] = global_send_raw[count];
        count++;
        globalA[i+1] = global_send_raw[count];
        count++;

        globalB[i] = 1.0;
        globalB[i+1] = 1.0;

        globalC[i] = global_send_raw[count];
        count++;
        globalC[i+1] = global_send_raw[count];
        count++;

        globalD[i] = global_send_raw[count];
        count++;
        globalD[i+1] = global_send_raw[count];
        count++;
      }

      thomas_algorithm(size_reduce_system, globalA, globalB, globalC, globalD, reduce_X);
    }

    /* -----------------------------------  Passo 5  -------------------------------------------------//
     * Usando una MPI_Scatter il processo con rango 0 invia le variabili di bordo ai processi corrispondenti.
     * X[0] e X[1] verranno inviati al processo P[0], X[2] e X[3] verranno inviati al processo P[1] ecc...
     */

    MPI_Scatter(reduce_X, 2, MPI_DOUBLE, X, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* -----------------------------------  Passo 6  -------------------------------------------------//
     * Ogni processo, conoscendo le variabili di bordo può calcolarsi le variabili interne usando la formula
     * X[i] = D[i] - A[i] * X[bordo_sinistro] - C[i] * X[bordo_destro];
     */

    double x_temp = X[1];
    X[1] = 0;
    X[block_size - 1] = x_temp;

    for (int i = 1; i < block_size - 1; i++) {
        X[i] = D[i] - A[i] * X[0] - C[i] * X[block_size - 1];
    }

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    /* -----------------------------------  Passo 7  -------------------------------------------------//
     * Stampo il tempo di esecuzione parallelo e sequenziale e effettuo il controllo con check_parallel_thomas
     */

    if (id == 0) {

        printf("Tempo parallelo con %d processori: %lf secondi\n", p, elapsed_time);

        int check = check_parallel_thomas(argv[1],argv[2],argv[3], argv[4], X, block_size, BLOCK_LOW(id, p, n));
        if (check == 0) {
            printf("Thomas_parallelo è corretto \n");
        }
    }

    free(reduce_X);
    free(X);
    free(A);
    free(B);
    free(C);
    free(D);

    MPI_Finalize();

    return 0;
}