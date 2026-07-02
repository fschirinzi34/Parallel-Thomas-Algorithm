# -------------------------------- create_tridiagonal_matrix.py -----------------------------#
# Script utilizzato per creare una matrice tridiagonale diagonalmente dominante di dimensione n
# Esempio di utilizzo: "python3 create_tridiagonal_matrix.py 10"

import random
import sys

if len(sys.argv) > 1:
    N = int(sys.argv[1])
else:
    print("Errore. Inserire la dimensione della matrice Tridiagonale (N) \n")

files = ['A.txt', 'B.txt', 'C.txt', 'D.txt']

random.seed()

with open(files[0], "w") as fa, \
     open(files[1], "w") as fb, \
     open(files[2], "w") as fc, \
     open(files[3], "w") as fd:

    fa.write(str(N) + "\n")
    fb.write(str(N) + "\n")
    fc.write(str(N) + "\n")
    fd.write(str(N) + "\n")


    for i in range(0, N):
        a = random.uniform(-10, 10)
        c = random.uniform(-10, 10)
        d = random.uniform(-10, 10)

        # Voglio che la matrice sia diagonalmente dominante, quindi voglio che b[i] >= |a[i]| + |c[i]|
        b = abs(a) + abs(c) + 1

        # Il primo valore di "a" deve essere 0
        if (i == 0):
            fa.write(str(0) + "\n")
            fc.write(str(c) + "\n")
            fd.write(str(d) + "\n")
            fb.write(str(b) + "\n")

        # L'ultimo elemento del vettore "c" deve essere 0
        elif(i == N-1):
            fa.write(str(a) + "\n")
            fc.write(str(0) + "\n")
            fd.write(str(d) + "\n")
            fb.write(str(b) + "\n")

        else:
            fa.write(str(a) + "\n")
            fc.write(str(c) + "\n")
            fd.write(str(d) + "\n")
            fb.write(str(b) + "\n")


print("File creati con successo \n")
