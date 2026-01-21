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

a = []
b = []
c = []
d = []

# Nell'ultima riga della matrice "a" deve essere 0

a.append(0)
c.append(random.uniform(-10, 10))
d.append(random.uniform(-10, 10))

# Voglio che la matrice sia diagonalmente dominante, quindi voglio che b[i] >= |a[i]| + |c[i]|
b.append(abs(a[0]) + abs(c[0]) + 1)


for i in range(1, N-1):
    a.append(random.uniform(-10, 10))
    c.append(random.uniform(-10, 10))
    d.append(random.uniform(-10, 10))

    b.append(abs(a[i]) + abs(c[i]) + 1)

# Nell'ultima riga della matrice "c" deve essere 0

a.append(random.uniform(-10, 10))
c.append(0)
d.append(random.uniform(-10, 10))

b.append(abs(a[N-1]) + abs(c[N-1]) + 1)

with open(files[0], "w") as f:
    f.write(str(N) + "\n")
    for i in range(N):
        f.write(str(a[i]) + "\n")

with open(files[1], "w") as f:
    f.write(str(N) + "\n")
    for i in range(N):
        f.write(str(b[i]) + "\n")

with open(files[2], "w") as f:
    f.write(str(N) + "\n")
    for i in range(N):
        f.write(str(c[i]) + "\n")

with open(files[3], "w") as f:
    f.write(str(N) + "\n")
    for i in range(N):
        f.write(str(d[i]) + "\n")


