# L'obiettivo di questo script è creare una matrice tridiagonale diagonalmente dominante.
import random

N = 100000000  # numero di elementi
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


