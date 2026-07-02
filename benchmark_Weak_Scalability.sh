#!/bin/bash

# Benchmark per valutare la Weak Scalability dell'algoritmo

# File output
OUTPUT="Weak_Scalability.txt"
> "$OUTPUT"

N=10000000

# Faccio crescere esponenzialmente il numero di processori utilizzati
for i in {0..3}
 do

   proc=$((2**i))
   new_N=$((N * proc))
   
   echo "Proc = $proc"

   python3 create_tridiagonal_matrix.py $new_N >> /dev/null

   echo -e "Caso con $proc processori e N = $new_N:\n" >> "$OUTPUT"

   for j in {1..10}
     do

     	t=$(mpirun -np $proc ./parallel_Thomas_Benchmark A.txt B.txt C.txt D.txt | awk '{print $(NF-1)}')

     	echo -e "$t\n" >> "$OUTPUT"

     done

  
 done

