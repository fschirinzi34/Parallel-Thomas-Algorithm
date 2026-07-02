#!/bin/bash

# Benchmark per valutare la Strong Scalability dell'algoritmo

OUTPUT="Strong_Scalability.txt"
> "$OUTPUT"

# Faccio crescere esponenzialmente il numero di processori utilizzati
for i in {0..3}
 do
   
   proc=$((2**i))
   
   echo "Proc = $proc"

   echo -e "Caso con $proc processori:\n" >> "$OUTPUT"

   for j in {1..10}
     do
     
     	t=$(mpirun -np $proc ./parallel_Thomas_Benchmark A.txt B.txt C.txt D.txt | awk '{print $(NF-1)}')
     	
     	echo -e "$t\n" >> "$OUTPUT"
     
     done
   
   
 done

