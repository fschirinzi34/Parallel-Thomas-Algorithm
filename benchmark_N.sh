#!/bin/bash

OUTPUT="Benchmark_N.txt"
> "$OUTPUT"

for i in {2..8}
 do
   
   echo "N = $((10**i))"
   
   echo -e "Caso con N = $((10**i)):\n" >> "$OUTPUT"
   python3 create_tridiagonal_matrix.py $((10**i)) > /dev/null
   
   for j in {1..10}
     do
     
     	t=$(mpirun -np 1 ./parallel_Thomas_Benchmark A.txt B.txt C.txt D.txt | awk '{print $(NF-1)}')
     	
     	echo "$t" >> "$OUTPUT"
     	
     
     done

    
   
 done
