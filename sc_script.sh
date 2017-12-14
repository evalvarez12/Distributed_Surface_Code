#!/bin/bash

distance="16"
topology="toric"

# Values
declare -a arr=("0.09" "0.092" "0.094" "0.096" "0.098" "0.10" "0.102"
                "0.104" "0.106" "0.108" "0.110")

for p in ${arr[@]}; do
  cmd="mpiexec -n 4 python one_point.py distance=$distance topology=$topology p=$p q=0 iterations=25000"
  echo $cmd
done
