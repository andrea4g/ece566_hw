#!/bin/bash

if [  -f ./mpitest.out ]; then
  rm ./mpitest.out
fi

cat /dev/null > out.txt

for i in {1..20}; do
  qsub mpitest.pbs
  while [ ! -f ./mpitest.out ]; do
    sleep 1
  done
  cat mpitest.out | cut -d " " -f4 | sed '/^\s*$/d' | sort -g | tail -1 >> out.txt 
  rm mpitest.out
done

qsub mpitest.pbs
