#!/bin/bash

for a in {2..5}; do
  for q in {2..4}; do 
    for i in {1..4}; do
      Nelement="$(( $a * $a ))";
      Kprocs="$(($a*$a*$a))";
      for j in {1..4}; do
        if [[ $j -gt 1 ]]
        then
          Nelement="$(( $Nelement * 10 ))";
        fi
        sed "s/Nelement/${Nelement}/g" pref.pbs > "test_${Nelement}_${Kprocs}_${q}.pbs";
        sed -i "s/Kprocs/${Kprocs}/g" "test_${Nelement}_${Kprocs}_${q}.pbs"
        sed -i "s/Exp/${q}/g" "test_${Nelement}_${Kprocs}_${q}.pbs"
      done
    done
  done
done
