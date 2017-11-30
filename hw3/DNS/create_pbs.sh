#!/bin/bash

for a in {2..5}; do
  for q in {2..4}; do
    Nelement=240
    for i in {1..4}; do
      Kprocs="$(($a*$a*$a))";
      sed "s/Nelement/${Nelement}/g" pref.pbs > "test_${Nelement}_${Kprocs}_${q}.pbs";
      sed -i "s/Kprocs/${Kprocs}/g" "test_${Nelement}_${Kprocs}_${q}.pbs"
      sed -i "s/Exp/${q}/g" "test_${Nelement}_${Kprocs}_${q}.pbs"
      Nelement="$(( $Nelement * 2 ))";
    done
  done
done
