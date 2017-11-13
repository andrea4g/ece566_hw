#!/bin/bash

Nelement="64";

for i in {1..6}; do
  Kprocs="4";
  for j in {1..3}; do
    for q in {1..5}; do
      sed "s/Nelement/${Nelement}/g" pref.pbs > "test_${Nelement}_${Kprocs}_${q}.pbs";
      sed -i "s/Kprocs/${Kprocs}/g" "test_${Nelement}_${Kprocs}_${q}.pbs"
      sed -i "s/Exp/${q}/g" "test_${Nelement}_${Kprocs}_${q}.pbs"
    done
    Kprocs="$(( $Kprocs * 4 ))";
  done
  Nelement="$(( $Nelement * 4 ))";
done

