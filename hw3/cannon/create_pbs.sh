#!/bin/bash

Nelement="42";
mult="2";

for i in {1..5}; do
  Kprocs="4";
  a="3";
  for j in {1..7}; do
    if [[ $a -ne 5 ]]
    then
      for q in {2..4}; do
        sed "s/Nelement/${Nelement}/g" pref.pbs > "test_${Nelement}_${Kprocs}_${q}.pbs";
        sed -i "s/Kprocs/${Kprocs}/g" "test_${Nelement}_${Kprocs}_${q}.pbs"
        sed -i "s/Exp/${q}/g" "test_${Nelement}_${Kprocs}_${q}.pbs"
      done
      Kprocs="$(( $a * $a ))";
    fi
    a="(( $a + 1 ))";
  done
  Nelement="$(( $Nelement * 2 ))";
done
