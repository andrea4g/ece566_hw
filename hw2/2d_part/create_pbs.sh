#!/bin/bash

Nelement=( "0" "42" "84" "126" "168" "210" "255" "294" "336");

for i in {1..8}; do
  Kprocs="4";
  a="3";
  for j in {1..7}; do
    if [[ $a -ne 5 ]]
    then
      sed "s/Nelement/${Nelement[i]}/g" pref.pbs > "test_${Nelement[i]}_${Kprocs}.pbs";
      sed -i "s/Kprocs/${Kprocs}/g" "test_${Nelement[i]}_${Kprocs}.pbs"
      Kprocs="$(( $a * $a ))";
    fi
    a="(( $a + 1 ))";
  done
done
