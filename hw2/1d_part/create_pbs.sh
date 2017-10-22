#!/bin/bash

Nelement="60";

for i in {1..7}; do
  Kprocs="2";
  for j in {1..10}; do
    sed "s/Nelement/${Nelement}/g" pref.pbs > "test_${Nelement}_${Kprocs}.pbs";
    sed -i "s/Kprocs/${Kprocs}/g" "test_${Nelement}_${Kprocs}.pbs"
    Kprocs="$(( $Kprocs + $Kprocs/2 ))";
  done
  Nelement="$(( $Nelement + $Nelement/3 ))";
done

