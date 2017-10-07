#!/bin/bash


Nelement="1024";

for i in {1..7}; do
  sed "s/Nelement/${Nelement}/g" pref.pbs > "test_$Nelement.pbs";
  sed -i "s/number_of_nodes/1/g" "test_$Nelement.pbs"
  sed -i "s/procs_per_node/1/g" "test_$Nelement.pbs"
  sed -i "s/Kprocs/16/g" "test_$Nelement.pbs"
  Nelement="$(( $Nelement * 4))";
done

