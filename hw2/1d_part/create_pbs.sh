#!/bin/bash

Nelement="49";

for i in {1..7}; do
  sed "s/Nelement/${Nelement}/g" pref.pbs > "test_$Nelement.pbs";
  sed -i "s/number_of_nodes/7/g" "test_$Nelement.pbs"
  sed -i "s/procs_per_node/7/g" "test_$Nelement.pbs"
  sed -i "s/Kprocs/49/g" "test_$Nelement.pbs"
  Nelement="$(( $Nelement + $Nelement/2 ))";
done

