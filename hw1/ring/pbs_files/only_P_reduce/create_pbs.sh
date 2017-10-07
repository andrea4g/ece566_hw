#!/bin/bash


Nelement="268435456";

for i in {1..7}; do
  kp=$(( $i*16 ));
  echo $kp
  sed "s/Nelement/${Nelement}/g" pref.pbs > "test_$kp.pbs";
  sed -i "s/number_of_nodes/$i/g" "test_$kp.pbs"
  sed -i "s/procs_per_node/16/g" "test_$kp.pbs"
  sed -i "s/Kprocs/$kp/g" "test_$kp.pbs"
done

