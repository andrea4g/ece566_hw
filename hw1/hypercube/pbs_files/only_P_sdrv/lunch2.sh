#!/bin/bash

filename="$(ls test_*)"

for file in $filename; do
  qsub $file
done

