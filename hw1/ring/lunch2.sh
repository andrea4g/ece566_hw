#!/bin/bash

filename="$(ls *.pbs)"

for file in $filename; do
  qsub $file
done

