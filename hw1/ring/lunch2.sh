#!/bin/bash

filename="$(ls test_*)"

for file in $filename; do
  echo $file
done

