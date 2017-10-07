#!/bin/bash


filename="$(ls *.out)";

cat /dev/null > "out.txt"
for file in $filename; do
  echo $file >> "out.txt";
  cat $file >> "out.txt";
done
