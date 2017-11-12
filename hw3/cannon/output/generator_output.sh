#!/bin/bash


filename="$(ls *.out)";

cat /dev/null > "out.txt"
for file in $filename; do
  echo -n "$file, " >> "out.txt";
  cat $file | head -1 >> "out.txt";
done
