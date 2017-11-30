#! /bin/bash

i=0
while read line
do
  filename[ $i ]="$line"
  n[ $i ]="${filename[ $i ]//[!0-9]/}"
  echo ${n[i]}
  (( i++ ))
done < <(ls atsp/)

