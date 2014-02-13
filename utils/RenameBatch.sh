#!/bin/bash
files=()

while read -r -d $'\0'; do
    files+=("$REPLY")
done < <(find $1 -iname "*.2s25motif" -print0)


for (( i = 0 ; i  < ${#files[@]}; i++ ))
do
        #echo ${files[$i]%.*}
        mv ${files[$i]%.*}3sproclog  ${files[$i]%.*}.3sproclog
        mv ${files[$i]%.*}.proclog  ${files[$i]%.*}.2s25proclog
done

