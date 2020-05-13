#!/bin/bash

echo "Usage: shows the number of subfolders in each folder located directly below the current one."
echo "Info: This is useful to know how many Octopus runs were generated in a subfolder."
for folder in $( find . -maxdepth 1 -type d )
do
  cd $folder 
  NumOfCases=$(find . -maxdepth 2 -mindepth 2 -type d | wc -l)
  echo -e "$NumOfCases\t$folder"
  cd - > /dev/null
done | sort -g
