#!/bin/bash

set -e

echo ""
echo "** Building up a file containing all names of subfolders."
echo "   How to use: Use the output listoffiles.txt for the input of plotDensity.py"

LocalPath=$( pwd )

du | awk '{ print $2 }' | awk 'BEGIN {FS="./"} { print $2 }' > listoffiles.tmp #TODO: THIS IS TOTALLY INEFFICIENT. 
sort -g listoffiles.tmp > listoffiles.tmp2
# cat listoffiles.tmp2 > listoffiles.tmp3
awk --assign localpath="${LocalPath}" '{ print localpath"/"$1 }' listoffiles.tmp2 > listoffiles.tmp3
tail -n+2 listoffiles.tmp3 > listoffiles.txt
rm listoffiles.tmp listoffiles.tmp2 listoffiles.tmp3

echo ""

echo "** Preparation of the TD-general files..."
cd ../td.general

echo "   Building laser 1 profile..."
awk '{ print $1, $2, $3 }' laser > laser.tmp
tail -n+7 laser.tmp > laser.txt
echo ""

echo "   Building laser 2 profile..."
awk '{ print $1, $2, $6 }' laser > laser.tmp
tail -n+7 laser2.tmp > laser2.txt
echo ""

echo "   Building current XYZ profile..."
awk '{ print $1, $2, $3, $4, $5 }' total_current > current.tmp
tail -n+5 current.tmp > current.txt
echo ""


cd - > /dev/null

echo "   TD-general files ready."
echo ""
