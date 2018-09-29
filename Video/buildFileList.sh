#!/bin/bash
source $QuantumLaPsources/02-Slab2D/functions.sh
set -e

echo ""
echo "** Building up a file containing all names of subfolders."
echo "   How to use: Use the output listoffiles.txt for the input of plotDensity.py"

LocalPath=$( pwd )

# Maybe construct the list of files which contain density.vtk only. 

list_folders_with_vtk () {
 find -name 'density.vtk'
}

# list_folder_with_vtk 
list_folders_with_vtk | awk 'BEGIN {FS="/"} { print $2 }' > listoffiles.tmp
sort -g listoffiles.tmp > listoffiles.tmp2
awk --assign localpath="${LocalPath}" '{ print localpath"/"$1 }' listoffiles.tmp2 > listoffiles.tmp3
# tail -n+2 listoffiles.tmp3 > listoffiles_long.txt
cp listoffiles.tmp3 listoffiles_long.txt
rm listoffiles.tmp listoffiles.tmp2 listoffiles.tmp3

echo "** Building a SHORT VERSION for pratical purposes..."
head -n1 listoffiles_long.txt  > listoffiles_short.txt
tail -n1 listoffiles_long.txt >> listoffiles_short.txt

echo ""

echo "** Preparation of the TD-general files..."
cd ../td.general

echo "   Building laser 1 profile..."
awk '{ print $1, $2, $3 }' laser > laser.tmp
tail -n+7 laser.tmp > laser.txt
echo ""

echo "   Building laser 2 profile..."
awk '{ print $1, $2, $6 }' laser > laser2.tmp
tail -n+7 laser2.tmp > laser2.txt
echo ""

echo "   Building current XYZ profile..."
awk '{ print $1, $2, $3, $4, $5 }' total_current > current.tmp
tail -n+5 current.tmp > current.txt
echo ""

cd - > /dev/null

echo "   TD-general files ready."
echo ""

cp listoffiles_short.txt listoffiles.txt
