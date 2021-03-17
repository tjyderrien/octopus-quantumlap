#!/bin/bash
source $QuantumLaPsources/02-Slab2D/functions.sh
set -e

echo ""
echo "** Building up a file containing all names of subfolders."
echo "   How to use: Use the output listoffiles.txt for the input of plotDensity.py"

LocalPath=$( pwd )

# Maybe construct the list of files which contain density.vtk only. 

list_folders_bruteforce () {
 # list the folders
 # for each folder
 for folder in $( list_folders )
 do
   cd $folder
   if [ -e "density.vtk" ]
   then
     if [ -e "geometry.vtk" ]
     then
       echo $folder
#        echo "Files well found in $folder."
     fi
   fi
   cd - >& /dev/null
 done > listoffiles.tmp
 #   check if density.vtk AND geometry.vtk are there. 
 #   if yes, we add the folder name to listoffiles.tmp
}

list_folders_with_vtk () {
 find -name 'geometry.vtk' > listoffiles.tmp
}

list_folders_bruteforce 

cat listoffiles.tmp | awk 'BEGIN {FS="/"} { print $2 }' > listoffiles.tmp1


# Sorting the files
sort -g listoffiles.tmp1 > listoffiles.tmp2
awk --assign localpath="${LocalPath}" '{ print localpath"/"$1 }' listoffiles.tmp2 > listoffiles.tmp3
# tail -n+2 listoffiles.tmp3 > listoffiles_long.txt
cp listoffiles.tmp3 listoffiles_long.txt
rm listoffiles.tmp listoffiles.tmp1 listoffiles.tmp2 listoffiles.tmp3

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
