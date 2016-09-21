#!/usr/bin/env bash

# This serial script links the epsmat and bsemat files.

echo "Linking epsmat files"
for directory in *-*/; do
  # Ignore the "1-epsilon/" directory and things that are not directories
  if [[ ( $directory == 1-* ) || ! -d $directory ]]; then continue; fi
  echo " Working on directory $directory"

  cd $directory
  ln -sf ../1-epsilon/epsmat.h5 ./
  ln -sf ../1-epsilon/eps0mat.h5 ./
  cd ../
done

if [[ -d 4-absorption ]]; then
cd 4-absorption/
ln -s ../3-kernel/bsemat.h5 ../4-absorption
cd ..
fi

echo
echo "Done!"
