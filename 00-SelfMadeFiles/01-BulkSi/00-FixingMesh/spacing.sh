#!/bin/bash
list="0.40 0.35 0.30 0.25"
export OCT_PARSE_ENV=1
for OCT_Spacing in $list
do
  export OCT_Spacing
  octopus_mpi >& out-$OCT_Spacing
  energy=`grep Total static/info`
  echo $OCT_Spacing $energy
  rm -rf restart
done
unset OCT_Spacing
