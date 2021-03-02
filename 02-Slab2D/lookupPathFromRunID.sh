#! /bin/bash

echo "From run ID, identify the path."
echo "Usage: source <script>"
echo "Declares functions to be used in the following sequence: "
echo "1. FetchPathFromID ( filename )"
echo "2. PlotEachRun ( dirnames.tmp  ) "
echo "3. Eventually disableProblematicRuns.sh <file>"

# Returns folder name from runID
FetchPathFromID() {
for line in $( cat $1 )
do
  PathOfTheRun=$( dirname $( export depth=4; find . -mindepth $depth -maxdepth $depth | grep $line | tail -n1 ) )
  echo $PathOfTheRun
done > dirnames.tmp && echo "List was generated in dirnames.tmp."
}

# Prepare a plot for each subfolder indicated into input file
PlotEachRun() {
counter=0
for line in $( cat $1 )
do
  counter=$( python -c "print( $counter + 1 )" )
  gnuplot -p -e "folder='$line'.'/td.general/'; set title folder; set xrange [:4000]; p folder.'laser' u 2:3 w l t 'Pulse 1', folder.'laser' u 2:6 w l t 'Pulse 2', folder.'n_ex' u 2:3 w l t 'Nexc'; print folder"
  echo $counter / $( wc -l $1 )
  read -s -n 1 -p "Press any key to continue . . ."
done
}
# FetchPathFromID ListOfCandidates.dat
# PlotEachRun dirnames.tmp
# /home/thibault.derrien/Documents/octopus_slabs/02-Slab2D/disableProblematicRuns.sh ListOfCandidates.dat
