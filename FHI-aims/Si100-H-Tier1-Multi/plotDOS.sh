#/bin/bash

echo "Plotting DOS for each directory..."

for i in `du | awk '{ print $2 }' | grep 'k-'`
do
	cd $i
	if [ -e KS_DOS_total.dat ]; then
		../../../aimsplot.py
	fi
	cd -
done
