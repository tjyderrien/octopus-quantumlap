#!/bin/bash
for i in $( find | grep 'elf' | grep 'vtk' | grep 'td.0000000' ); do ls -lhtr $i; done >& DatesOfELF.txt
