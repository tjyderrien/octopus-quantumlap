#!/bin/bash
sudo apt-get install qiv
for x in $(ls *.eps); do echo $x; convert $x ${x/eps/png}; done
qiv 800.0nm-Phase0.0+400.0nm-Phase*png
