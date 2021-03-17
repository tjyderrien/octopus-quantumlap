#!/bin/bash -l

# This serial script links the WFN files from the MF to the BGW folder.

echo "Linking pseudopotentials and charge densities"
for directory in *-*/; do
  if [[ $directory == 1* ]]; then
    continue
  fi
  echo " Working on directory $directory"
  cd $directory
  mkdir -p restart/gs 2> /dev/null
  cd restart/gs
  ln -sf ../../../1-scf/restart/gs/density.obf .
  ln -sf ../../../1-scf/restart/gs/grid .
  cd ../../..
done
echo

echo "Linking WFNs to BGW directory"
cd ../2-bgw/

echo " Working on directory 1-epsilon"
cd 1-epsilon/
ln -sf ../../1-mf/2.1-epsilon/static/WFN WFN
ln -sf ../../1-mf/2.2-epsilon_q/static/WFN WFNq
cd ../

echo " Working on directory 2-sigma"
cd 2-sigma/
ln -sf ../../1-mf/2.1-epsilon/static/WFN WFN_inner
ln -sf ../../1-mf/2.1-epsilon/static/vxc.dat .
ln -sf ../../1-mf/2.1-epsilon/static/RHO .
cd bandstructure/
ln -sf ../../../1-mf/2.1-epsilon/static/WFN WFN_co
ln -sf ../../../1-mf/3-inteqp/static/WFN WFN_fi
ln -sf ../../../1-mf/3-inteqp/static/WFN WFNq_fi
cd ../
cd ../

if [[ -d 3-kernel ]]; then
echo " Working on directory 3-kernel"
cd 3-kernel/
ln -sf ../../1-mf/2.1-epsilon/static/WFN WFN_co
cd ../
fi

if [[ -d 4-absorption ]]; then
echo " Working on directory 4-absorption"
cd 4-absorption/
ln -sf ../../1-mf/2.1-epsilon/static/WFN WFN_co
ln -sf ../../1-mf/4.1-absorption/static/WFN WFN_fi
ln -sf ../../1-mf/4.2-absorption_q/static/WFN WFNq_fi
cd ../
fi

echo
echo "Done!"
