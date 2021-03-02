#!/usr/bin/env python

# This script generates BN_bandstructure.pdf from the output of inteqp
# using numpy and matplotlib. Parameters are hard-coded into the script.

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt

f = open('bandstructure.dat', 'r')
f_inp = f.readlines()
f.close()

nhead = 2
ntot = len(f_inp) - nhead
bndmin = int(f_inp[nhead].split()[1]) - 1
bndmax = int(f_inp[nhead + ntot - 1].split()[1]) - 1
nbnd = bndmax - bndmin + 1
nkpt = int(ntot / nbnd)
kpoint = np.zeros((nkpt))
eqp  = np.zeros((nbnd, nkpt))
elda = np.zeros((nbnd, nkpt))
klen = 0.0
for ik in range(nkpt):
   ktext = f_inp[nhead + ik][14:50].split()
   knew = [float(ktext[0]), float(ktext[1]), float(ktext[2])]
   if ik > 0:
      klen += np.sqrt( (knew[0]-kold[0])**2 + (knew[1]-kold[1])**2 + (knew[2]-kold[2])**2 )
   kold = knew
   kpoint[ik] = klen
   for ib in range(nbnd):
      eqp[ib][ik]  = float(f_inp[nhead + ik + ib * nkpt].split()[6])
      elda[ib][ik] = float(f_inp[nhead + ik + ib * nkpt].split()[5])

efermi = -3.0
emin = efermi - 12.0
emax = efermi + 15.0
kmin = kpoint[0]
kmax = kpoint[nkpt - 1]

matplotlib.rc('figure', figsize = (8.0, 6.0))
matplotlib.rc('axes', linewidth = 1.5)
matplotlib.rc('lines', linewidth = 1.5)
matplotlib.rc('font', size = 18.0)
matplotlib.rc('xtick.major', size = 0.0, pad = 8.0)
matplotlib.rc('xtick.minor', size = 0.0, pad = 8.0)
matplotlib.rc('ytick.major', size = 6.0, pad = 8.0)
matplotlib.rc('ytick.minor', size = 3.0, pad = 8.0)

kindex = [0, 4, 9]
klabel = ['K', 'X', '$\Gamma$']
kvalue = []
for i in range(len(kindex)):
   kvalue.append(kpoint[kindex[i]])
myxticks = plt.xticks(kvalue, klabel)

for ib in range(nbnd):
   myplot = plt.plot(kpoint, eqp[ib], 'b-')

for ib in range(nbnd):
   myplot = plt.plot(kpoint, elda[ib], 'g-')

myplot = plt.plot(kpoint, eqp[0], 'b-', label="GW")
myplot = plt.plot(kpoint, elda[0], 'g-', label="LDA")

for i in range(len(kindex)):
   if kvalue[i] != kmin and kvalue[i] != kmax:
      myplot = plt.plot([kvalue[i], kvalue[i]], [emin, emax], 'k-')
myplot = plt.plot([kmin, kmax], [efermi, efermi], 'k-')

mygca = plt.gca()
for i in mygca.get_xticklines() + mygca.get_yticklines():
   i.set_markeredgewidth(1.5)

myaxis = plt.axis([kmin, kmax, emin, emax])
myxlabel = plt.xlabel('Wavevector')
myylabel = plt.ylabel('Energy (eV)')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=2, mode="expand", borderaxespad=0.)

plt.savefig('BN_bands.pdf', bbox_inches='tight')
