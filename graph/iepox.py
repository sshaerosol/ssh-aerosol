#!/usr/bin/env python3
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
dir1 = '../results/aerosol-iepox/aero/'
dir2 = '../results/aerosol-iepox1/aero/'
dir3 = '../results/aerosol-iepox2/aero/'
#dir4 = '../results/aerosol-platt-nh3-isopa-oligo-inorg/aero/'


Nt= 24 * 60
#####################################################################


font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

shape = [Nt]
org1=np.zeros(shape)
iepox1=np.zeros(shape)
mt1=np.zeros(shape)
sulf1=np.zeros(shape)
nit1=np.zeros(shape)

org2=np.zeros(shape)
iepox2=np.zeros(shape)
mt2=np.zeros(shape)
sulf2=np.zeros(shape)
nit2=np.zeros(shape)

org3=np.zeros(shape)
iepox3=np.zeros(shape)
mt3=np.zeros(shape)
sulf3=np.zeros(shape)
nit3=np.zeros(shape)

with open (dir1+'/Organic.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  org1[j] = org1[j] + float(string_conc1[j])

with open (dir1+'/../gas/IEPOX.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  iepox1[j] = iepox1[j] + float(string_conc1[j])

with open (dir1+'/PBiMT_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  mt1[j] = mt1[j] + float(string_conc1[j])

with open (dir1+'/PBiNO3_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  nit1[j] = nit1[j] + float(string_conc1[j])

with open (dir1+'/PBiSULF_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  sulf1[j] = sulf1[j] + float(string_conc1[j])



with open (dir2+'/Organic.txt') as finit :
        string_conc2 = finit.read().splitlines()
        for j in range(Nt) :
                  org2[j] = org2[j] + float(string_conc2[j])

with open (dir2+'/../gas/IEPOX.txt') as finit :
        string_conc2 = finit.read().splitlines()
        for j in range(Nt) :
                  iepox2[j] = iepox2[j] + float(string_conc2[j])

with open (dir2+'/PBiMT_1.txt') as finit :
        string_conc2 = finit.read().splitlines()
        for j in range(Nt) :
                  mt2[j] = mt2[j] + float(string_conc2[j])

with open (dir2+'/PBiNO3_1.txt') as finit :
        string_conc2 = finit.read().splitlines()
        for j in range(Nt) :
                  nit2[j] = nit2[j] + float(string_conc2[j])

with open (dir2+'/PBiSULF_1.txt') as finit :
        string_conc2 = finit.read().splitlines()
        for j in range(Nt) :
                  sulf2[j] = sulf2[j] + float(string_conc2[j])

      
with open (dir3+'/Organic.txt') as finit :
        string_conc3 = finit.read().splitlines()
        for j in range(Nt) :
                  org3[j] = org3[j] + float(string_conc3[j])

with open (dir3+'/../gas/IEPOX.txt') as finit :
        string_conc3 = finit.read().splitlines()
        for j in range(Nt) :
                  iepox3[j] = iepox3[j] + float(string_conc3[j])

with open (dir3+'/PBiMT_1.txt') as finit :
        string_conc3 = finit.read().splitlines()
        for j in range(Nt) :
                  mt3[j] = mt3[j] + float(string_conc3[j])

with open (dir3+'/PBiNO3_1.txt') as finit :
        string_conc3 = finit.read().splitlines()
        for j in range(Nt) :
                  nit3[j] = nit3[j] + float(string_conc3[j])

with open (dir3+'/PBiSULF_1.txt') as finit :
        string_conc3 = finit.read().splitlines()
        for j in range(Nt) :
                  sulf3[j] = sulf3[j] + float(string_conc3[j])            
                  

time=range(60,60*(Nt+1),60)
time=numpy.array(time)
time=time/3600.0

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
fig = plt.figure(1,figsize = (15,15))
plt.ylim(ymin=0.)
plt.xlim(xmin=0.)
plt.xlim(xmax=time[Nt-1])
plt.ylim(ymax=17.)

plt.plot(time, org1,'k',label='SOA (no NH$_3$)')
#plt.plot(time, iepox1,'0.5',label='gas-phase IEPOX')
plt.plot(time, mt1,'g',label='Particulate methyltetrols')
plt.plot(time, sulf1,'r',label='Particulate organosulfate')
#plt.plot(time, nit1,'b',label='Particulate organonitrate')

plt.plot(time, org2,'--k',label='SOA for NH$_3$/SO$_4$=1')
plt.plot(time, mt2,'--g')
plt.plot(time, sulf2,'--r')

plt.plot(time, org3,'-.k',label='SOA for NH$_3$/SO$_4$=2')
plt.plot(time, mt3,'-.g')
plt.plot(time, sulf3,'-.r')

#plt.plot(time, nit1,'--b')

#plt.plot(time, org3,'-.k',label='SOA for NH3=3ug/m3')

plt.ylabel("Concentrations ($\mu g.m^{-3}$)", fontsize=18)
plt.xlabel("time [h]", fontsize=18)
plt.legend(loc=2)
plt.savefig("iepox.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)
