#!/usr/bin/env python
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
dir1 = '../results/aerosol-platt-nh3-isopa/aero/'
dir2 = '../results/aerosol-platt-nh3-isopa-oligo/aero/'

Nt= 5 * 60
#####################################################################


font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

shape = [Nt]
inorg=np.zeros(shape)
isop1=np.zeros(shape)
isop2=np.zeros(shape)
oligo2=np.zeros(shape)

with open (dir1+'/PBiPER_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop1[j] = isop1[j] + float(string_conc1[j])
with open (dir1+'/PBiDER_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop1[j] = isop1[j] + float(string_conc1[j])
with open (dir1+'/PBiMT_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop1[j] = isop1[j] + float(string_conc1[j])
with open (dir1+'/PBiNIT3_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop1[j] = isop1[j] + float(string_conc1[j])
with open (dir1+'/PBiMGA_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop1[j] = isop1[j] + float(string_conc1[j])
with open (dir1+'/PBiNGA_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop1[j] = isop1[j] + float(string_conc1[j])
              
                      

with open (dir2+'/PBiPER_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop2[j] = isop2[j] + float(string_conc1[j])
with open (dir2+'/PBiDER_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop2[j] = isop2[j] + float(string_conc1[j])
with open (dir2+'/PBiMT_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop2[j] = isop2[j] + float(string_conc1[j])
with open (dir2+'/PBiNIT3_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop2[j] = isop2[j] + float(string_conc1[j])
with open (dir2+'/PBiMGA_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop2[j] = isop2[j] + float(string_conc1[j])
with open (dir2+'/PBiNGA_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop2[j] = isop2[j] + float(string_conc1[j])

           
with open (dir2+'/POligoBiPER_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   oligo2[j] = oligo2[j] + float(string_conc1[j])
with open (dir2+'/POligoBiDER_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   oligo2[j] = oligo2[j] + float(string_conc1[j])
with open (dir2+'/POligoBiMT_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   oligo2[j] = oligo2[j] + float(string_conc1[j])
with open (dir2+'/POligoBiNIT3_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   oligo2[j] = oligo2[j] + float(string_conc1[j])
with open (dir2+'/POligoBiMGA_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   oligo2[j] = oligo2[j] + float(string_conc1[j])
with open (dir2+'/POligoBiNGA_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   oligo2[j] = oligo2[j] + float(string_conc1[j])

time=range(60,60*(Nt+1),60)
time=numpy.array(time)
time=time/3600.0

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
fig = plt.figure(1,figsize = (15,15))
plt.ylim(0.,150.)

plt.plot(time, isop1,'k',label='Without oligomerization')
plt.plot(time, isop2+oligo2,'b',label='With oligomerization')
plt.plot(time, isop2,'b-.',label='Monomers')
plt.plot(time, oligo2,'b--',label='Oligomers')
plt.ylabel("Isoprene SOA ($\mu g.m^{-3}$)", fontsize=18)
plt.xlabel("time [h]", fontsize=18)
plt.legend(loc=2)
plt.savefig("platt-oligo.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)
