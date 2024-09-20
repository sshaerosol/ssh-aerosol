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
dir4 = '../results/aerosol-iepox-iso/aero/'
dir5 = '../results/aerosol-iepox1-iso/aero/'
dir6 = '../results/aerosol-iepox2-iso/aero/'
dir7 = '../results/aerosol-iepox-dyn/aero/'
dir8 = '../results/aerosol-iepox-iso-dyn/aero/'
#dir4 = '../results/aerosol-platt-nh3-iso
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

org4=np.zeros(shape)
iepox4=np.zeros(shape)
mt4=np.zeros(shape)
sulf4=np.zeros(shape)
nit4=np.zeros(shape)

org5=np.zeros(shape)
iepox5=np.zeros(shape)
mt5=np.zeros(shape)
sulf5=np.zeros(shape)
nit5=np.zeros(shape)
org6=np.zeros(shape)
iepox6=np.zeros(shape)
mt6=np.zeros(shape)
sulf6=np.zeros(shape)
nit6=np.zeros(shape)

org7=np.zeros(shape)
iepox7=np.zeros(shape)
mt7=np.zeros(shape)
sulf7=np.zeros(shape)
nit7=np.zeros(shape)

org8=np.zeros(shape)
iepox8=np.zeros(shape)
mt8=np.zeros(shape)
sulf8=np.zeros(shape)
nit8=np.zeros(shape)



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

with open (dir4+'/Organic.txt') as finit :
        string_conc4 = finit.read().splitlines()
        for j in range(Nt) :
                  org4[j] = org4[j] + float(string_conc4[j])

with open (dir4+'/../gas/IEPOX.txt') as finit :
        string_conc4 = finit.read().splitlines()
        for j in range(Nt) :
                  iepox4[j] = iepox4[j] + float(string_conc4[j])

with open (dir4+'/PBiMT_1.txt') as finit :
        string_conc4 = finit.read().splitlines()
        for j in range(Nt) :
                  mt4[j] = mt4[j] + float(string_conc4[j])

with open (dir4+'/PBiNO3_1.txt') as finit :
        string_conc4 = finit.read().splitlines()
        for j in range(Nt) :
                  nit4[j] = nit4[j] + float(string_conc4[j])

with open (dir4+'/PBiSULF_1.txt') as finit :
        string_conc4 = finit.read().splitlines()
        for j in range(Nt) :
                  sulf4[j] = sulf4[j] + float(string_conc4[j])            
                  
with open (dir5+'/Organic.txt') as finit :
        string_conc5 = finit.read().splitlines()
        for j in range(Nt) :
                  org5[j] = org5[j] + float(string_conc5[j])

with open (dir5+'/../gas/IEPOX.txt') as finit :
        string_conc5 = finit.read().splitlines()
        for j in range(Nt) :
                  iepox5[j] = iepox5[j] + float(string_conc5[j])

with open (dir5+'/PBiMT_1.txt') as finit :
        string_conc5 = finit.read().splitlines()
        for j in range(Nt) :
                  mt5[j] = mt5[j] + float(string_conc5[j])

with open (dir5+'/PBiNO3_1.txt') as finit :
        string_conc5 = finit.read().splitlines()
        for j in range(Nt) :
                  nit5[j] = nit5[j] + float(string_conc5[j])

with open (dir5+'/PBiSULF_1.txt') as finit :
        string_conc5 = finit.read().splitlines()
        for j in range(Nt) :
                  sulf5[j] = sulf5[j] + float(string_conc5[j])            


with open (dir6+'/Organic.txt') as finit :
        string_conc6 = finit.read().splitlines()
        for j in range(Nt) :
                  org6[j] = org6[j] + float(string_conc6[j])

with open (dir6+'/../gas/IEPOX.txt') as finit :
        string_conc6 = finit.read().splitlines()
        for j in range(Nt) :
                  iepox6[j] = iepox6[j] + float(string_conc6[j])

with open (dir6+'/PBiMT_1.txt') as finit :
        string_conc6 = finit.read().splitlines()
        for j in range(Nt) :
                  mt6[j] = mt6[j] + float(string_conc6[j])

with open (dir6+'/PBiNO3_1.txt') as finit :
        string_conc6 = finit.read().splitlines()
        for j in range(Nt) :
                  nit6[j] = nit6[j] + float(string_conc6[j])

with open (dir6+'/PBiSULF_1.txt') as finit :
        string_conc6 = finit.read().splitlines()
        for j in range(Nt) :
                  sulf6[j] = sulf6[j] + float(string_conc6[j])            

with open (dir7+'/Organic.txt') as finit :
        string_conc7 = finit.read().splitlines()
        for j in range(Nt) :
                  org7[j] = org7[j] + float(string_conc7[j])

with open (dir7+'/../gas/IEPOX.txt') as finit :
        string_conc7 = finit.read().splitlines()
        for j in range(Nt) :
                  iepox7[j] = iepox7[j] + float(string_conc7[j])

with open (dir7+'/PBiMT_1.txt') as finit :
        string_conc7 = finit.read().splitlines()
        for j in range(Nt) :
                  mt7[j] = mt7[j] + float(string_conc7[j])

with open (dir7+'/PBiNO3_1.txt') as finit :
        string_conc7 = finit.read().splitlines()
        for j in range(Nt) :
                  nit7[j] = nit7[j] + float(string_conc7[j])

with open (dir7+'/PBiSULF_1.txt') as finit :
        string_conc7 = finit.read().splitlines()
        for j in range(Nt) :
                  sulf7[j] = sulf7[j] + float(string_conc7[j])

with open (dir8+'/Organic.txt') as finit :
        string_conc8 = finit.read().splitlines()
        for j in range(Nt) :
                  org8[j] = org8[j] + float(string_conc8[j])

with open (dir8+'/../gas/IEPOX.txt') as finit :
        string_conc8 = finit.read().splitlines()
        for j in range(Nt) :
                  iepox8[j] = iepox8[j] + float(string_conc8[j])

with open (dir8+'/PBiMT_1.txt') as finit :
        string_conc8 = finit.read().splitlines()
        for j in range(Nt) :
                  mt8[j] = mt8[j] + float(string_conc8[j])

with open (dir8+'/PBiNO3_1.txt') as finit :
        string_conc8 = finit.read().splitlines()
        for j in range(Nt) :
                  nit8[j] = nit8[j] + float(string_conc8[j])

with open (dir8+'/PBiSULF_1.txt') as finit :
        string_conc8 = finit.read().splitlines()
        for j in range(Nt) :
                  sulf8[j] = sulf8[j] + float(string_conc8[j])     

time=range(60,60*(Nt+1),60)
time=numpy.array(time)
time=time/3600.0

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#fig = plt.figure(1,figsize = (15,15))
plt.ylim(ymin=0.)
plt.xlim(xmin=0.)
plt.xlim(xmax=time[Nt-1])
plt.ylim(ymax=17.)

plt.plot(time, org1,'k',label='SOA (no NH$_3$)')
plt.plot(time, org2,'--k',label='SOA for NH$_3$/SO$_4$=1')
plt.plot(time, org3,'-.k',label='SOA for NH$_3$/SO$_4$=2')

plt.plot(time, org4,'r',label='SOA (no NH$_3$) with Isorropia')
plt.plot(time, org5,'--r',label='SOA for NH$_3$/SO$_4$=1 with Isorropia')
plt.plot(time, org6,'-.r',label='SOA for NH$_3$/SO$_4$=2 with Isorropia')

plt.ylabel("Concentrations ($\mu g.m^{-3}$)", fontsize=14)
plt.xlabel("time [h]", fontsize=14)
plt.title("SOA mass")
plt.legend(loc=2,fontsize=10,ncol=2)
plt.savefig("iepox.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#fig = plt.figure(1,figsize = (15,15))
plt.ylim(ymin=0.)
plt.xlim(xmin=0.)
plt.xlim(xmax=time[Nt-1])
plt.ylim(ymax=17.)

plt.plot(time, sulf1,'k',label='SOA (no NH$_3$)')
plt.plot(time, sulf2,'--k',label='SOA for NH$_3$/SO$_4$=1')
plt.plot(time, sulf3,'-.k',label='SOA for NH$_3$/SO$_4$=2')

plt.plot(time, sulf4,'r',label='SOA (no NH$_3$) with Isorropia')
plt.plot(time, sulf5,'--r',label='SOA for NH$_3$/SO$_4$=1 with Isorropia')
plt.plot(time, sulf6,'-.r',label='SOA for NH$_3$/SO$_4$=2 with Isorropia')

plt.ylabel("Concentrations ($\mu g.m^{-3}$)", fontsize=14)
plt.xlabel("time [h]", fontsize=14)
plt.title("Organosulfate mass")
plt.legend(loc=2,fontsize=10,ncol=2)
plt.savefig("iepox-bisulf.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)

#fig = plt.figure()
##ax = fig.add_subplot(1,1,1)
##fig = plt.figure(1,figsize = (15,15))
#plt.ylim(ymin=0.)
#plt.xlim(xmin=0.)
#plt.xlim(xmax=time[Nt-1])
#plt.ylim(ymax=17.)

#plt.plot(time, nit1,'k',label='SOA (no NH$_3$)')
#plt.plot(time, nit2,'--k',label='SOA for NH$_3$/SO$_4$=1')
#plt.plot(time, nit3,'-.k',label='SOA for NH$_3$/SO$_4$=2')

#plt.plot(time, nit4,'r',label='SOA (no NH$_3$) with Isorropia')
#plt.plot(time, nit5,'--r',label='SOA for NH$_3$/SO$_4$=1 with Isorropia')
#plt.plot(time, nit6,'-.r',label='SOA for NH$_3$/SO$_4$=2 with Isorropia')

#plt.ylabel("Concentrations ($\mu g.m^{-3}$)", fontsize=14)
#plt.xlabel("time [h]", fontsize=14)
#plt.title("Organonitrate mass")
#plt.legend(loc=2,ncol=2,fontsize=10)
#plt.savefig("iepox-binit.png",bbox_inches='tight')
#plt.tight_layout()
#plt.close(fig)

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#fig = plt.figure(1,figsize = (15,15))
plt.ylim(ymin=0.)
plt.xlim(xmin=0.)
plt.xlim(xmax=time[Nt-1])
plt.ylim(ymax=17.)

plt.plot(time, mt1,'k',label='SOA (no NH$_3$)')
plt.plot(time, mt2,'--k',label='SOA for NH$_3$/SO$_4$=1')
plt.plot(time, mt3,'-.k',label='SOA for NH$_3$/SO$_4$=2')

plt.plot(time, mt4,'r',label='SOA (no NH$_3$) with Isorropia')
plt.plot(time, mt5,'--r',label='SOA for NH$_3$/SO$_4$=1 with Isorropia')
plt.plot(time, mt6,'-.r',label='SOA for NH$_3$/SO$_4$=2 with Isorropia')

plt.ylabel("Concentrations ($\mu g.m^{-3}$)", fontsize=18)
plt.xlabel("time [h]", fontsize=18)
plt.title("Methyltetrols mass")
plt.legend(loc=2,fontsize=10,ncol=2)
plt.savefig("iepox-bimt.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#fig = plt.figure(1,figsize = (10,12))
plt.ylim(ymin=0.)
plt.xlim(xmin=0.)
plt.xlim(xmax=time[Nt-1])
plt.ylim(ymax=17.)

plt.plot(time, org1,'k',label='Equilibrium with soap_inorg=1')
plt.plot(time, org4,'r',label='Equilibrium with soap_inorg=0')
plt.plot(time, org7,'--k',label='Dynamic with soap_inorg=1')
plt.plot(time, org8,'--r',label='Dynamic with soap_inorg=0')

plt.ylabel("Concentrations ($\mu g.m^{-3}$)", fontsize=14)
plt.xlabel("time [h]", fontsize=14)
plt.title("SOA (no NH$_3$)")
plt.legend(loc="lower right",fontsize=10)
plt.savefig("iepox-config.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)
