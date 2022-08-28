#!/usr/bin/env python3
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
direc = ['../results/cond-caco3-dyn/aero/','../results/cond-caco3-dyn-hno3/aero/','../results/cond-caco3-eq/aero/','../results/cond-caco3-eq-hno3/aero/']

Nt= 180 #5 * 60
#####################################################################


font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

Nsim=len(direc)

shape = [Nt,Nsim]
total=np.zeros(shape)
co3=np.zeros(shape)
no3=np.zeros(shape)
water=np.zeros(shape)
ca=np.zeros(shape)

time=range(10,10*(Nt+1),10)
time=numpy.array(time)
time=time/3600.0

fig = plt.figure()

nbins=15

#ax = fig.add_subplot(1,1,1)
fig = plt.figure(1,figsize = (15,15))
#plt.ylim(0.,80.)

type_line=['','-.','--',':']

for isim in range(0,Nsim):
        shape = [Nt]
        total=np.zeros(shape)
        co3=np.zeros(shape)
        no3=np.zeros(shape)
        water=np.zeros(shape)
        ca=np.zeros(shape)
        dir1=direc[isim]
        for ibins in range(0,nbins):
                with open (dir1+'/PCO3_'+str(ibins+1)+'.txt') as finit :
                        string_conc1 = finit.read().splitlines()
                        for j in range(Nt) :
                                co3[j] = co3[j] + float(string_conc1[j])
                with open (dir1+'/PNO3_'+str(ibins+1)+'.txt') as finit :
                        string_conc1 = finit.read().splitlines()
                        for j in range(Nt) :
                                no3[j] = no3[j] + float(string_conc1[j])
                with open (dir1+'/PH2O_'+str(ibins+1)+'.txt') as finit :
                        string_conc1 = finit.read().splitlines()
                        for j in range(Nt) :
                                water[j] = water[j] + float(string_conc1[j])
                with open (dir1+'/PCa_'+str(ibins+1)+'.txt') as finit :
                        string_conc1 = finit.read().splitlines()
                        for j in range(Nt) :
                                ca[j] = ca[j] + float(string_conc1[j])                  

        if (isim==0):
                plt.plot(time, co3+no3+water+ca,'k'+type_line[isim],label='Total dynamic - no HNO3')
        elif (isim==1):
                plt.plot(time, co3+no3+water+ca,'k'+type_line[isim],label='Total dynamic with HNO3')
        elif (isim==2):
                plt.plot(time, co3+no3+water+ca,'k'+type_line[isim],label='Total equilibrium - no HNO3')
        else:
                plt.plot(time, co3+no3+water+ca,'k'+type_line[isim],label='Total equilibrium with HNO3')      
        plt.plot(time, no3,'r'+type_line[isim],label='Nitrate')
        plt.plot(time, co3,'g'+type_line[isim],label='Carbonates')
        plt.plot(time, water,'b'+type_line[isim],label='Water')
        
plt.ylabel("$\mu g.m^{-3}$", fontsize=18)
plt.xlabel("time [h]", fontsize=18)
plt.legend(prop={'size': 6})
plt.savefig("caco3.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)


