#!/usr/bin/env python3
# -*- coding: iso-8859-1  -*-
import matplotlib
matplotlib.use('Agg')
import math,sys, optparse, glob, os, re, string
from datetime import *
from math import *
from pylab import *
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
# Turn interactive plotting off
plt.ioff()

sizebin = 1
visco2=[]
visco3=[]
visco4=[]
visco5=[]
visco6=[]
visco7=[]
visco8=[]

t=[]
t0=0
deltat=600.
#namefic2 = '../results/aiomfac-visc/oc0/aero/'
namefic2 = '../results/aiomfac-visc/oc1/aero/'
namefic3 = '../results/aiomfac-visc/oc2/aero/'
namefic4 = '../results/aiomfac-visc/oc3/aero/'
namefic5 = '../results/aiomfac-visc/oc4/aero/'
namefic6 = '../results/aiomfac-visc/oc4b/aero/'
namefic7 = '../results/aiomfac-visc/oc4c/aero/'
namefic8 = '../results/aiomfac-visc/oc4d/aero/'

species1 = 'PPOAmP_1.txt'
#species2 = 'PSOAlP_1.txt'

font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

with open (namefic2+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]       
        values[i] = float(values[i])*1000   
        visco2.append(float(values[i]))
        t.append(t0)
        t0=t0+deltat*1./3600

with open (namefic3+species1) as f1 :
        values = f1.read().splitlines()
        
for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]         
        values[i] = float(values[i])*1000   
        visco3.append(float(values[i]))
        
with open (namefic4+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]       
        values[i] = float(values[i])*1000 
        visco4.append(float(values[i]))
        
with open (namefic5+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]        
        values[i] = float(values[i])*1000   
        visco5.append(float(values[i]))

with open (namefic6+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]        
        values[i] = float(values[i])*1000 
        visco6.append(float(values[i]))

with open (namefic7+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]        
        values[i] = float(values[i])*1000 
        visco7.append(float(values[i]))

with open (namefic8+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco8.append(float(values[i]))

#ylim(0, 0.5)
#semilogx(t,visco2,'b',label='10$^{-18}$') #m$^2$ s$^{-1}$')
semilogx(t,visco2,'b',label='O=1') #m$^2$ s$^{-1}$')
semilogx(t,visco3,'c',label='O=2') #m$^2$ s$^{-1}$')
semilogx(t,visco4,'g',label='O=3') #m$^2$ s$^{-1}$')
semilogx(t,visco5,'y',label='O=4') #m$^2$ s$^{-1}$')
#semilogx(t,visco8,'m',label='C=16, 0=4, RH=30%') #m$^2$ s$^{-1}$')
#semilogx(t,visco7,'r',label='C=16, 0=4, RH=50%') #m$^2$ s$^{-1}$')
#semilogx(t,visco6,'k',label='C=16, 0=4, RH=70%') #m$^2$ s$^{-1}$')
#semilogx(t,visco8,'k',label='10$^{-24}$') #m$^2$ s$^{-1}$')
title("Condensation on C=16 and RH=10% with changing O",fontsize=16)
xlabel('Time (h)')
ylabel('POAmP concentrations (ng m$^{-3}$)')
legend(loc ='best')        
savefig("aiomfac-visc-oxygens.png")

figure()
#semilogx(t,visco2,'b',label='10$^{-18}$') #m$^2$ s$^{-1}$')
#semilogx(t,visco2,'b',label='C=16, O=1, RH=10%') #m$^2$ s$^{-1}$')
#semilogx(t,visco3,'c',label='C=16, O=2, RH=10%') #m$^2$ s$^{-1}$')
#semilogx(t,visco4,'g',label='C=16, 0=3, RH=10%') #m$^2$ s$^{-1}$')
semilogx(t,visco5,'y',label='RH=10%') #m$^2$ s$^{-1}$')
semilogx(t,visco8,'m',label='RH=30%') #m$^2$ s$^{-1}$')
semilogx(t,visco7,'r',label='RH=50%') #m$^2$ s$^{-1}$')
semilogx(t,visco6,'k',label='RH=70%') #m$^2$ s$^{-1}$')
#semilogx(t,visco8,'k',label='10$^{-24}$') #m$^2$ s$^{-1}$')
title("Condensation on C=16 and O=4 with changing RH",fontsize=16)
xlabel('Time (h)')
ylabel('POAmP concentrations (ng m$^{-3}$)')
legend(loc ='best')        
savefig("aiomfac-visc-rh.png")



