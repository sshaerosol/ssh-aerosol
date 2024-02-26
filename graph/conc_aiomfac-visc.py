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
visco9=[]
visco10=[]
visco11=[]
visco12=[]
visco13=[]
visco14=[]
visco15=[]
visco16=[]
visco17=[]
visco18=[]
visco19=[]
visco20=[]
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
namefic9 = '../results/aiomfac-visc/oc2co/aero/'
namefic10 = '../results/aiomfac-visc/oc2ooh/aero/'
namefic11 = '../results/aiomfac-visc/oc2acd/aero/'
namefic12 = '../results/aiomfac-visc/oc2no3/aero/'
namefic13 = '../results/aiomfac-visc/oc2a/aero/'
namefic14 = '../results/aiomfac-visc/oc2b/aero/'
namefic15 = '../results/aiomfac-visc/oc2c/aero/'
namefic16 = '../results/aiomfac-visc/oc2d/aero/'
namefic17 = '../results/aiomfac-visc/oc2e/aero/'
namefic18 = '../results/aiomfac-visc/oc2f/aero/'
namefic19 = '../results/aiomfac-visc/oc2g/aero/'
namefic20 = '../results/aiomfac-visc/oc2h/aero/'

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

with open (namefic9+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco9.append(float(values[i]))

with open (namefic10+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco10.append(float(values[i]))

with open (namefic11+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco11.append(float(values[i]))

with open (namefic12+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco12.append(float(values[i]))

with open (namefic13+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco13.append(float(values[i]))

with open (namefic14+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco14.append(float(values[i]))

with open (namefic15+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco15.append(float(values[i]))

with open (namefic16+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco16.append(float(values[i]))

with open (namefic17+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco17.append(float(values[i]))

with open (namefic18+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco18.append(float(values[i]))

with open (namefic19+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco19.append(float(values[i]))

with open (namefic20+species1) as f1 :
        values = f1.read().splitlines()

for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]             
        values[i] = float(values[i])*1000 
        visco20.append(float(values[i]))
        
#ylim(0, 0.5)
#semilogx(t,visco2,'b',label='10$^{-18}$') #m$^2$ s$^{-1}$')
semilogx(t,visco2,'b',label='OH=1') #m$^2$ s$^{-1}$')
semilogx(t,visco3,'c',label='OH=2') #m$^2$ s$^{-1}$')
semilogx(t,visco4,'g',label='OH=3') #m$^2$ s$^{-1}$')
semilogx(t,visco5,'y',label='OH=4') #m$^2$ s$^{-1}$')
#semilogx(t,visco8,'m',label='C=16, 0=4, RH=30%') #m$^2$ s$^{-1}$')
#semilogx(t,visco7,'r',label='C=16, 0=4, RH=50%') #m$^2$ s$^{-1}$')
#semilogx(t,visco6,'k',label='C=16, 0=4, RH=70%') #m$^2$ s$^{-1}$')
#semilogx(t,visco8,'k',label='10$^{-24}$') #m$^2$ s$^{-1}$')
title("Condensation on C=16 and RH=10%\n with several OH groups",fontsize=16)
xlabel('Time (h)',fontsize=12)
ylabel('POAmP concentrations (ng m$^{-3}$)',fontsize=12)
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
title("Condensation on C=16 and 4 OH groups\n with changing RH",fontsize=16)
xlabel('Time (h)',fontsize=12)
ylabel('POAmP concentrations (ng m$^{-3}$)',fontsize=12)
legend(loc ='best')        
savefig("aiomfac-visc-rh.png")


figure()
#semilogx(t,visco2,'b',label='10$^{-18}$') #m$^2$ s$^{-1}$')
#semilogx(t,visco2,'b',label='C=16, O=1, RH=10%') #m$^2$ s$^{-1}$')
#semilogx(t,visco3,'c',label='C=16, O=2, RH=10%') #m$^2$ s$^{-1}$')
#semilogx(t,visco4,'g',label='C=16, 0=3, RH=10%') #m$^2$ s$^{-1}$')
semilogx(t,visco9,'b',label='1 alcohol + ketone') #m$^2$ s$^{-1}$')
semilogx(t,visco3,'c',label='2 alcohol') #m$^2$ s$^{-1}$')
semilogx(t,visco11,'m',label='1 alcohol + 1 acid') #m$^2$ s$^{-1}$')
semilogx(t,visco10,'r',label='1 alcohol + 1 hydroperoxide') #m$^2$ s$^{-1}$')
semilogx(t,visco12,'k',label='1 alcohol + 1 nitrate') #m$^2$ s$^{-1}$')
#semilogx(t,visco8,'k',label='10$^{-24}$') #m$^2$ s$^{-1}$')
title("Condensation on C=16 and RH=10%\n with different functional groups",fontsize=16)
xlabel('Time (h)',fontsize=12)
ylabel('POAmP concentrations (ng m$^{-3}$)',fontsize=12)
legend(loc ='best')        
savefig("aiomfac-visc-compo.png")

figure()
semilogx(t,visco13,'lightblue',label='C=15') #m$^2$ s$^{-1}$')
semilogx(t,visco14,'b',label='C=16') #m$^2$ s$^{-1}$')
semilogx(t,visco3,'c',label='C=17') #m$^2$ s$^{-1}$')
semilogx(t,visco15,'y',label='C=18') #m$^2$ s$^{-1}$')
semilogx(t,visco16,'g',label='C=19') #m$^2$ s$^{-1}$')
semilogx(t,visco17,'purple',label='C=20') #m$^2$ s$^{-1}$')
semilogx(t,visco18,'m',label='C=21') #m$^2$ s$^{-1}$')
semilogx(t,visco19,'r',label='C=22') #m$^2$ s$^{-1}$')
semilogx(t,visco20,'k',label='C=23') #m$^2$ s$^{-1}$')
#semilogx(t,visco8,'k',label='10$^{-24}$') #m$^2$ s$^{-1}$')
title("Condensation on C=16 and RH=10%\n with different functional groups",fontsize=16)
xlabel('Time (h)',fontsize=12)
ylabel('POAmP concentrations (ng m$^{-3}$)',fontsize=12)
legend(loc ='best')        
savefig("aiomfac-visc-carbons.png")

