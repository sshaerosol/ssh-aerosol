#!/usr/bin/env python
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
namefic2 = '../results/visco/visco2_poalp/aero/'
namefic3 = '../results/visco/visco3_poalp/aero/'
namefic4 = '../results/visco/visco4_poalp/aero/'
namefic5 = '../results/visco/visco5_poalp/aero/'
namefic6 = '../results/visco/visco6_poalp/aero/'
namefic7 = '../results/visco/visco7_poalp/aero/'
namefic8 = '../results/visco/visco8_poalp/aero/'


species1 = 'PPOAlP_1.txt'
species2 = 'PSOAlP_1.txt'

font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

with open (namefic2+species1) as f1 :
        values = f1.read().splitlines()
with open (namefic2+species2) as f2 :
        values2 = f2.read().splitlines()
for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]        
        tmp2 = values2[i].split('   ', -1)[0]      
        values[i] = float(values[i]) + float(values2[i])  
        visco2.append(float(values[i]))
        t.append(t0)
        t0=t0+deltat*1./3600

with open (namefic3+species1) as f1 :
        values = f1.read().splitlines()
with open (namefic3+species2) as f2 :
        values2 = f2.read().splitlines()
for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]        
        tmp2 = values2[i].split('   ', -1)[0]      
        values[i] = float(values[i]) + float(values2[i])  
        visco3.append(float(values[i]))
        
with open (namefic4+species1) as f1 :
        values = f1.read().splitlines()
with open (namefic4+species2) as f2 :
        values2 = f2.read().splitlines()
for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]        
        tmp2 = values2[i].split('   ', -1)[0]      
        values[i] = float(values[i]) + float(values2[i])  
        visco4.append(float(values[i]))
        
with open (namefic5+species1) as f1 :
        values = f1.read().splitlines()
with open (namefic5+species2) as f2 :
        values2 = f2.read().splitlines()
for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]        
        tmp2 = values2[i].split('   ', -1)[0]      
        values[i] = float(values[i]) + float(values2[i])  
        visco5.append(float(values[i]))

with open (namefic6+species1) as f1 :
        values = f1.read().splitlines()
with open (namefic6+species2) as f2 :
        values2 = f2.read().splitlines()
for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]        
        tmp2 = values2[i].split('   ', -1)[0]      
        values[i] = float(values[i]) + float(values2[i])  
        visco6.append(float(values[i]))

with open (namefic7+species1) as f1 :
        values = f1.read().splitlines()
with open (namefic7+species2) as f2 :
        values2 = f2.read().splitlines()
for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]        
        tmp2 = values2[i].split('   ', -1)[0]      
        values[i] = float(values[i]) + float(values2[i])  
        visco7.append(float(values[i]))

with open (namefic8+species1) as f1 :
        values = f1.read().splitlines()
with open (namefic8+species2) as f2 :
        values2 = f2.read().splitlines()
for i in range(len(values)) :
        tmp = values[i].split('   ', -1)[0]        
        tmp2 = values2[i].split('   ', -1)[0]      
        values[i] = float(values[i]) + float(values2[i])  
        visco8.append(float(values[i]))

ylim(4.9, 10)
semilogx(t,visco2,'b',label='10$^{-18}$')
semilogx(t,visco3,'c',label='10$^{-19}$')
semilogx(t,visco4,'g',label='10$^{-20}$')
semilogx(t,visco5,'y',label='10$^{-21}$')
semilogx(t,visco6,'m',label='10$^{-22}$')
semilogx(t,visco7,'r',label='10$^{-23}$')
semilogx(t,visco8,'k',label='10$^{-24}$')
title("Kp $\simeq$ 1")
xlabel('Time (s)')
ylabel('Organic concentrations ($\mu$g m$^{-3}$)')
legend(loc ='best')        
plt.tight_layout()
savefig("visco_poalp.png")



