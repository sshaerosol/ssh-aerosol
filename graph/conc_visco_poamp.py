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
namefic2 = '../results/visco/visco2_poamp/aero/'
namefic3 = '../results/visco/visco3_poamp/aero/'
namefic4 = '../results/visco/visco4_poamp/aero/'
namefic5 = '../results/visco/visco5_poamp/aero/'
namefic6 = '../results/visco/visco6_poamp/aero/'
namefic7 = '../results/visco/visco7_poamp/aero/'
namefic8 = '../results/visco/visco8_poamp/aero/'


species1 = 'PPOAmP_1.txt'
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

ylim(4.9, 5.5)
semilogx(t,visco2,'b',label='10$^{-18}$') #m$^2$ s$^{-1}$')
semilogx(t,visco3,'c',label='10$^{-19}$') #m$^2$ s$^{-1}$')
semilogx(t,visco4,'g',label='10$^{-20}$') #m$^2$ s$^{-1}$')
semilogx(t,visco5,'y',label='10$^{-21}$') #m$^2$ s$^{-1}$')
semilogx(t,visco6,'m',label='10$^{-22}$') #m$^2$ s$^{-1}$')
semilogx(t,visco7,'r',label='10$^{-23}$') #m$^2$ s$^{-1}$')
semilogx(t,visco8,'k',label='10$^{-24}$') #m$^2$ s$^{-1}$')
title("Kp $\simeq$ 0.01")
xlabel('Time (h)')
ylabel('Organic concentrations ($\mu$g m$^{-3}$)')
legend(loc ='best')        
savefig("visco_poamp.png")



