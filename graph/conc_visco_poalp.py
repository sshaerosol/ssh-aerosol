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
from atmopy import *
import matplotlib.pyplot as plt
# Turn interactive plotting off
plt.ioff()

sizebin = 1
visco1=[]
visco2=[]
visco3=[]
visco4=[]
visco5=[]
t=[]
t0=0
deltat=600.
namefic1 = '../results/visco/visco1_poalp/aero/'
namefic2 = '../results/visco/visco2_poalp/aero/'
namefic3 = '../results/visco/visco3_poalp/aero/'
namefic4 = '../results/visco/visco4_poalp/aero/'
namefic5 = '../results/visco/visco5_poalp/aero/'
namefic6 = '../results/visco/visco6_poalp/aero/'
namefic7 = '../results/visco/visco7_poalp/aero/'
namefic8 = '../results/visco/visco8_poalp/aero/'


species1 = 'PPOAlP_1.txt'
species2 = 'PSOAlP_1.txt'
with open (namefic1+species1) as f1 :
	values = f1.read().splitlines()
with open (namefic1+species2) as f2 :
	values2 = f2.read().splitlines()
for i in range(len(values)) :
	tmp = values[i].split('   ', -1)[0]        
	tmp2 = values2[i].split('   ', -1)[0]      
        values[i] = float(values[i]) + float(values2[i])  
        visco1.append(float(values[i]))
        t.append(t0)
        t0=t0+deltat*1./3600

with open (namefic2+species1) as f1 :
	values = f1.read().splitlines()
with open (namefic2+species2) as f2 :
	values2 = f2.read().splitlines()
for i in range(len(values)) :
	tmp = values[i].split('   ', -1)[0]        
	tmp2 = values2[i].split('   ', -1)[0]      
        values[i] = float(values[i]) + float(values2[i])  
        visco2.append(float(values[i]))

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

semilogx(t,visco1,'0.5')
semilogx(t,visco2,'b')
semilogx(t,visco3,'c')
semilogx(t,visco4,'g')
semilogx(t,visco5,'y')
semilogx(t,visco6,'m')
semilogx(t,visco7,'r')
semilogx(t,visco8,'k')
savefig("visco_poalp.png")



