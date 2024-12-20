#!/usr/bin/env python3
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
dir2 = '../results/aerosol-vocox/ref-h2onoautox/aero/'
dir3 = '../results/aerosol-vocox/ref-h2onoautox-dNO2/aero/'
dir4 = '../results/aerosol-vocox/ref-h2onoautox-mNO2/aero/'
dir5 = '../results/aerosol-vocox/ref-expl/aero/'
dir6 = '../results/aerosol-vocox/ref-expl-dNO2/aero/'
dir7 = '../results/aerosol-vocox/ref-expl-mNO2/aero/'

Nt= 5 * 60
#####################################################################


font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

shape = [Nt]
tol = np.zeros(shape)
toldNO2 = np.zeros(shape)
tolmNO2 = np.zeros(shape)
tolcompl = np.zeros(shape)
tolcompldNO2 = np.zeros(shape)
tolcomplmNO2 = np.zeros(shape)

species_tol ={'PBiA0D','PBiA1D','PBiA2D','PBiA3D', 'PBiNIT'}
for species in species_tol:
   with open (dir2+'/'+species+'_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  tol[j] = tol[j] + float(string_conc1[j])

for species in species_tol:
   with open (dir3+'/'+species+'_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  toldNO2[j] = toldNO2[j] + float(string_conc1[j])

for species in species_tol:
   with open (dir4+'/'+species+'_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  tolmNO2[j] = tolmNO2[j] + float(string_conc1[j])

#species_tol ={'PLIMALOOH','PNAPINAOOH','PNBPINAOOH','PBPINAOOH','PC98OOH','PC106OOH','PC920PAN','PC97OOH','PC719OOH','PNC102OOH','PLIMALNO3','PNORLIMOOH','PPINIC','PC812OOH','PC98NO3','PH3C25C6PAN','PC621OOH','PC813OOH','PC813NO3','PC516OOH','PC10H14O9','PC10H14O11','PC20H30O13'}
species_tol={'PNOPINONE','PNAPINAOOH','PPINAL','PAPINBNO3','PNAPINBOOH','PAPINANO3','PNC101CO','PAPINAOOH','PAPINBOH','PAPINBOOH','PAPINBCO','PAPINCOOH','PAPINCNO3','PAPINCOH','PNBPINAOOH','PBPINBNO3','PNBPINBOOH','PNC91CHO','PBPINANO3','PBPINAOOH','PBPINAOH','PBPINBOOH','PC918CHO','PBPINCOOH','PBPINCNO3','PBPINCOH','PNLIMOOH','PLIMAL','PLIMBNO3','PLIMAOOH','PLIMANO3','PLIMAOH','PLIMBOOH','PLIMBCO','PLIMCOOH','PLIMCNO3','PLIMKET','PLIMCOH','PC107OOH','PC107OH','PC109OOH','PC109CO','PC109OH','PPINONIC','PC96OOH','PC96NO3','PC96OH','PNORPINAL','PHCC7CO','PNOPINAOOH','PNOPINANO3','PC9DC','PNOPINAOH','PNOPINBOOH','PNOPINBNO3','PNOPINBCO','PNOPINBOH','PNOPINCOOH','PNOPINCNO3','PNOPINCOH','PNOPINDOOH','PNOPINDCO','PNOPINDOH','PLIMALAOOH','PC729CHO','PLIMALACO','PLIMALAOH','PLIMALBOOH','PLIMALBCO','PLIMALBOH','PC923OOH','PC923NO3','PC923OH','PNORLIMAL','PLIMONONIC','PNC101OOH','PPERPINONIC','PC10PAN2','PC720OOH','PC720NO3','PC720OH','PNC91CO3H','PNC91PAN','PC619CO','PC8BCOOH','PC8BCNO3','PC918CO3H','PC918PAN','PLMLKET','PC923CO3H','PC923PAN','PPINALOOH','PPINALNO3','PPINALOH','PC108OOH','PC108NO3','PC108OH','PC89CO2H','PC89CO3H','PC89PAN','PC920CO3H','PHOPINONIC','PC920PAN','PC920OOH','PC97OOH','PC97OH','PC85CO3H','PC9PAN2','PC85OOH','PHYPERACET','PC719OOH','PC719NO3','PC716OH','PC719OH','PC918OOH','PC918NO3','PC918OH','PC9DCOOH','PC9DCNO3','PC9DCCO','PC9DCOH','PC915OOH','PC915NO3','PC88CHO','PC915OH','PC917OOH','PC917NO3','PC917OH','PNLIMALOOH','PNLIMALOH','PLIMALOOH','PLIMALNO3','PLIMALOH','PC626CHO','PC729CO2H','PC729CO3H','PC729PAN','PC822CO2H','PC822CO3H','PC822PAN','PC924OOH','PC622CHO','PC924CO','PC924OH','PC816CO3H','PC816PAN','PNORLIMOOH','PC817CO','PC816OOH','PNLMKAOOH','PLMKBNO3','PLMKAOOH','PLMKANO3','PLMKAOH','PLMKBOOH','PLMKBCO','PCO235C6CHO','PNC102OOH','PCO13C4CHO','PC512CO2H','PC512CO3H','PC512PAN','PC89OOH','PC89NO3','PC89OH','PC926OOH','PCO25C6CHO','PC926OH','PC817CO3H','PKLIMONONIC','PC817PAN','PC817OOH','PC817NO3','PC817OH','PNC826OOH','PNC826OH','PC826OOH','PC826NO3','PC826OH','PC729OOH','PC729NO3','PC816CO','PLMLKAOOH','PLMLKACO','PLMLKAOH','PLMLKBOOH','PLMLKBCO','PLMLKBOH','PC106OOH','PC106NO3','PC106OH','PC717OOH','PC717NO3','PC717OH','PC811CO3H','PPINIC','PC811PAN','PC921OOH','PC98OOH','PC98NO3','PC98OH','PC86OOH','PC919OOH','PC919NO3','PC919OH','PC914OOH','PC914CO','PC914OH','PC916OOH','PC916NO3','PC916OH','PC88CO2H','PC88CO3H','PC88PAN','PC88OOH','PC88CO','PC88OH','PC512OOH','PC512NO3','PC512OH','PC619OOH','PC619OH','PC626CO2H','PC626CO3H','PC626PAN','PC626OOH','PC626NO3','PC511CHO','PC517CHO','PC735OOH','PC735OH','PC822OOH','PC822NO3','PC822OH','PC823CO3H','PLIMONIC','PC823PAN','PC925OOH','PC818CO','PC622CO2H','PC622CO3H','PC622PAN','PMEKAOOH','PMEKANO3','PC235C6CO3H','PC7PAN3','PCHOC3COOOH','PCHOC3COPAN','PNC71OOH','PC811OOH','PC811NO3','PC721CHO','PC811OH','PC413COOOH','PC614CO','PH3C25CCO2H','PH3C25CCO3H','PH3C25C6PAN','PH3C25C6OOH','PH3C25C5CHO','PH3C25C6OH','PHCOCH2CO2H','PHCOCH2CO3H','PCO123C5CHO','PCO13C3CO2H','PC810OOH','PC810NO3','PC810OH','PC818OOH','PC818OH','PC727CO3H','PC727PAN','PNC728OOH','PNC728OH','PC728OOH','PC728NO3','PC728OH','PC622OOH','PC622NO3','PC518CHO','PC622OH','PC823OOH','PC823NO3','PC823CO','PC823OH','PC819OOH','PMACROOH','PMACRNO3','PMACROH','PMACROHOOH','PMACRNB','PC3MDIALOH','PC731CO2H','PC731CO3H','PC731PAN','PCO235C6OOH','PCO235C5CHO','PC4CODIAL','PC716OOH','PNC71CO','PC922OOH','PC614OOH','PC614NO3','PC614OH','PC511OOH','PC620OOH','PC515CHO','PC620OH','PC87CO2H','PC87CO3H','PC87PAN','PC616OOH','PCO12C4CHO','PC616OH','PC718CO2H','PC718CO3H','PC718PAN','PC513OOH','PC513CO','PC513OH','PCO25C6CO2H','PCO25C6CO3H','PC627PAN','PC627OOH','PC627OH','PCO2C4GLYOX','PC727OOH','PC727CO','PC511CO3H','PC511PAN','PC517CO2H','PC517CO3H','PC517PAN','PC517OOH','PC517NO3','PC517OH','PHMVKBCHO','PC628OOH','PCO13C4OH','PC628OH','PC824OOH','PC624CHO','PC824CO','PC824OH','PMAE','PCONM2CHO','PHOCO3C4OOH','PHO14CO2C4','PHO1CO3CHO','PCO2C3CO3H','PCO2C3PAN','PHYETHO2H','PETHOHNO3','PETHGLY','PC812OOH','PC812OH','PC721CO3H','PNORPINIC','PC721PAN','PC721OOH','PBIACETOOH','PBIACETOH','PHOCH2CO2H','PHOCH2CO3H','PPHAN','PH3C2C4CO2H','PH3C2C4CO3H','PH3C2C4PAN','PHMVKAOOH','PHMVKANO3','PCO2H3CHO','PHO12CO3C4','PC87OOH','PC87CO','PC87OH','PCHOC2CO2H','PCHOC2CO3H','PCHOC2PAN','PC718OOH','PC718NO3','PC617CHO','PC718OH','PC514OOH','PC514NO3','PC514OH','PC820OOH','PC518CO2H','PC518CO3H','PC518PAN','PNC623OOH','PNC623OH','PC623OOH','PC623NO3','PC623OH','PISOPDOOH','PISOPDNO3','PHCOC5','PISOPDOH','PC825OOH','PC825CO','PC825OH','PCO2C4CO2H','PCO2C4CO3H','PC5PAN2','PMACRNCO2H','PMACRNCO3H','PMACRNPAN','PMACRNBCO2H','PMACRNBCO3H','PMACRNBPAN','PCHOMOHCO3H','PCHOMOHPAN','PC731OOH','PC731NO3','PC731OH','PC732CO3H','PKLIMONIC','PC732PAN','PCO23C4CO3H','PC5PAN9','PC312COCO3H','PC312COPAN','PALCOCH2OOH','PNC72OOH','PC621OOH','PH1C23C4CHO','PHCOCO2H','PC515CO3H','PC515PAN','PC515OOH','PCHOC2H4OOH','PHOC2H4CO2H','PHOC2H4CO3H','PC3PAN1','PC821OOH','PHMVKBCO2H','PHMVKBCO3H','PHMVKBPAN','PMVKNO3','PHMVKBOOH','PC520OOH','PHOCH2COCHO','PC520OH','PIEPOXB','PC519CHO','PC624CO2H','PC624CO3H','PC624PAN','PCONM2CO2H','PCONM2CO3H','PCONM2PAN','PIPRHOCO2H','PIPRHOCO3H','PC4PAN5','PCOHM2CO2H','PCOHM2CO3H','PCOHM2PAN','PC732OOH','PC732NO3','PC732CO','PC732OH','PC813OOH','PC813NO3','PC813OH','PC722OOH','PCO2H3CO3H','PC4PAN6','PC515CO','PC615CO2H','PC615CO3H','PC615PAN','PC617CO2H','PC617CO3H','PC617PAN','PC618CO2H','PC618CO3H','PC618PAN','PC617OOH','PC615CO','PC617OH','PC618OOH','PC67CHO','PCO2N3CHO','PIEB1CHO','PIEB4CHO','PINDOOH','PINB1NO3','PINDOH','PISOPCOOH','PISOPCNO3','PHC4ACHO','PC5HPALD2','PHC4CCHO','PISOPAOH','PC59OOH','PNC730OOH','PNC730OH','PC730OOH','PC730NO3','PC730OH','PC624OOH','PC624NO3','PC624CO','PC624OH','PCH3COCO2H','PCH3COPAN','PC733OOH','PC733CO','PC733OH','PNC61CO3H','PNC6PAN1','PC615OOH','PC615OH','PC57OOH','PC57NO3','PC57OH','PC58AOOH','PC58ANO3','PINDHPCHO','PINB1NACHO','PINB1NBCHO','PINDHCHO','PIEPOXC','PC537OOH','PDHPMPAL','PC5PACALD2','PC4MDIAL','PIEC2OOH','PC519CO2H','PC519CO3H','PC519PAN','PC629OOH','PHO1CO34C5','PC629OH','PC734OOH','PC734CO','PC734OH','PC516OOH','PC44OOH','PHC23C4CO3H','PH1C23C4PAN','PH1C23C4OOH','PCO1M22CO2H','PCO1M22CO3H','PCO1M22PAN','PCO2N3CO3H','PCO2N3PAN','PC4M2ALOHNO3','PC47CHO','PC4MALOHOOH','PC4M2AL2OH','PIECCHO','PC527OOH','PC527NO3','PINCOOH','PINCNO3','PINCCO','PINCOH','PC4CO2OOH','PC3MDIALOOH','PPXYFUONE','PHPC52OOH','PHC4CCO2H','PHC4CCO3H','PC5PAN19','PC5PACALD1','PC57AOOH','PC519OOH','PHO13CO4C5','PC625OOH','PC625OH','PH1CO23CHO','PHC4ACO2H','PHC4ACO3H','PC5PAN17','PC58OOH','PC58NO3','PC58OH','PINCNCHO','PINCGLYOX','PC535OOH','PC3MCODBPAN','PMC3ODBCO2H','PMC3CODBPAN','PC531OOH','PC58NO3CO2H','PC58NO3CO3H','PC58NO3PAN','PC531CO','PC31CO3H','PC31PAN','PC533OOH','PC10H14O6','PC10H15O5NO3','PC10H14O7','PC10H15O6NO3','PC10H14O8','PC10H15O7NO3','PC10H14O9','PC10H15O8NO3','PC10H14O10','PC10H15O9NO3','PC10H14O11','PC10H15O10NO3','PC10H14O12','PC10H15O11NO3','PC10H14O13','PC10H15O12NO3','PC10H16O4iso1','PC10H16O5iso1','PC10H16O6iso1','PC10H16O7iso1','PC10H16O8iso1','PC10H16O9iso1','PC10H16O10','PC10H16O11','PC10H16O12','PC10H16O13','PC10H16O14','PC10H16O3','PC10H14O3','PC20H30O13','PC10H14O4','PC10H14O5','PC10H16O5iso2','PC10H17O4NO3','PC10H16O6iso2','PC10H17O5NO3','PC10H16O7iso2','PC10H17O6NO3','PC10H16O8iso2','PC10H17O7NO3','PC10H16O9iso2','PC10H17O8NO3','PC10H18O5','PC10H18O6','PC10H18O7','PC10H18O8','PC10H18O9','PC10H18O10','PC10H18O4','PC10H16O4iso2','PC20H34O10'}
for species in species_tol:
   with open (dir5+'/'+species+'_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  tolcompl[j] = tolcompl[j] + float(string_conc1[j])

for species in species_tol:
   with open (dir6+'/'+species+'_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  tolcompldNO2[j] = tolcompldNO2[j] + float(string_conc1[j])

for species in species_tol:
   with open (dir7+'/'+species+'_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  tolcomplmNO2[j] = tolcomplmNO2[j] + float(string_conc1[j])


time=range(60,60*(Nt+1),60)
time=numpy.array(time)
time=time/3600.0

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
fig = plt.figure(1,figsize = (15,15))
plt.ylim(0.,1.5)

#plt.plot(time, totsoa,'k-',label='total')
plt.plot(time, tol,color='k',label='H$^2$O - No autox.- Ref')
plt.plot(time, tolmNO2,'m',label='H$^2$O - No autox.- NO$_2$ x 2')
plt.plot(time, toldNO2,'b',label='H$^2$O - No autox.- NO$_2$ / 2')
plt.plot(time, tolcompl,color='k', linestyle='dashed',label='Expl - Ref')
plt.plot(time, tolcomplmNO2,'m',linestyle='dashed',label='Expl - NO$_2$ x 2')
plt.plot(time, tolcompldNO2,'b',linestyle='dashed',label='Expl - NO$_2$ / 2')
plt.ylabel("$\mu g.m^{-3}$", fontsize=18)
plt.xlabel("time [h]", fontsize=18)
plt.legend(loc=2)
plt.savefig("particles-vocox-monoterp-h2o-noautox.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)

