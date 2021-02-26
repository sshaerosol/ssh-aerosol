import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import *
from matplotlib import rc
import matplotlib.pyplot as plt
from numpy import *
from pylab import *
from matplotlib.ticker import LogLocator, LogFormatter 
from matplotlib import colors 
from matplotlib import ticker
# Turn interactive plotting off
plt.ioff()


################### input parameters #################################"
pcase = '../results/nucl_coupl/'
tag_fig = 'nucl'
sizebin_ssh = 50
time_ssh = 60
pas_temps=1 #time step in min
#####################################################################
number = np.zeros((sizebin_ssh,time_ssh))
cell_diam = np.zeros(sizebin_ssh)
temps = arange(time_ssh)

for j in range(sizebin_ssh) :
        with open (pcase+'/number/NUMBER_' + str(j+1) + '.txt') as finit :
                values = finit.read().splitlines()
        for i in range(time_ssh):
             number[j][i] = float(values[i])
        with open (pcase+'/diameter/DIAMETER_' + str(j+1) + '.txt') as finit :
                values = finit.read().splitlines()
        cell_diam[j] = float(values[0])

# on repasse en echelle log pour les sections alors
# a*log(section1)+b=1 et a*log(sectionN)+b=N
a=(1-(sizebin_ssh))/(log(0.001)-log(10))
b=1-a*log(0.001)
num_y=a*log([0.01,0.1,1])+b
diam_arrondi=['0.01','0.1','1']

tmp1 = -1
tmp2 = -1
deltalogd=[]
diam=[]
with open ('cond_ref/number.init') as f1 :
        values = f1.read().splitlines()
for i in range(len(values)) :
        if i <= (sizebin_ssh + 2) : continue # sizebound + title + star line
        tmp2 = values[i].split('   ', -1)[0]
        if tmp1 == tmp2 : continue
        else :
             diam.append(tmp2.replace('D','E'))
             tmp1 = tmp2
diam.append(10.)
#### diameter
for j in range(len(diam)-1) :
        deltalogd.append(math.log10(float(diam[j+1]))-math.log10(float(diam[j])))

maxnumber = 0.0
for j in range(sizebin_ssh):
   for i in range(time_ssh):
#             number[j][i] = number[j][i]/deltalogd[j]
           if(number[j][i] < 1e3):
                number[j][i] = 1e3
           if(number[j][i] > maxnumber):
                maxnumber = number[j][i] 
print (number.min(),maxnumber,'minmaxnumber')
      
maxnumber = 6.e12      
print (number.min(),maxnumber,'minmaxnumber')

font = {'size'   : 16}
matplotlib.rc('font', **font)

#mftst = 20
#mfts = 20
#mftw = 'bold'

#--- graphique ---#
fig = plt.figure()
#ax = fig.add_subplot(111)

title('Number distribution') # with '+str(int(sizebin_ssh))+' sections')
#cmap='jet'
cmap = plt.get_cmap('PuBuGn')
norm = colors.LogNorm(1e5,maxnumber)
# make an image plot
plt.imshow(number,interpolation='bilinear',aspect='auto',origin='lower',cmap=cmap, norm=norm)
nombre_mn=int(len(temps)*pas_temps)
plt.xlim(0,nombre_mn)
#yticks(num_y, ("0.01","0.1","1"))
plt.ylim(0,a*log(7)+b)
yticks(num_y-(a*log(0.001)+b),diam_arrondi)
nb_ticks = 6
x_set_ticks=[0]
for i in arange(1,nb_ticks,1): 
  x_set_ticks.append(i*nombre_mn/nb_ticks)

xticks(x_set_ticks,x_set_ticks,rotation=0) 
# make a color bar
cb=colorbar(cmap=cmap, norm=norm)
#cb.set_ticks([1e-2,1e0,1e2,1e4,1e6])
#cb.set_ticklabels([r'0', r'$1 \times 10^{-3}$',r'$1 \times 10^{-1}$',r'$10$',r'$1 \times 10^{3}$'])
#cb.set_label(r'd N/d log d_p (${10^3 cm^{-3}}$)')

tick_range = np.logspace(3,13,num=13-3+1, base=10, dtype="int")
cb.set_ticks(tick_range)
xlabel(r"time ($min$)")
ylabel(r"diameter (${\mu m}$)")
plt.tight_layout()
savefig('fig_banana.png')


