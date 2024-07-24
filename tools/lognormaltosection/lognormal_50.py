import numpy as np
from scipy.special import erf

import matplotlib.pyplot as plt

def lognormal_N_dn(bounds, N, dn_mean, sigma):
    """
    Lognormal distribution density for number concentration, computed from total number and mean (number) diameter.
    S.I. units

    n(d) = \frac{N}{\sqrt{2\pi} \log{\sigma} d} \exp{-\frac{\log{d-\bar{d}_n}}{2 \log{\sigma}}}
    """
    d = np.sqrt(bounds[1:]*bounds[:-1])
    norm = N / (np.sqrt(2*np.pi) * np.log(sigma))
    x = norm * 1/d * np.exp(-0.5*((np.log(d)-np.log(dn_mean))/np.log(sigma))**2)
    return x

def lognormal_V_dv(bounds, V, dv_mean, sigma):
    """
    Lognormal distribution density for number concentration, computed from total volume and mean (volume) diameter.
    S.I. units

    n(d) = \frac{N}{\sqrt{2\pi} \log{\sigma} d} \exp{-\frac{\log{d-\bar{d}_n}}{2 \log{\sigma}}}
    """
    d = np.sqrt(bounds[1:]*bounds[:-1])
    norm = 6./np.pi * V / dv_mean**3 * np.exp(4.5 * np.log(sigma)**2) / (np.sqrt(2*np.pi) * np.log(sigma))
    x = norm / d * np.exp(-0.5 * ((np.log(d)-np.log(dv_mean)+3*np.log(sigma)**2)/np.log(sigma))**2)
    return x

def int_lognormal_N_dn(bounds, N, dn_mean, sigma):
    """
    Cumulative lognormal distribution for number concentration, computed from total number and mean (number) diameter.
    S.I. units

    N(d_-, d_+) = [\frac{1}{2} + \frac{1}{2} \erf{\frac{\log{d} - \log{\bar{d}_n}}{\sqrt{2} \log{\sigma}}}]_{d_-}^{d_+}
    """
    norm = 0.5 * N
    c = norm * erf((np.log(bounds)-np.log(dn_mean))/(np.log(sigma)*np.sqrt(2.)))
    x = c[1:]-c[:-1]
    return x

def int_lognormal_V_dv(bounds, V, dv_mean, sigma):
    """
    Cumulative lognormal distribution for number concentration, computed from total volume and mean (volume) diameter.
    S.I. units

    N(d_-, d_+) = [\frac{V}{2} \frac{6}{\pi \bar{d}_v^3} \exp{\frac{9}{2} \log{\sigma}^2} (1+\erf{\frac{\log{d} - \log{\bar{d}_v} + 3\log{\sigma}^2}{\sqrt{2} \log{\sigma}}})]_{d_-}^{d_+}
    """
    norm = 0.5 * 6/np.pi * V / dv_mean**3 * np.exp(4.5 * np.log(sigma)**2)
    c = norm * erf((np.log(bounds)-np.log(dv_mean)+3*np.log(sigma)**2)/(np.log(sigma)*np.sqrt(2)))
    x = c[1:]-c[:-1]
    return x

# Ref : C. Seigneur et al.,
#       Simulation of Aerosol Dynamics: A Comparative Review of Mathematical Models,
#       1986, Aerosol Science and Technology

# Urban lognormal parametrization

dv_n = 3.8e-8
dv_a = 3.2e-7
dv_c = 5.7e-6

sigma_n = 1.8
sigma_a = 2.16
sigma_c = 2.21

V_n = 0.63e-12
V_a = 38.4e-12
V_c = 30.8e-12

# Geometrical spacing between 1 nm and 10 \mu m
nbound=51
bounds = 10 ** np.linspace(-9, -5, num=51)
#print (bounds)

# Comparaison between density interpolation and actual cumulative distribution

#density = lognormal_V_dv(bounds, V_a, dv_a, sigma_a)
#cumulative = int_lognormal_V_dv(bounds, V_a, dv_a, sigma_a)

#fig = plt.figure()
#plt.plot(np.sqrt(bounds[:-1]*bounds[1:]), cumulative, linestyle="solid", label="Cumulative")
#plt.plot(np.sqrt(bounds[:-1]*bounds[1:]), density * (bounds[1:]-bounds[:-1]), linestyle="dashed", label="Density")
#plt.xscale("log")
#plt.yscale("log")
#plt.legend()

#plt.savefig("lognormal.png", dpi=200)

# Comparison with SSHaerosol coagulation test case initial condition

density = lognormal_V_dv(bounds, V_n, dv_n, sigma_n) + lognormal_V_dv(bounds, V_a, dv_a, sigma_a) + lognormal_V_dv(bounds, V_c, dv_c, sigma_c)
cumulative = int_lognormal_V_dv(bounds, V_n, dv_n, sigma_n) + int_lognormal_V_dv(bounds, V_a, dv_a, sigma_a) + int_lognormal_V_dv(bounds, V_c, dv_c, sigma_c)

cumulative_n = int_lognormal_V_dv(bounds, V_n, dv_n, sigma_n)
cumulative_a = int_lognormal_V_dv(bounds, V_a, dv_a, sigma_a)
cumulative_c = int_lognormal_V_dv(bounds, V_c, dv_c, sigma_c)

#print(cumulative)
#for i in range(nbound-1):
#  print(cumulative[i]*np.sqrt(bounds[i]*bounds[i+1])**3 * 1.84*np.pi/6.*1e+12)

with open('init_num_coag50.dat') as f:
    f.readline()
    rawdata = f.readline()
    data = rawdata.split("\t")

reference = [float(d) for d in data[1:]]

fig = plt.figure()
plt.plot(np.sqrt(bounds[:-1]*bounds[1:]), reference, linestyle="solid", color="black", label="Reference")
plt.plot(np.sqrt(bounds[:-1]*bounds[1:]), cumulative, linestyle="dashed", label="Cumulative distribution")
plt.plot(np.sqrt(bounds[:-1]*bounds[1:]), cumulative_n, linestyle="dashed", alpha=0.3, label="Nucleation")
plt.plot(np.sqrt(bounds[:-1]*bounds[1:]), cumulative_a, linestyle="dashed", alpha=0.3, label="Accumulation")
plt.plot(np.sqrt(bounds[:-1]*bounds[1:]), cumulative_c, linestyle="dashed", alpha=0.3, label="Coarse")
plt.xscale("log")
plt.yscale("log")
plt.ylim(1e1, 1e11)
plt.legend()

plt.savefig("urban.png", dpi=200)

# Hazy lognormal parametrization

dv_n = 4.4e-8
dv_a = 2.4e-7
dv_c = 6.0e-6

sigma_n = 1.2
sigma_a = 1.8
sigma_c = 2.2

V_n = 0.09e-12
V_a = 5.8e-12
V_c = 25.9e-12

cumulative = int_lognormal_V_dv(bounds, V_n, dv_n, sigma_n) + int_lognormal_V_dv(bounds, V_a, dv_a, sigma_a) + int_lognormal_V_dv(bounds, V_c, dv_c, sigma_c)

#print(cumulative)
#for i in range(nbound-1):
#  print(bounds[i])
#  print(cumulative[i]*np.sqrt(bounds[i]*bounds[i+1])**3 * 1.84*np.pi/6.*1e+12)

cumulative_n = int_lognormal_V_dv(bounds, V_n, dv_n, sigma_n)
cumulative_a = int_lognormal_V_dv(bounds, V_a, dv_a, sigma_a)
cumulative_c = int_lognormal_V_dv(bounds, V_c, dv_c, sigma_c)

fig = plt.figure()
plt.plot(np.sqrt(bounds[:-1]*bounds[1:]), cumulative, linestyle="dashed", label="Cumulative distribution")
plt.plot(np.sqrt(bounds[:-1]*bounds[1:]), cumulative_n, linestyle="dashed", alpha=0.3, label="Nucleation")
plt.plot(np.sqrt(bounds[:-1]*bounds[1:]), cumulative_a, linestyle="dashed", alpha=0.3, label="Accumulation")
plt.plot(np.sqrt(bounds[:-1]*bounds[1:]), cumulative_c, linestyle="dashed", alpha=0.3, label="Coarse")
plt.xscale("log")
plt.yscale("log")
plt.ylim(1e1, 1e10)
plt.legend()

plt.savefig("hazy.png", dpi=200)
