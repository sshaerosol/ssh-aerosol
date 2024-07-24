import upsample
import numpy as np

M = 1.
rho = 1.

dmin = 1e-1
dmax = 1e1

print('Arbitrary unit system')
print('>>> Parameters')
print('Total mass :', M)
print('Density :', rho)
print('Diameter range : from {:.2e} to {:.2e}'.format(dmin, dmax))

print('\nChecking method consistency with one section')

M_1, N_1 = upsample.M_N_conservative_upsample(M, rho, np.array([dmin, dmax]))

print('Total mass', M_1[0])
print('Total number', N_1[0])

print('\nSpliting section at geometric mean diameter in 2 logspaced sections')

M_log2, N_log2 = upsample.M_N_conservative_upsample_logspace(M, rho, dmin, dmax, 2)

print('Mass concentrations', M_log2)
print('Sum of mass concentrations', np.sum(M_log2))
print('Number concentrations', N_log2)
print('Sum of number concentrations', np.sum(N_log2))

print('\nSpliting section in 3 logspaced sections')

M_log3, N_log3 = upsample.M_N_conservative_upsample_logspace(M, rho, dmin, dmax, 3)

print('Mass concentrations', M_log3)
print('Sum of mass concentrations', np.sum(M_log3))
print('Number concentrations', N_log3)
print('Sum of number concentrations', np.sum(N_log3))

print('\nSpliting section in 12 logspaced sections')

M_log12, N_log12 = upsample.M_N_conservative_upsample_logspace(M, rho, dmin, dmax, 12)

print('Mass concentrations', M_log12)
print('Sum of mass concentrations', np.sum(M_log12))
print('Number concentrations', N_log12)
print('Sum of number concentrations', np.sum(N_log12))
