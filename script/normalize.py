import h5py as h5
import numpy as np
from scipy import *
# Read Baryon density
f = h5.File('./ic.enzo/GridDensity', 'r')
dataset = f['GridDensity']
dens = np.array(dataset[0])

baryonfrac = 0.0441/0.268
# renormalize so mean is 1
ave = dens.mean()
dens = dens/ave*(1-baryonfrac)
print dens.shape
print 'average = ',ave
print 'new average = ', dens.mean()

# write out new density to new file
#f1 = h5.File('./ic.enzo/GridDensity.new', 'w')
#new_dens = f1.create_dataset('GridDensity.new', data=dens)
#for a in dataset.attrs.keys():
#   new_dens.attrs.create(a, dataset.attrs[a])
#f1.close()

# write out FDM density to new file
f1 = h5.File('./ic.enzo/GridFDMDensity', 'w')
new_dens = f1.create_dataset('GridFDMDensity', data=dens)
for a in dataset.attrs.keys():
   new_dens.attrs.create(a, dataset.attrs[a])
f1.close()

# unit and constants
BoxLength = 10.0
hbar = 1.05457266e-27
mass_unit = 0.1
mass = 1e-22*1.6021772e-12/(2.99792458e10**2)*mass_unit 
print mass_unit

InitialRedshift = 100.
HubbleConstantNow = 0.704
OmegaMatterNow = 0.268

a0 = 1./(1+InitialRedshift)
InitialTime = 0.81651316219217
LengthUnits = 3.085678e24*BoxLength/HubbleConstantNow/(1 + InitialRedshift)
TimeUnits = 2.519445e17/sqrt(OmegaMatterNow)/HubbleConstantNow/(1 + InitialRedshift)**1.5
acoef = 1.5**(1./3.)*a0
coef = hbar/mass*TimeUnits/(LengthUnits**2)
print coef

# solve poisson equation for theta
N = dens.shape[-1]
# use fft to solve theta
LHS = -2./3.*(dens-1)/InitialTime/coef

klhs = np.fft.fftn(LHS)
G1d = N * np.fft.fftfreq(N) * 2*pi
kx, ky, kz = meshgrid(G1d, G1d, G1d, indexing='ij')

G2 = kx**2 + ky**2 + kz**2
G2[0,0,0] = 1.

thetak = klhs/(-G2)
thetak[0,0,0] = 0
theta = real(np.fft.ifftn(thetak))
#calculate wave function
repsi = sqrt(dens)*cos(theta)
impsi = sqrt(dens)*sin(theta)
#repsi = sqrt(dens)
#impsi = sqrt(dens)*0
# write out to new file
f1 = h5.File('./ic.enzo/GridRePsi', 'w')
new_dens = f1.create_dataset('GridRePsi', data=repsi)
for a in dataset.attrs.keys():
    new_dens.attrs.create(a, dataset.attrs[a])
f1.close()

f1 = h5.File('./ic.enzo/GridImPsi', 'w')
new_dens = f1.create_dataset('GridImPsi', data=impsi)
for a in dataset.attrs.keys():
    new_dens.attrs.create(a, dataset.attrs[a])
f1.close()
f.close()
