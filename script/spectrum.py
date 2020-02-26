import matplotlib
matplotlib.use('pdf')
import h5py as h5
import numpy as np
from scipy import *
import pylab as plt
# Read Baryon density
f1 = h5.File('./ic.enzo/GridRePsi', 'r')
dataset = f1['GridRePsi']
re = np.array(dataset)

f2 = h5.File('./ic.enzo/GridImPsi', 'r')
dataset = f2['GridImPsi']
im = np.array(dataset)
f1.close()
f2.close()

n = re.shape[0]
print n
dens = re**2+im**2

print dens.mean()
print 'real mean',re.mean()
print np.mean(re**2)
print 'imag mean',im.mean()
print np.mean(im**2)

fd = np.fft.fftn(dens)/n**3
fre = np.fft.fftn(re/re.mean())/n**3
fim = np.fft.fftn(im/im.mean())/n**3 

nx, ny, nz = dens.shape
kx = np.fft.fftfreq(nx)*nx
ky = np.fft.fftfreq(ny)*ny
kz = np.fft.fftfreq(nz)*nz

kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
k = sqrt(kx3d**2 + ky3d**2 + kz3d**2)*2*pi/20
kbins=np.logspace(-1, 3, 100)
Pmass, bin_edges = np.histogram(k, weights=abs(fd)**2,  bins=kbins)
Pmassr, bin_edges = np.histogram(k, weights=abs(fre)**2, bins=kbins)
Pmassi, bin_edges = np.histogram(k, weights=abs(fim)**2, bins=kbins)
bincenter = 0.5*(bin_edges[1:]+bin_edges[:-1])
kwidth = (bin_edges[1:] - bin_edges[0:-1])
Pmass /= 4*pi*bincenter**2*kwidth
Pmassr /= 4*pi*bincenter**2*kwidth
Pmassi /= 4*pi*bincenter**2*kwidth
#plt.imshow (re[256,:,:], interpolation='nearest')
plt.imshow (arctan(im[256,:,:]/re[256,:,:]), interpolation='nearest')
plt.colorbar()
k3 = bincenter**3/(2*pi**2)

#plt.loglog(bincenter,Pmass*k3,label='Density')
#plt.loglog(bincenter,Pmassr*k3,label=r'$Re(\Psi)$')
#plt.loglog(bincenter,Pmassi*k3,label=r'$Im(\Psi)$')
#plt.ylim([1e-8,1e0])
#plt.xlabel(r'$k$')
#plt.ylabel(r'$k^3 P(k)$')
#plt.legend(loc='best')

plt.savefig('img.pdf')
