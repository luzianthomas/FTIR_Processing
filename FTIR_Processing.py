## simple processing of time domain infrared signals recorded on a FTIR spectrometer

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

filename1 = 'Irdata/BAlk.RIFG.dpt'
filename2 = 'Irdata/BAlk.SIFG.dpt'
dataBG = np.loadtxt(filename1, delimiter=',',skiprows=0)   # background time domain signal
dataS = np.loadtxt(filename2, delimiter=',',skiprows=0)    # sample time domain signal


## raw BG data plotted (mirror goes once in both directions -> ZPD crossed twice -> effectively 2 spectra taken)
fig = plt.figure()
plt.plot(dataBG[:,0],dataBG[:,1])
plt.xlabel('time')
plt.ylabel('signal')
plt.show()
fig.savefig('FTIR_timedomain.pdf')


## maximum taken (one of the ZPD) and N points around it plotted
N = pow(2,14)
ix = np.argmax(dataBG[:,1])
lowi = int(ix-N/2)
upperi = int(ix+N/2-1)
IFGR = dataBG[lowi:upperi,1]
IFGS = dataS[lowi:upperi,1]

fig = plt.figure()
plt.plot(IFGR)
plt.plot(IFGS)
plt.xlabel('index')
plt.ylabel('signal')
ax = plt.gca()
ax.set_xlim([8000,8400])
ax.set_ylim([-0.1,0.1])
plt.show()


## fft of N points around ZPD, absolute of resulting array is symmetric because it is a real signal (only half taken)
## middle point is spectrometer specific frequency (16707.63 cm^-1 in this case)
FT_BG = abs(np.fft.fft(IFGR,n=N))[1:int(N/2)]
FT_S = abs(np.fft.fft(IFGS,n=N))[1:int(N/2)]

length = len(FT_BG)
nu = np.linspace(start=0,stop=16707.63,num=length)

## Intensity spectra bg and sample overlayed
fig = plt.figure()
plt.plot(nu,FT_BG)
plt.plot(nu,FT_S)
plt.xlabel('wavenumber (cm-1)')
plt.ylabel('Intensity I')
ax = plt.gca()
ax.set_xlim([4000,500])
#ax.set_ylim([-0.1,0.1])
plt.show()
fig.savefig('FTIR_intensity.pdf')


## transmission spectrum
transm = FT_S/FT_BG


fig = plt.figure()
plt.plot(nu,100*transm)
plt.xlabel('wavenumber (cm-1)')
plt.ylabel('%Transmission')
ax = plt.gca()
ax.set_xlim([4000,500])
ax.set_ylim([0,100])
plt.show()
fig.savefig('FTIR_transmission.pdf')



## absorbance spectrum

absorb = -np.log10(transm)

## region for peak picking
nu_region = nu[np.argmin(abs(nu-500)):np.argmin(abs(nu-4000))]
absorb_region = absorb[np.argmin(abs(nu-500)):np.argmin(abs(nu-4000))]
## peak picking based on prominence
peaks_idx = find_peaks(absorb_region, prominence=max(absorb_region)*0.1)

fig = plt.figure()
plt.plot(nu_region,absorb_region)
plt.plot(nu_region[peaks_idx[0]], absorb_region[peaks_idx[0]], 'x', color="red")
plt.xlabel('Wavenumber (cm-1)')
plt.ylabel('Absorbance A')
ax = plt.gca()
ax.set_xlim([4000,500])
ax.set_ylim([-0.02,1])
plt.show()
fig.savefig('FTIR_absorbance.pdf')

## peak list
print(f'peak list:\nnu [cm-1]\tAbsorbance [A]')
for d in reversed(peaks_idx[0]):
    print(f'{nu_region[d]:.5f}\t{absorb_region[d]:.2f}')




