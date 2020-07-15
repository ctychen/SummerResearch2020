import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import ccdproc
from pydl.pydlutils import bspline
from astropy import table
from scipy import interpolate
from astropy import units
from calib_utils import *

st=fits.open(r'C:\Users\claire\Documents\_ research _\notebooks\spec_reduction_lecture_student\source\feige34_005.fits.gz')[1].data
wav_st = st['WAVELENGTH']
flux_st = st['FLUX']
exptime = 240
airmass = 1.03

aa = np.loadtxt(r'C:\Users\claire\Documents\_ research _\notebooks\spec_reduction_lecture_student\source\J1624+6259.txt')

spectra = []

prefix_processed = r'C:\Users\claire\Documents\_ research _\notebooks\spec_reduction_lecture_student\processed_data\\'
color_select = 'r'
stand_id = 'Feige34'


st_spec = np.loadtxt(prefix_processed+stand_id+'_'+color_select+'.txt')[:,2]
solution  = np.loadtxt(prefix_processed+'solution_manual_r.txt')

pix = np.arange(st_spec.size)
wave = np.polyval(solution.T, pix)


sensfunc_obs, ext_corr = sensfunc(wave, st_spec, wav_st, flux_st, exptime, airmass)
zem =1.72#1.994 #1.72
for idx in np.arange(2):


    sci_spec  = np.loadtxt(prefix_processed+'science_spec_'+color_select+'_'+str(idx)+'.txt')
    senstot = sensfunc_obs * ext_corr
    flam = sci_spec * senstot / 1800.0

    plt.plot(wave/(1+zem),flam*10**17)
    plt.show()

    spectra.append(flam*10**16)
print('Done')
plt.plot(wave/(1+zem), np.mean(np.array(spectra),axis=0),label = 'Lick')
plt.plot(aa[:,0]/(1+zem), aa[:,1],label = 'J1624+6259')
plt.xlabel('Wavelength')
plt.ylabel(r'Flux ($10^{-17} ergs/s/cm^2/\AA$)')
plt.legend(loc='best')
plt.show()
exit()
#r$\mathrm{\AA}$
