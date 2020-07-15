import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from pydl.pydlutils import bspline
from astropy import table
from scipy import interpolate
from astropy import units


def extinction_corr(wave, airmass):


    extinct = table.Table.read('./source/lick_extin.txt', comment='#', format='ascii',
                               names=('iwave', 'mag_ext'))
    wave_tab = table.Column(np.array(extinct['iwave']) * units.AA, name='wave')
    extinct.add_column(wave_tab)
    extinct_final = extinct[['wave', 'mag_ext']]

    f_mag_ext = interpolate.interp1d(extinct['wave'], extinct['mag_ext'], bounds_error=False,
                                     fill_value=0.)
    mag_ext = f_mag_ext(wave)

    flux_corr = 10.0 ** (0.4 * mag_ext * airmass)

    return flux_corr


def standard_sensfunc(wave, flux, ivar, flux_std):

    # Create copy of the arrays to avoid modification
    wave_obs = wave.copy()
    flux_obs = flux.copy()
    ivar_obs = ivar.copy()


    # Removing outliers
    logflux_obs = 2.5 * np.log10(flux_obs)
    logflux_std = 2.5 * np.log10(flux_std)

    logflux_obs_test = logflux_obs[np.isfinite(logflux_obs) ]

    plt.plot(wave_obs,logflux_std-logflux_std.min(), color = 'black', label='Template')
    plt.plot(wave_obs[np.isfinite(logflux_obs) ],logflux_obs_test-logflux_obs_test.min(), color = 'red', label='Data')
    plt.xlabel('wavelength')
    plt.ylabel('flux')
    plt.legend(loc='best')    
    plt.show()


    magfunc = logflux_std - logflux_obs

    msk_fit_sens = np.isfinite(ivar_obs) & np.isfinite(logflux_obs) & np.isfinite(logflux_std) #& (~index_line)

    #plt.plot(wave_obs[msk_fit_sens], magfunc[msk_fit_sens])
    #plt.xlabel('wavelength')
    #plt.ylabel('ratio')
    #plt.show()

    curve = bspline.iterfit(wave_obs[msk_fit_sens], magfunc[msk_fit_sens], nord=5,bkspace=50,maxiter=2,low=0.2,hi=0.01)[0]
    spline = curve.value(wave_obs)[0]

    sensfunc = 10.0 ** (0.4 * magfunc)

    plt.plot(wave_obs[msk_fit_sens], sensfunc[msk_fit_sens], label= 'Data')
    plt.plot(wave_obs,  10.0 ** (0.4 * spline),color='red',label='Fit')
    plt.xlabel('wavelength')
    plt.legend(loc='best')
    plt.show()

    return 10.0 ** (0.4 * spline)


def sensfunc(wave, spec, template_wav, template_flux, exptime, airmass):
    '''
    this function fits a spline funtion to the 2.5*log10(flux_std/flux_count)
    '''
    flux_star = spec/exptime
    ivar = 1/spec*exptime**2

    #extinction corr:
    ext_corr = extinction_corr(wave * units.AA, airmass)

    flux_star  = flux_star*ext_corr
    ivar = ivar / ext_corr**2

    # Interpolate the standard star onto the current set of observed wavelengths
    flux_true = interpolate.interp1d(template_wav, template_flux, bounds_error=False,
                                     fill_value='extrapolate')(wave)

    sensfunc= standard_sensfunc(wave, flux_star, ivar, flux_true)
    return sensfunc,ext_corr
