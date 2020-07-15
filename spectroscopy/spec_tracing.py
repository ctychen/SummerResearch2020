import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from pydl.pydlutils import bspline
from lmfit import Model
import ccdproc


def tracing_fit(pix_arr, flux_arr, peak_pix, sigma_init,peak_value):

    def gaussian(x, amp, cen, wid):

        return amp * np.exp(-(x-cen)**2 / (2*wid**2))


    gmodel = Model(gaussian)
    result = gmodel.fit(flux_arr, x=pix_arr, amp=peak_value, cen=peak_pix, wid=sigma_init)
    print(result.best_values)
    return result.best_fit, result.best_values


#load dome data:
prefix = r'C:\Users\claire\Documents\_ research _\notebooks\spec_reduction_lecture_student\data_for_class\\'
prefix_processed = r'C:\Users\claire\Documents\_ research _\notebooks\spec_reduction_lecture_student\processed_data\\'
obs_info = np.loadtxt(prefix+'obs_info.txt',dtype='str',delimiter= ",")

file_name = obs_info[:,0]
obstype = obs_info[:,1]
exptime = obs_info[:,2]
obj_name = obs_info[:,3]
x_size = obs_info[:,4]
y_size = obs_info[:,5]
color_name = obs_info[:,6]


#start to process bias_data
color_select = 'r'
stand_id = 'Feige34'
#start to process flat d
bias_data = fits.open(prefix_processed+'bias_'+color_select+'.fits')[0].data
flat_data = fits.open(prefix_processed+'flat_'+color_select+'.fits')[0].data
index_science = (obj_name == stand_id)  & (color_name == color_select)


print(file_name[index_science])

sigma =2.0
edge_l = 25
edge_r = 1185
edge_low = 55
edge_up = 230


for file_id in file_name[index_science]:

    science_data = fits.open(file_id)[0].data[edge_low:edge_up, edge_l:edge_r]
    science = (science_data - bias_data[edge_low:edge_up, edge_l:edge_r])/(flat_data+10**-18)

    plt.figure(figsize=(12,6))
    plt.subplot(1,2,1)

    plt.imshow(science, origin = 'lower', vmin = 0, vmax = 1000, aspect='auto', cmap="jet")

    science, mask = ccdproc.cosmicray_lacosmic(science, sigclip=5)
    plt.subplot(1,2,2)
    #plt.clf()
    plt.imshow(science, origin = 'lower', vmin = 0, vmax = 1000, aspect='auto', cmap="jet")
    plt.colorbar()
    plt.show()


    plt.figure(figsize=(12,6))
    plt.imshow(science, origin = 'lower', vmin = 0, vmax = 1000, aspect='auto', cmap="jet")
    plt.colorbar()
    plt.show()

    pix_array = np.arange(science.shape[0], dtype = np.float32)

    trace_pair = np.zeros((3, science.shape[1]))
    #130, 160, 142
    #80:200, 142, 130, 160 fit
    #bg 90:110, 180,190
    plt.plot(science[:,int(science.shape[1]/2)])
    plt.show()


    for idx in np.arange(science.shape[1]):

        #find peaks:
        #peak 87
        print(idx)
        # find max val in the fitting region to help get bounds
        idx_max = np.argmax(science[67:107,idx])
        bg = np.concatenate( (science[25:50,idx].reshape((-1,1)), science[125:250,idx].reshape((-1,1)) ), axis=0)

        mask_arr_bg = np.zeros(science.shape[0])
        mask_arr_bg[27:57] =1.0
        mask_arr_bg[117:147] =1.0
        y_bg = science[:,idx]*mask_arr_bg

        pix_arr_bg = np.arange(science.shape[0])
        index_bg  = y_bg > 0.0
        pars_bg = np.polyfit(pix_arr_bg[index_bg], y_bg[index_bg], 3)
        func_bg = np.poly1d(pars_bg)
        y_fit_bg = func_bg(pix_arr_bg)


        pix_max = pix_array[67:107][idx_max]
        peak_max = science[67:107,idx][idx_max]
        fit, values = tracing_fit(pix_array[57:117], science[57:117,idx] - y_fit_bg[57:117], pix_max, 2.0, peak_max)

        if idx>science.shape[1]/2*100:
            plt.plot(pix_array, science[:,idx])
            plt.plot(pix_array[57:117],fit+y_fit_bg[57:117])
            plt.axvline(values['cen'] - sigma*values['wid'],color='red')
            plt.axvline(values['cen'] + sigma*values['wid'],color='red')
            plt.show()


        trace_pair[0,idx] = values['cen'] - sigma*values['wid']
        trace_pair[1,idx] = values['cen'] + sigma*values['wid']
        trace_pair[2,idx] = values['cen']

        #print(trace_pair[0,:],trace_pair[1,:])
    #index_good = np.abs(trace_pair[0,10:1190] - np.median(trace_pair[0,10:1190]) ) < 10.0

    pars0 = np.polyfit(np.arange(science.shape[1]), trace_pair[0,:], 3)
    func0 = np.poly1d(pars0)


    plt.figure(figsize=(12,6))
    plt.subplot(1,2,1)
    plt.plot(np.arange(science.shape[1]), trace_pair[0,:])
    plt.plot(np.arange(science.shape[1]),func0(np.arange(science.shape[1])))
    plt.xlabel('Dispersion Direction')
    plt.ylabel('Left Boundary')

    pars1 = np.polyfit(np.arange(science.shape[1]), trace_pair[1,:], 3)
    func1 = np.poly1d(pars1)

    plt.subplot(1,2,2)
    plt.plot(np.arange(science.shape[1]), trace_pair[1,:])
    plt.plot(np.arange(science.shape[1]),func1(np.arange(science.shape[1])))
    plt.xlabel('Dispersion Direction')
    plt.ylabel('Right Boundary')
    plt.show()

    #plt.clf()
    plt.imshow(science, origin = 'lower', vmin = 0, vmax = 1000, aspect='auto', cmap="jet")
    plt.plot(np.arange(science.shape[1]), func0(np.arange(science.shape[1])),color='white')
    plt.plot(np.arange(science.shape[1]), func1(np.arange(science.shape[1])),color='white')
    plt.colorbar()
    plt.show()

    lt_arr = func0(np.arange(science.shape[1]))
    rt_arr = func1(np.arange(science.shape[1]))

    final_spec = np.zeros(science.shape[1])

    for kk in np.arange(science.shape[1]):

        pix_arr = np.arange(science.shape[0])

        lt= np.argmin( np.abs(lt_arr[kk]-pix_arr) )
        rt= np.argmin( np.abs(rt_arr[kk]-pix_arr) )

        mask_arr = np.zeros(science.shape[0])
        mask_arr[lt-20:lt] =1.0
        mask_arr[rt:rt+20] =1.0

        y = science[:,kk]*mask_arr

        pix_arr = np.arange(science.shape[0])
        index  = y > 0.0
        #plt.plot(pix_arr, y)
        #plt.show()
        #print(pix_arr.shape, y.size)
        pars = np.polyfit(pix_arr[index], y[index], 3)
        func = np.poly1d(pars)
        y_fit = func(pix_arr)
        if kk > 100*science.shape[1]/2:
            print(kk)
            plt.plot(pix_arr[lt-20:rt+20], science[lt-20:rt+20,kk], color='blue')
            plt.plot(pix_arr[index], y[index],color='red', marker = '.', linestyle='None')
            plt.plot(pix_arr[lt-20:rt+20], y_fit[lt-20:rt+20],color='black')
            plt.axvline(pix_arr[lt])
            plt.axvline(pix_arr[rt])
            plt.show()

        final_spec[kk] = np.sum((science[:,kk] - y_fit)[lt:rt])


    plt.plot(final_spec)
    plt.show()
    out_file = prefix_processed+stand_id+'_'+color_select+'.txt'
    np.savetxt(out_file, np.array([lt_arr, rt_arr, final_spec]).T)
