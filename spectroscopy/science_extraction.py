import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import ccdproc

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


color_select = 'r'
extract_id = 'J1624+6259'
#extract_id = 'J0819+4808'
stand_id = 'Feige34'
edge_l = 25
edge_r = 1185
edge_low = 55
edge_up = 230

bias_data = fits.open(prefix_processed+'bias_'+color_select+'.fits')[0].data
flat_data = fits.open(prefix_processed+'flat_'+color_select+'.fits')[0].data

index_science = (obj_name == extract_id) & (color_name == color_select)

trace_cube = np.loadtxt(prefix_processed+stand_id+'_'+color_select+'.txt')
lt_arr= trace_cube[:,0]
rt_arr= trace_cube[:,1]

print(file_name[index_science])


for idx, file_id in enumerate(file_name[index_science]):

    science_data = fits.open(file_id)[0].data[edge_low:edge_up, edge_l:edge_r]
    science = (science_data - bias_data[edge_low:edge_up, edge_l:edge_r])/(flat_data+10**-18)
    science, mask = ccdproc.cosmicray_lacosmic(science, sigclip=5)



    plt.figure(figsize=(12,6))
    plt.imshow(science, origin = 'lower', vmin = 0, vmax = 1000, aspect='auto', cmap="jet")
    plt.plot(np.arange(science.shape[1]), lt_arr,color='red')
    plt.plot(np.arange(science.shape[1]), rt_arr,color='red')
    plt.colorbar()
    plt.show()


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
        pars = np.polyfit(pix_arr[index], y[index], 3)
        func = np.poly1d(pars)
        y_fit = func(pix_arr)
        if kk > 100*science.shape[1]/2:

            plt.plot(pix_arr[lt-20:rt+20], science[lt-20:rt+20,kk], color='blue')
            plt.plot(pix_arr[index], y[index],color='red')
            plt.plot(pix_arr[lt-20:rt+20], y_fit[lt-20:rt+20],color='black')
            plt.show()

        final_spec[kk] = np.sum((science[:,kk] - y_fit)[lt:rt])


    plt.plot(final_spec)
    plt.show()

    out_file = 'science_spec_'+color_select+'_'+str(idx)+'.txt'
    np.savetxt(prefix_processed+out_file, np.array([final_spec]).T)
