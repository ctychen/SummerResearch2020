import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from pydl.pydlutils import bspline

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
index_dome1 = (np.char.find(obj_name, 'dome') != -1) & (color_name == color_select)  \
              & (x_size == '250') & (y_size == '1232')

# define trimmed region
edge_l = 25
edge_r = 1185
edge_low = 55
edge_up = 230

#import bias data
bias_data = fits.open(prefix_processed+'bias_'+color_select+'.fits')[0].data
dome1_cube = []
print(file_name[index_dome1])
for file_id in file_name[index_dome1]:
    science_data = fits.open(file_id)[0].data
    science_raw = science_data - bias_data
    dome1_cube.append(science_raw)
    #plt.clf()
    #plt.imshow(science_raw, origin = 'lower', vmin = 0, vmax = science_raw.max())
    #plt.colorbar()
    #plt.show()
#does the median (or mean) combine
dome1_frame = np.median(np.asarray(dome1_cube), axis=0)
dome1_file = prefix_processed+'dome_master_'+color_select+'.fits'
hdu = fits.PrimaryHDU(dome1_frame)
hdu.writeto(dome1_file,overwrite=True)
print('Done')




print('Start to normalize the master flat.')
dome_data = fits.open(prefix_processed+'dome_master_'+color_select+'.fits')[0].data
#trim the data:
dome_data = dome_data[edge_low:edge_up, edge_l:edge_r]
flat_field =np.zeros(dome_data.shape)
for x_id in np.arange(dome_data.shape[0]):
    #print(x_id)
    pixel_arr = np.arange(dome_data.shape[1],dtype=np.float32)
    flux_count = dome_data[x_id,:].reshape(-1)


    curve = bspline.iterfit(pixel_arr, flux_count, nord=5,bkspace=50,maxiter=2,low=0.2,hi=0.01)[0]

# use spline to fit lamp spectrum data so it can be removed
    spline = curve.value(pixel_arr)[0]
    flat_field[x_id,:] = flux_count/spline

    # if x_id > 10000.0:
    if x_id == 100.0: # shows only 1 row (row #100)
        plt.plot(pixel_arr, flux_count,color='black')
        plt.plot(pixel_arr, spline,color='red')
        plt.show()

        plt.plot(pixel_arr, flux_count/spline,color='red')
        plt.show()


plt.imshow(flat_field, origin = 'lower', vmin = 0, vmax = 1.5, aspect='auto', cmap="jet")
plt.colorbar()
plt.show()
flat_file = prefix_processed+'flat_'+color_select+'.fits'
hdu = fits.PrimaryHDU(flat_field)
hdu.writeto(flat_file,overwrite=True)
