from astropy.io import fits
import matplotlib.pyplot as plt


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
index_bias = (obj_name == 'bias')  & (color_name == color_select) & (x_size == '250') & (y_size == '1232')
bias_cube = []

print(file_name[index_bias])
print('There are %s bias data' % file_name[index_bias].size)
for file_id in file_name[index_bias]:
    bias_cube.append(fits.open(file_id)[0].data)

bias_frame = np.median(np.asarray(bias_cube), axis=0)
out_file = prefix_processed+'bias_'+color_select+'.fits'
hdu = fits.PrimaryHDU(bias_frame)
hdu.writeto(out_file,overwrite=True)
print('Done')
