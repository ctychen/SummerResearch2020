import numpy as np
import os
from astropy.io import fits


prefix = r'C:\Users\claire\Documents\_ research _\notebooks\spec_reduction_lecture_student\data_for_class\\'

obstype = []
exptime = []
obj_name = []
x_size = []
y_size = []
color_name = []
file_array = []
type_arr = []

print(os.listdir(prefix))

for file_name in os.listdir(prefix):
    if file_name.endswith(".fits"):

        file_array.append(prefix+file_name)

        a=fits.open(prefix+file_name, ignore_missing_end=True)

        try:
            print('this is exp:'+str(a[0].header['EXPTIME']))
            exptime.append(a[0].header['EXPTIME'])
            obj_name.append(a[0].header['OBJECT'])
            color_name.append(file_name[0])
        except:
            print('not found')
            exptime.append(-999.0)
            obj_name.append('non_file')
            color_name.append(file_name[0])

        try:
            print('this is obj:'+a[0].header['OBSTYPE'])
            obstype.append(a[0].header['OBSTYPE'])
        except:
            print('not found')
            obstype.append('NAN')

        try:
            print(a[0].header['NAXIS1'],a[0].header['NAXIS2'])
            x_size.append(a[0].header['NAXIS2'])
            y_size.append(a[0].header['NAXIS1'])
        except:
            print('not found')
            x_size.append(-999.0)
            y_size.append(-999.0)

results = np.array([file_array,obstype, exptime, obj_name, x_size, y_size, color_name])
np.savetxt(prefix+'obs_info.txt', results.T, delimiter= ",", fmt="%s")
