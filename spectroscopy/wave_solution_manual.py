import numpy as np
import matplotlib.pyplot as plt


prefix_processed = r'C:\Users\claire\Documents\_ research _\notebooks\spec_reduction_lecture_student\processed_data\\'
spectra = []

for idx in np.arange(3):

    out_file = prefix_processed+'arc_spec_r_'+str(idx)+'.txt'
    spectra.append(np.loadtxt(out_file))


spectrum = np.mean(np.array(spectra),axis=0)
template = np.loadtxt('kast_template_r.txt')
plt.figure(figsize=(18,6))
plt.subplot(2,1,1)
plt.plot(np.arange(spectrum.size), spectrum)
plt.subplot(2,1,2)
plt.plot(template[:,0], template[:,1])
plt.show()

# need to find these manually from plot, find the values tha correspond
# have to find this for around 10-15 lines to fit well
line_pix = np.array([113, 562, 292, 406, 438, 459, 783, 955])
line_wav = np.array([5461, 6506.5, 5875.5, 6143, 6217.8, 6267.1, 7031.8, 7438.3])

# line_pix = array([113., 247., 256., 282., 292., 322., 335., 358., 377., 386., 406.,
#        415., 438., 459., 475., 488., 508., 517., 561., 572., 600., 633.,
#        650., 739., 783., 797., 843., 873., 955.])
#
# line_wav =np.array([ 5461.1, 5769.2, 5790.3, 5851.9,  5875.4,\
#  5945.8,   5976.3, 6030.1, 6074.5, 6095.6,\
#  6142.4, 6163.5, 6217.3,  6266.6, 6304.1,\
#  6334.7, 6381.7, 6402.9, 6506.7, 6532.7,\
#  6599.0, 6677.2, 6717.5, 6928.7, 7033.0,\
#  7066.1, 7175.0, 7245.6, 7438.0])



#
'''
for kk in np.arange(line_wav.size):

    print(line_wav[kk], line_pix[kk])
    plt.figure(figsize=(18,6))
    plt.subplot(2,1,1)
    plt.plot(np.arange(spectrum.size), spectrum)
    plt.axvline(line_pix[kk],color='red')
    plt.subplot(2,1,2)
    plt.plot(template[:,0], template[:,1])
    plt.axvline(line_wav[kk], color='blue')
    plt.show()
'''

polydeg = 3
best_p = np.polyfit(line_pix, line_wav, polydeg)
err = np.polyval(best_p, line_pix) - line_wav
print(best_p)
labels = 'RMS = %0.3f' % np.std(err)
plt.figure(figsize=(12,6))
plt.subplot(1,2,1)
plt.plot(line_pix, line_wav, color='red', marker='.', linestyle='None')
plt.plot(line_pix, np.polyval(best_p, line_pix), color='black')
plt.xlabel('Pixel')
plt.ylabel('Wavelength')
plt.subplot(1,2,2)
plt.plot(line_pix, err, color='red', marker='.', linestyle='None',label = labels)
plt.xlabel('Dispersion direction')
plt.ylabel('Residuals')
plt.legend(loc='best')
plt.show()

wav_ma = np.polyval(best_p, np.arange(spectrum.size))
plt.plot(wav_ma, spectrum,color = 'red',label='Lick Data')
plt.plot(template[:,0], template[:,1], color = 'black', label='Template')
plt.xlabel('wavelength')
plt.ylabel('Counts')
plt.legend(loc='best')
plt.show()

np.savetxt(prefix_processed+'solution_manual_r.txt',best_p.T)
