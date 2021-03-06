# -*- coding: utf-8 -*-
"""PL1_advanced_notebook.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1jrRAJwAUQdxXiPfnm7zCsz9nU34NaUdf

# Plotting: Images
"""

from sklearn.datasets import load_sample_image
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = [8.0, 6.0]

i_flower = load_sample_image('flower.jpg')

for i in range(3):
  plt.figure()
  plt.imshow(i_flower[:, :, i])

i_flower = load_sample_image('flower.jpg')

for i in range(3):
  plt.figure()
  plt.imshow(i_flower[:, :, i])

i_flower.shape

np.shape(i_flower)

from PIL import Image
import requests
from io import BytesIO
import numpy as np

url = 'https://upload.wikimedia.org/wikipedia/commons/4/4f/Black_hole_-_Messier_87_crop_max_res.jpg'
response = requests.get(url, stream=True)
im = np.array(Image.open(response.raw))

plt.figure()
plt.imshow(im[:, :, i])

im.shape

import matplotlib.pyplot as plt
plt.figure()
plt.plot(im[2000, :, 0], '-r')
plt.plot(im[2000, :, 1] * 1, '-g')
plt.plot(im[2000, :, 2] * 1, '-b')

"""# Plotting: line & marker, colors & styles"""

x = np.arange(20)

y=[]
y.append(x ** 1.0)
y.append(x ** 1.1)
y.append(x ** 1.2)

plt.figure()
plt.plot(x, y[0], marker='.', linestyle='-', color='r')
plt.plot(x, y[1], '*--b')
plt.plot(x, y[2], '^-.g')

plt.figure()
plt.plot(x, y[1], 'b*-')

plt.plot(x, y[1], 'b*--')

plt.plot(x, y[0], marker='.', linestyle='-', color=[0.66, 0.5, 0.44])

"""# Plotting: Subplots"""

data_arr1 = np.arange(0, 110, 10) ** 2

data_arr2 = np.arange(0, 110, 10) * 0.5



plt.subplot(2,1,1)

plt.plot(data_arr1)


plt.subplot(2,1,2)

plt.plot(data_arr2)

data_arr1 = np.arange(0, 110, 10) ** 2

data_arr2 = np.arange(0, 110, 10) * 0.5



plt.subplot(1,2,1)

plt.plot(data_arr1)


plt.subplot(1,2,2)

plt.plot(data_arr2)

"""# Plotting: Legend, title, labels, & gridlines"""

#plt.style.use('dark_background')
plt.style.use('default')
plt.rcParams['figure.figsize'] = [8.0, 6.0]


data_arr = np.arange(0, 110, 10)

data_arr1 = data_arr ** 2

data_arr2 = data_arr * 0.5



plt.figure(figsize=(8, 8))

plt.subplot(2,1,1)

plt.title('Data squared')

plt.plot(data_arr1, '-r', label='Our data squared')

plt.plot(data_arr, '-', color=[0, 1, 0], label='Our original data')

plt.xlabel('X')

plt.ylabel('data$^2$')   # <---- The dollar sign "$" turns on latex interpreter for mathematical symbols/writing

plt.grid(True)

plt.legend()


plt.subplot(2,1,2)

plt.title('Data halved')

plt.plot(data_arr2, '-c', label='Our data divided by 2')

plt.plot(data_arr, color=[0, 1, 0], label='Our original data')

plt.xlabel('X')

plt.ylabel('data / 2')

plt.grid(True)

plt.legend()

"""# Fitting

Fitting is a general term that ususally means to obtain a mathematical decription of the data you are trying to fit. This mathematical description is typically some time of line if you are working with 2D data or a surface fit if its 3D data. With the typical types of fitting (polynomials and splines), you are granted the ability to predict. Prediciton can be done by extrapolation or interpolation. These two tools are discussed further in later sections.

**Polynomials**
"""

from numpy import polyfit, polyval
import numpy as np
import matplotlib.pyplot as plt

# Create some data
x = np.arange(100)

m = 2.3

b = 55

y = m * x + b

# Perform the 1-degree polynomial fit
pfc = np.polyfit(x, y, deg=1)

# Extract coefficients
m_solved = pfc[0]
b_solved = pfc[1]

# Evaluate the fit at some points
yfit = np.polyval(pfc, x)

# Plot all of this data
plt.figure(figsize=(10,10))
plt.plot(x, y, '.k')
plt.plot(x, yfit, '-r')

# Show the original & solved parameters
print(b, b_solved)
print(m, m_solved)



# Define initial parameters

x = np.arange(100)

noise = np.random.random(len(x)) * 15

m = 2.3

b = 55

y = m * x + b

y = y + noise

import matplotlib.pyplot as plt

# polyfit(), polyval()

plt.plot(x, y, '.k')

# Peform 1-degree polynomial fit

# .....


# Show solved coefficients vs original


# Plot fit and original data



time_myr = (np.arange(100)) # Days

time_myr = np.repeat(time_myr, 4)

noise = np.random.normal(0.0, 1450.5, len(time_myr))

mass = noise + 70.44 * time_myr**1.5 + 2.7

mass = mass + 4000

plt.plot(time_myr, mass, '.b', label='Measured data')

plt.xlabel('Time [Myr]')
plt.ylabel('Mass [M$_{sun}$]')
plt.legend()
plt.title('Mass of a black hole over time')

pfc = np.polyfit(time_myr, mass, deg=7)

time_myr2 = np.linspace(time_myr[0], time_myr[-1] + 100, 50)

plt.plot(time_myr, mass, '.k', label='Measured data')
plt.plot(time_myr2, np.polyval(pfc, time_myr2), '-r')

time_fit = np.linspace(time_myr[0], time_myr[-1], 100)


p1_coef = np.polyfit(time_myr, mass, 1)

p3_coef = np.polyfit(time_myr, mass, 3)


mass_p1 = np.polyval(p1_coef, time_fit)

mass_p3 = np.polyval(p3_coef, time_fit)


plt.plot(time_fit, mass_p1, '--r', label='Poly Deg 1 Fit', linewidth=3)
plt.plot(time_fit, mass_p3, '-', color=[0,1,0], label='Poly Deg 3 Fit', linewidth=5)

plt.plot(time_myr, mass, '.b', label='Measured data')


plt.xlabel('Time [Myr]')
plt.ylabel('Mass [M$_{sun}$]')
plt.legend()
plt.title('Mass of a black hole over time')

"""**Splines**"""

from scipy.interpolate import interp1d

# Create some data
x = np.arange(50)

y = np.sin(2*np.pi/25 * x)

# Perform the cubic spline fit
yf_pp = interp1d(x, y, kind='cubic')

# Evaluate the fit at some points
yf = yf_pp(x)

# Plot the original data and the fit
plt.figure(figsize=(10,6))
plt.plot(x, y, '.-k')
plt.plot(x, yf, '-r')

# Things to consider: Spline will not deal well with nonsmooth aor sparse data

# Create some data
x = np.arange(50)

y = np.sin(2*np.pi/30 * x)

y[y > 0.25] = 0.25

# Perform the cubic spline fit
yf_pp = interp1d(x[::5], y[::5], kind='cubic', fill_value='extrapolate')

# Evaluate the fit at some points
yf = yf_pp(x)

# Plot the original data and the fit
plt.figure(figsize=(10,6))
plt.plot(x, y, '.-k')
plt.plot(x, yf, '-r')
plt.plot(x[::5], y[::5], '*g', markersize=25)



"""**Smoothing Splines**"""

!pip install csaps

import csaps

from csaps import UnivariateCubicSmoothingSpline as ss

x = np.arange(50)

y = np.sin(2*np.pi/10 * x)


# pfc = np.polyfit(x, y, deg=14)

x2 = np.linspace(x[0], x[-1], 500)

## ss is for smoothing splines, range from 0-1 for smooth param
yf_pp = ss(x, y, smooth=0.8)

y2 = yf_pp(x2)


plt.figure(figsize=(10,6))
plt.plot(x, y, '.-k')
plt.plot(x2, y2, '-r')
plt.plot(x2, ss(x, y, smooth=0.99)(x2), '-')
plt.plot(x2, ss(x, y, smooth=0.5)(x2), '-')
plt.plot(x2, ss(x, y, smooth=0.1)(x2), '-m')

x = np.arange(100)

noise = np.random.random(len(x)) * 35

m = 2.3

b = 55

y = m * x + b

y = y + noise

x2 = np.linspace(x[0], x[-1], 500)

plt.figure(figsize=(15,9))
plt.plot(x, y, '.k')

# for s in [0.99, 0.5, 0.1, 1e-5]:
plt.plot(x2, ss(x, y, smooth=1e-5)(x2), '-', label=str(s))

plt.legend()

x = np.arange(100)

noise = np.random.random(len(x)) * 35

m = 2.3

b = 55

y = m * x + b

y = y + noise

x2 = np.linspace(x[0], x[-1], 500)

plt.figure(figsize=(15,9))
plt.plot(x, y, '.k')

for s in [0.99, 0.5, 0.1, 1e-5]:
  plt.plot(x2, ss(x, y, smooth=s)(x2), '-', label=str(s))

plt.legend()



"""**Extrapolation and Interpolation**

Typically, all functions that allow us to fit data will also alow us to extrapolate and interpolation, but not always. Also, some functions reqire us to specify if extrapolation is allowed.

Spline Extap/Interp
"""

time_sec = np.linspace(0, 10, num=11, endpoint=True)

velocity = np.cos(-time_sec**2/9.0)

print('Time = ', time_sec)
print('velocity = ', velocity)

plt.figure(figsize=(15,9))
plt.plot(time_sec, velocity, 'or', label='Measured data')


plt.plot([3.5, 3.5], [-1, 1], '--g', label='Missing data we want to know')
plt.plot([6.4, 6.4], [-1, 1], '--g')
plt.plot([8.6, 8.6], [-1, 1], '--g')

plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.legend()
plt.title('Velocity of a ball over time')

from scipy.interpolate import interp1d


cubic_interpolant = interp1d(time_sec, velocity, kind='cubic')


times_i_want_to_know = np.array([3.5, 6.4, 8.6])

velocity_approx = cubic_interpolant( times_i_want_to_know )


plt.figure(figsize=(15,9))
plt.plot(time_sec, velocity, 'or', label='Measured data')

plt.plot([3.5, 3.5], [-1, 1], '--g', label='Missing data we want to know')
plt.plot([6.4, 6.4], [-1, 1], '--g')
plt.plot([8.6, 8.6], [-1, 1], '--g')

plt.plot(times_i_want_to_know, velocity_approx, '*y', label='Approximated Velocities', markersize=15)

plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.legend()
plt.title('Velocity of a ball over time')

cubic_interpolant = interp1d(time_sec, velocity, kind='cubic')


times_i_want_to_know = np.arange(0, 10, 0.01)

velocity_approx = cubic_interpolant( times_i_want_to_know )


plt.figure(figsize=(15,9))
plt.plot(time_sec, velocity, 'or', label='Measured data')

plt.plot(times_i_want_to_know, velocity_approx, '-y', label='Approximated Velocities')

plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.legend()
plt.title('Velocity of a ball over time')

cubic_interpolant = interp1d(time_sec, velocity, kind='cubic', fill_value='extrapolate')


times_i_want_to_know = np.arange(0, 10.0, 0.01)

velocity_approx = cubic_interpolant( times_i_want_to_know )


plt.figure(figsize=(15,9))
plt.plot(time_sec, velocity, 'or', label='Measured data')

plt.plot(times_i_want_to_know, velocity_approx, '-y', label='Approximated Velocities')

for tf in np.linspace(10, 10.8, 8):
  plt.plot(tf, cubic_interpolant( tf ), '.')

plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.legend()
plt.title('Velocity of a ball over time')

"""Polynomial Extrap/interp"""

time_myr = (np.arange(100)) # Days

time_myr = np.repeat(time_myr, 4)

noise = np.random.normal(0.0, 1450.5, len(time_myr))

mass = noise + 70.44 * time_myr**1.5 + 2.7

mass = mass + 4000

plt.plot(time_myr, mass, '.b', label='Measured data')

plt.xlabel('Time [Myr]')
plt.ylabel('Mass [M$_{sun}$]')
plt.legend()
plt.title('Mass of a black hole over time')

pfc = np.polyfit(time_myr, mass, deg=2)

time_myr2 = np.linspace(time_myr[0], time_myr[-1] + 100, 50)

plt.plot(time_myr, mass, '.k', label='Measured data')
plt.plot(time_myr2, np.polyval(pfc, time_myr2), '-r')

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 50, 100)

pha = 5
per = 15

y = 4 * np.sin(pha + 2*np.pi / per * x)

plt.plot(x, y, '.-k')

#count inflections for degrees of freedom

pfc = np.polyfit(x, y, 9)

y2 = np.polyval(pfc, x)
plt.plot(x, y, '.-k')
plt.plot(x, y2, '.-r')

np.random.normal(size=len(y))

yn = y + np.random.normal(size=len(y))*10

plt.plot(x, y, '.-k')
plt.plot(x, yn, '.-r')

t = np.arange(100)

d = np.exp(t/10)

plt.plot(t,d)

from scipy.interpolate import interp1d

a = interp1d(t, d, kind='cubic')

# evaluates at same time as original data
df = a(t)
plt.plot(t,d, '.-r', linewidth = 10)
plt.plot(t, df, '_-b')
