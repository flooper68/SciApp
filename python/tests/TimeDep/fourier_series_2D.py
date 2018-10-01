"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

# First set up the figure, the axis, and the plot element we want to animate

delta = 0.1
x = y = np.arange(-3.0, 3.0, delta)
X, Y = np.meshgrid(x, y)

x_2 = y_2 = range(len(x))
X_2, Y_2 = np.meshgrid(x_2, y_2)

def myFun(x,y):
    x_0 = 0
    y_0 = 0
    # return np.real(np.exp(1j*(0.5*x + 0)))
    return np.exp( -((x-x_0)**2 + (y-y_0)**2)/0.3)

Z = np.vectorize(myFun)(X,Y)

ft = np.fft.fftshift(np.fft.fft2(Z))
freq = np.fft.fftshift(np.fft.fftfreq(Z.shape[-1]))

fig = plt.figure()
im = plt.imshow(Z, cmap=cm.coolwarm,
               origin='lower', extent=[-3, 3, -3, 3])

fig2 = plt.figure()
im2 = plt.imshow(np.abs(ft), cmap=cm.coolwarm, extent=[freq[0], freq[-1], freq[0], freq[-1]]
                 ,vmax=abs(ft).max(), vmin=-abs(ft).max())

fig3 = plt.figure()
im3 = plt.imshow(np.angle(ft), cmap=cm.coolwarm, extent=[freq[0], freq[-1], freq[0], freq[-1]])

ft = np.fft.fft2(Z)


def myFun2(x,y):
    sum = 0
    for i in range(freq.shape[0]):
        for j in range(freq.shape[0]):
            tosum = ft[j, i] * np.exp(1j*2*np.pi*(j*y/freq.shape[0] + i*x/freq.shape[0]))/(freq.shape[0])**2

            sum += tosum
    return np.real(sum)


Z_fft= np.vectorize(myFun2)(X_2, Y_2)

fig5 = plt.figure()
im5 = plt.imshow(np.real(Z_fft),cmap=cm.coolwarm,
               origin='lower', extent=[-3, 3, -3, 3])

plt.show()

