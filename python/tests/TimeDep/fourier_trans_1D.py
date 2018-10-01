import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import *

dt = 1
N = 300
T = dt * N
dw = 2*np.pi/T
frequencies = fftfreq(N) * N * dw

omega = 0.03
data = np.zeros(N)
points = np.linspace(-100, 100, N)
for x in range(N):
    data[x] = 2 * np.exp(-points[x]**2/20) * np.exp(1j*2*points[x]) #np.exp(1j*omega*x)

# apply FFT on samples
waves = fft(data)

res = np.zeros(N)
for k in range(0, N):
    val = 0.0
    for n in range(0, len(waves)):
        # val += waves[n]*np.exp(1.j * 2*np.pi * n * k / len(waves)) / len(waves)
        val += waves[n] * np.exp(1.j * frequencies[n] * k) / len(waves)
    res[k] = val.real

#np implementation
res2 = np.fft.ifft(waves)

plt.figure()
plt.plot(points, data, 'b') # original
plt.plot(points, res,'r--') # my implementation
plt.plot(points, res2,'c-.') # np implementation


frames = 200

list = []
for x in range(frames):
    res = np.zeros(N)
    for k in range(0, N):
        val = 0.0
        for n in range(0, len(waves)):
            # val += waves[n]*np.exp(1.j * 2*np.pi * n * k / len(waves)) / len(waves)
            val += waves[n] * np.exp(1.j * frequencies[n] * (x - k)) / len(waves)
        res[k] = val.real
    list.append(res)

fig = plt.figure()
ax = plt.axes()
line, = ax.plot(points, res)



def animate(i):
    line.set_data(points, list[i])
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, frames=frames, interval=100, blit=True)
plt.show()

