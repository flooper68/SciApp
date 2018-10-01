import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as fft

start = 0
stop = 100
N = 20
dt = (stop-start)/N

time_axes = np.linspace(start, stop, N)

time_values = np.exp(1j*0.05*time_axes)#2 * np.exp(-time_axes**2/40) * np.exp(1j*500*time_axes)

waves = fft.fft(np.real(time_values))

T = dt * N
dw = 2*np.pi/T
frequencies = fft.fftfreq(N) * N * dw

def inverse_ft(k):
    val = 0.0
    for n in range(len(waves)):
        # val += waves[n] * np.exp(1.j * 2 * np.pi * n * k / len(waves)) / len(waves)
        print(time_axes[k])
        val += waves[n] * np.exp(1.j * frequencies[n] * time_axes[k]) / len(waves)
    return val.real

res = np.zeros(N)
for k in range(0, N):
    val = 0.0
    for n in range(0, len(waves)):
        # https://dsp.stackexchange.com/a/510/25943
        # val += waves[n]*np.exp(1.j * 2*np.pi * n * k / len(waves)) / len(waves)
        val += waves[n] * np.exp(1.j * frequencies[n] * k) / len(waves)
    res[k] = val.real


for x in range(len(waves)):
    if np.abs(waves)[x] == np.abs(waves.max()):
        print(frequencies[x])

plt.plot(frequencies, np.abs(waves))
plt.figure()
plt.plot(time_axes, res, 'b--', time_axes, np.real(time_values), 'r-', time_axes, fft.ifft(waves), 'g-.')
plt.show()
