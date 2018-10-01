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
from matplotlib import animation

# First set up the figure, the axis, and the plot element we want to animate

delta = 0.1
x = y = np.arange(-3.0, 3.0, delta)
X, Y = np.meshgrid(x, y)
Z = 10*np.real(np.exp(-1j*X*2))

fig = plt.figure()
im = plt.imshow(Z, interpolation='bilinear', cmap=cm.coolwarm,
               origin='lower', extent=[-3, 3, -3, 3],
               vmax=abs(Z).max(), vmin=-abs(Z).max())

angle = -40
k_0 = 10
n = 4

k_y = k_0*np.sin(np.deg2rad(angle))
k_z1 = k_0*np.cos(np.deg2rad(angle))
k_z2 = np.sqrt((k_0*n)**2 - k_y**2)

frames = 100
list = []

amp = 10
r_coef = 0.5
t_coef = 0.5

ref_angle = np.rad2deg(np.arcsin(np.sin(np.deg2rad(angle))/n))

print(ref_angle)

def myFun(x,y,i):
    if y > 0:
        if y < -1/np.sin(np.deg2rad(-angle)) * (x - 0.5) and y < 1/np.sin(np.deg2rad(-angle))* (x + 1.5):
            return amp * np.real(np.exp(1j * (k_y * x + k_z1 * y + 0.1 * i))) + amp * r_coef * np.real(np.exp(1j * (k_y * x - k_z1 * y + 0.1 * i)))
        if y < -1/np.sin(np.deg2rad(-angle)) * (x - 0.5):
            return amp * np.real(np.exp(1j * (k_y * x + k_z1 * y + 0.1 * i))) + amp * np.real(np.exp(1j * (k_y * x + k_z1 * y + 0.2 * i)))+ amp * np.real(np.exp(1j * (k_y * x + k_z1 * y + 0.3 * i)))
        elif y < 1/np.sin(np.deg2rad(-angle)) * (x + 1.5):
            return amp * r_coef * np.real(np.exp(1j * (k_y * x - k_z1 * y + 0.1 * i)))
        else:
            return 0
    elif y < 0 and y < -1/np.sin(np.deg2rad(-ref_angle)) * (x - 0.5):
        return amp * t_coef * np.real(np.exp(1j*(k_y * x - k_z2 * y + 0.1*i)))
    else:
        return 0

for i in range(frames):
    list.append(np.vectorize(myFun)(X,Y,i))

# initialization function: plot the background of each frame
def init():
    im.set_data(Z)
    return im,

# animation function.  This is called sequentially
def animate(i):
    im.set_data(list[i])
    return im,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=frames, interval=20, blit=True)


plt.show()

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
# anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

