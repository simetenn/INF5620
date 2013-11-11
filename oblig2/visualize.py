import pylab as p
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
from matplotlib import cm
from scitools.std import movie
import glob, os, sys

import time

path = "./plots/"

#Remove old plotfiles
for filename in glob.glob(path+"wave*.png"):
    os.remove(filename)

#Load arrays
u = p.np.load("u.npy")
h = p.np.load("h.npy")
y = p.np.load("y.npy")
x = p.np.load("x.npy")

T = p.shape(u)[2]

X, Y = p.meshgrid(x, y)

#Plot for each timestep
for t in xrange(T):
    fig = p.figure()
    ax = fig.add_subplot(111, projection='3d')
    #Plot the wave
    surf = ax.plot_surface(X, Y, u[:,:,t], rstride=1, cstride=1, cmap=cm.get_cmap("winter_r"),#color="blue",
                            linewidth=0, antialiased=True, shade=True, 
                            vmin=p.amin(u),vmax=p.amax(u))
    #Plot the bottom
    surf2 = ax.plot_surface(X, Y, h, rstride=1, cstride=1, color="brown",
                                linewidth=0, antialiased=True, shade=True,
                                alpha=0.6)
    ax.set_zlim3d(p.amin(h), p.amax(u))
    fig.colorbar(surf)
    p.savefig(path + "wave%05d"  % (t) + ".png")
    

#Create the movie
movie(path + "wave*.png", encoder="convert", fps=10, output_file="wave.gif")
