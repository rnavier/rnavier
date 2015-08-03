import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline as intpl

#Compare a slice of 2d data with 1d

# Set the angle of slice of 2d plot
phideg = 0
t1d = np.loadtxt("sod1d_1dslice_x_0000_0000.out")
ngr = len(t1d)/2;
t1d.shape = (2, ngr,20)
t = np.loadtxt("sod2d_xy_2dslice_xy_0000.out")
ngr = np.sqrt(len(t)/2);
t.shape=(2,ngr,ngr,20)
x = t[1,:,0,1]
y = t[1,0,:,2]
phi = np.pi/180.*phideg
xr = x/np.cos(phi) 
yr = x*np.tan(phi)
e = t[1,:,:,4]
f = intpl(x,y,e)

plt.xlim([-1,1])
plt.ylim([0,30])
plt.xlabel("x")
plt.xlabel("e")
plt.title("Energy")
fr = np.zeros(len(x))
for i in range(len(x)):
    fr[i] = f(x[i],yr[i])

plt.plot(xr,fr)
plt.plot(t1d[1,:,1],t1d[1,:,4])
plt.show()
