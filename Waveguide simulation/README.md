# **波導管電場磁場模擬與視覺化**

## **I. 參數設定**
```
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from mpl_toolkits.mplot3d import Axes3D

#----------------------------------------
#Start of parameter setting(SI-unit)

#Shape of Waveguide(Rectangle) L:z-dir WD:y-dir H:z-dir
L = 0.5; WD = 0.05; H = 0.05

#The Number of points for oberservation on each dir, delta: thickness of the waveguide
numberx = 10; numbery = 10; numberz = 30; delta = 0

#Wave Property in vaccum
f = 5*10**9; w = 2*math.pi*f; C = 3*10**8; mu = 4*math.pi*(10**-7)

#WaveGuide's Mode(mn) & Field to Simulate
mode = 'TE'; m = 1; n = 0
Bfield = True
Efield = True

#Parameters for Observation
size = (9, 12)
fig = pylab.figure(figsize=size)
ax = Axes3D(fig)

'''
#ax.view_int(deg1,deg2) - deg1:inclination; deg2:rotaion(on x-y plane) - unit:deg
#if there is nothing in the parentheses of ax.view_init(), system default(30,-60).

'''
ax.view_init(0, 90)

#Range of display on z-dir (Where to Zoom In)
Range = (0.4, 0.5)
```
