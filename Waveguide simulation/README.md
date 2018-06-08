# **波導管電場磁場模擬與視覺化程式說明**

## **I. 環境設定**
```
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from mpl_toolkits.mplot3d import Axes3D
```
## **II. 參數設定(SI-unit)**
1.波導管性質
- L：z方向長度、WD：y方向長度、H：x方向長度、delta：波導管管壁厚度。註：z方向為波傳遞方向
```
L = 0.5; WD = 0.05; H = 0.05; delta = 0
```

2.入射波物理性質
- f：入射波頻率；w：入射波角頻率；C:真空中光速；mu：真空介電係數。
```
f = 5*10**9; w = 2*math.pi*f; C = 3*10**8; mu = 4*math.pi*(10**-7)。
```

3.模擬參數調整
- 
```
mode = 'TE'; m = 1; n = 0
Bfield = True
Efield = True
```
- number(x、y、z)： 三方向上欲觀測的數目點，總數目點為三者之積。
```
numberx = 10; numbery = 10; numberz = 30; 
```

```
#Parameters for Observation
size = (9, 12)
fig = pylab.figure(figsize=size)
ax = Axes3D(fig)
ax.view_init(0, 90)
```
```
#Range of display on z-dir (Where to Zoom In)
Range = (0.4, 0.5)
```
