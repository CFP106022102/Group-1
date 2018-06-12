# **波導管電場磁場模擬與視覺化程式說明**
## **0. 檔案使用**
- 如有多個檔案，請見檔案備註，*upload序數最末者*，即為應使用之檔案。

## **I. 環境設定**
```
import numpy as np
import matplotlib.pylab as pylab
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import yt
```
## **II. 參數設定(SI-unit)**
1.波導管性質
- L：z方向長度、WD：y方向長度、H：x方向長度、delta：波導管管壁厚度。註：z方向為波傳遞方向。
```
L = 0.5; WD = 0.05; H = 0.05; delta = 0
```

2.入射波物理性質
- f：入射波頻率；w：入射波角頻率；C:真空中光速；mu：真空磁導率。
```
f = 5*10**9; w = 2*math.pi*f; C = 3*10**8; mu = 4*math.pi*(10**-7)。
```

3.選擇個波導管種類
- mode：TE/TM；mn：mode調整(m、n為非負數之整數)。
```
mode = raw_input("Choose the mode of Waveguide (TE or TM)")
m = raw_input("Type the value of m (0~9)")
n = raw_input("Type the value of n (0~9)")
m = int(m)
n = int(n)
```
- number(x、y、z)： 三方向上欲觀測的數目點，總數目點為三者之積。
```
numberx = 100; numbery = 100; numberz = 100; 
```
- size：呈現之圖片大小(水平，垂直)；ax.view_int：欲觀察3D圖像角度(傾角, 旋轉)，單位：deg。
```
size = (9, 12)
fig = pylab.figure(figsize=size)
ax = Axes3D(fig)
ax.view_init(0, 90)
```
- Range：欲顯示在圖片上之Z軸部分(可作為放大效果使用)。
```
Range = (0.4, 0.5)
```
## **III. 開始模擬**
- 對步驟II填入參數進行檢查。

    (1)入射波頻率須超過底限頻率。

    (2)不得出現TM10或TM01模式。
```
#----------------------------------------
#Check Before Simulation

#TM-mode Check
stop = False

if mode == 'TM' and (m == 0 or n == 0):
    print('Error:No TM10(01) mode!')
    stop = True
    
#Cutoff frequency Check
h2 = (m*math.pi/H)**2+(n*math.pi/WD)**2
cutoff = (C*math.sqrt((m*math.pi/H)**2+(n*math.pi/WD)**2))/(2*math.pi)

if ((w/C)**2-h2) <= 0:
    print('Cutoff frequency:',cutoff,'Hz')
    print('Error: Under Cutoff frequency!')
    print('Instruction: Please use HIGHER freqncy.')
elif stop == True:
    print('Stop simulating...')
else:
    print('Cutoff frequency:',cutoff,'Hz')
    print('Start simulating...')

Gamma = 1j * math.sqrt((w/C)**2-h2) #Gamma: i*kz

#End of Check Before Simulation
#----------------------------------------
```

- 開始計算電場與磁場。
```
#Constructing Space
x=np.linspace(0+delta, H-delta , numberx)
y=np.linspace(0+delta, WD-delta, numbery)
z=np.linspace(0+delta, L-delta , numberz)

X, Y, Z = np.meshgrid(x, y, z)

#TE-mode Simulation
if mode == 'TE':
    
    B0 = 1
    alpha = math.pi/4
    BT = B0*(math.cos(alpha)+1j*math.sin(alpha))
    ex = np.exp(-Gamma*Z)
    
    Bz0 = ( BT * np.cos(m/H*math.pi*X) * np.cos(n/WD*math.pi*Y) * ex).real
    Bx0 = ( BT * (Gamma/h2)*(m*math.pi/H ) * np.sin(m*math.pi/H*X) * np.cos(n/WD*math.pi*Y) * ex).real
    By0 = ( BT * (Gamma/h2)*(n*math.pi/WD) * np.cos(m*math.pi/H*X) * np.sin(n/WD*math.pi*Y) * ex).real
    
    Ez0 = 0 * X 
    Ex0 = ( BT * (( 1j)*w*mu/h2) * (n*math.pi/WD) * np.cos(m*math.pi/H*X) * np.sin(n/WD*math.pi*Y) * ex).real
    Ey0 = ( BT * ((-1j)*w*mu/h2) * (m*math.pi/H ) * np.sin(m*math.pi/H*X) * np.cos(n/WD*math.pi*Y) * ex).real
    
#TM-mode Simulation
elif mode == 'TM':
    
    E0 = 1
    alpha = math.pi/4
    ET = E0*(math.cos(alpha)+1j*math.sin(alpha))
    ex = np.exp(-Gamma*Z)

    Ez0 = ( ET * np.sin(m/H*math.pi*X) * np.sin(n/WD*math.pi*Y) * ex).real
    Ex0 = ( ET * (-Gamma/h2) * (m*math.pi/H) * np.cos(m*math.pi/H*X) * np.sin(n/WD*math.pi*Y) * ex).real
    Ey0 = ( ET * (-Gamma/h2) * (m*math.pi/WD) * np.sin(m*math.pi/H*X) * np.cos(n/WD*math.pi*Y) * ex).real

    Bz0 = 0 * X 
    Bx0 = ( ET * (( 1j)*w*mu/h2) * (n*math.pi/WD) * np.sin(m*math.pi/H*X) * np.cos(n/WD*math.pi*Y) * ex).real
    By0 = ( ET * ((-1j)*w*mu/h2) * (m*math.pi/H ) * np.cos(m*math.pi/H*X) * np.sin(n/WD*math.pi*Y) * ex).real
```

- 將向量場繪至三維立體圖上。
```
#Field simulation
Bq = ax.quiver(X, Y, Z, Bx0, By0, Bz0, length=0.005, normalize=True, alpha=0.6, color='red',label='B-field')
print('Quivering B-field...')
Eq = ax.quiver(X, Y, Z, Ex0, Ey0, Ez0, length=0.005, normalize=True, alpha=0.6, color='blue',label='E-field')
print('Quivering E-field...')
```

- 繪圖參數。
```
#Visulization
m = str(m); n = str(n)
plt.title('Vector-Field Visualization - '+mode+m+n+' mode')
plt.axis('equal')

ax.set_xlim3d(0, H )
ax.set_ylim3d(0, WD)
ax.set_zlim3d(Range)
ax.set_xlabel('X axis(m)')
ax.set_ylabel('Y axis(m)')
ax.set_zlabel('Z axis(m)')
ax.legend(loc='best')
```
- yt繪製並存取圖形。
```
#B-field cross-section
mag = Bz0**2 + Bx0**2 + By0**2
data = dict( Bz = (Bz0),By = (By0),Bx = (Bx0), Magnitude = (mag))
print data.keys()
bbox = np.array([[0, 0.05], [0, 0.5],[0, 0.5]])
ds = yt.load_uniform_grid(data, data['Magnitude'].shape , length_unit="Mpc", bbox = bbox)

slc = yt.SlicePlot(ds, "y", ['Magnitude'],center = (0.025,0.25,0.25))
slc.set_log('Magnitude', False)
slc.save('mode '+mode+m+n+' B-field '+'x-z cross-section')
slc = yt.SlicePlot(ds, "x", ['Magnitude'],center = (0.025,0.25,0.25))
slc.set_log('Magnitude', False)
slc.save('mode '+mode+m+n+' B-field '+'y-z cross-section')

#E-field cross-section
mag2 = Ez0**2 + Ex0**2 + Ey0**2
data = dict( Ez = (Ez0),Ey = (Ey0),Ex = (Ex0), Magnitude = (mag2))
print data.keys()
bbox = np.array([[0, 0.05], [0, 0.5],[0, 0.5]])
ds = yt.load_uniform_grid(data, data['Magnitude'].shape , length_unit="Mpc", bbox = bbox)

slc = yt.SlicePlot(ds, "y", ['Magnitude'],center = (0.025,0.25,0.25))
slc.set_log('Magnitude', False)
slc.save('mode '+mode+m+n+' E-field '+'x-z cross-section')
slc = yt.SlicePlot(ds, "x", ['Magnitude'],center = (0.025,0.25,0.25))
slc.set_log('Magnitude', False)
slc.save('mode '+mode+m+n+' E-field '+'y-z cross-section')
```
- 程式模擬完成。
```
print('End of Simulation')
```
