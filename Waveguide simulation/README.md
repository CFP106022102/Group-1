# **波導管電場磁場模擬與視覺化程式說明**
## **0. 檔案使用**
- 如有重複檔案，請見檔案備註，upload序數最末者，即為應使用之檔案。

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
- L：z方向長度、WD：y方向長度、H：x方向長度、delta：波導管管壁厚度。註：z方向為波傳遞方向。
```
L = 0.5; WD = 0.05; H = 0.05; delta = 0
```

2.入射波物理性質
- f：入射波頻率；w：入射波角頻率；C:真空中光速；mu：真空介電係數。
```
f = 5*10**9; w = 2*math.pi*f; C = 3*10**8; mu = 4*math.pi*(10**-7)。
```

3.模擬參數調整
- mode：TE/TM；mn：mode調整(m、n為非負數之整數)；B/F field：選擇欲模擬之場(B：磁場、E：電場)。
```
mode = 'TE'; m = 1; n = 0
Bfield = True
Efield = True
```
- number(x、y、z)： 三方向上欲觀測的數目點，總數目點為三者之積。
```
numberx = 10; numbery = 10; numberz = 30; 
```
- size：呈現之圖片大小(長，寬)；ax.view_int：欲觀察3D圖像角度(傾角, 旋轉)，單位：deg。
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

    (1)頻率須超過底限頻率。

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
    
    Bz = ( BT * np.cos(m/H*math.pi*X) * np.cos(n/WD*math.pi*Y) * ex).real
    Bx = ( BT * (Gamma/h2)*(m*math.pi/H ) * np.sin(m*math.pi/H*X) * np.cos(n/WD*math.pi*Y) * ex).real
    By = ( BT * (Gamma/h2)*(n*math.pi/WD) * np.cos(m*math.pi/H*X) * np.sin(n/WD*math.pi*Y) * ex).real
    
    Ez = 0 * X 
    Ex = ( BT * (( 1j)*w*mu/h2) * (n*math.pi/WD) * np.cos(m*math.pi/H*X) * np.sin(n/WD*math.pi*Y) * ex).real
    Ey = ( BT * ((-1j)*w*mu/h2) * (m*math.pi/H ) * np.sin(m*math.pi/H*X) * np.cos(n/WD*math.pi*Y) * ex).real
    
#TM-mode Simulation
elif mode == 'TM':
    
    E0 = 1
    alpha = math.pi/4
    ET = E0*(math.cos(alpha)+1j*math.sin(alpha))
    ex = np.exp(-Gamma*Z)

    Ez = ( ET * np.sin(m/H*math.pi*X) * np.sin(n/WD*math.pi*Y) * ex).real
    Ex = ( ET * (-Gamma/h2) * (m*math.pi/H) * np.cos(m*math.pi/H*X) * np.sin(n/WD*math.pi*Y) * ex).real
    Ey = ( ET * (-Gamma/h2) * (m*math.pi/WD) * np.sin(m*math.pi/H*X) * np.cos(n/WD*math.pi*Y) * ex).real

    Bz = 0 * X 
    Bx = ( ET * (( 1j)*w*mu/h2) * (n*math.pi/WD) * np.sin(m*math.pi/H*X) * np.cos(n/WD*math.pi*Y) * ex).real
    By = ( ET * ((-1j)*w*mu/h2) * (m*math.pi/H ) * np.cos(m*math.pi/H*X) * np.sin(n/WD*math.pi*Y) * ex).real
```

- 將向量場繪至三維立體圖上。
```
#B-Field simulation
if (Bfield & ~stop) == True:
    Bq = ax.quiver(X, Y, Z, Bx, By, Bz, length=0.005, normalize=True, alpha=0.6, arrow_length_ratio=0.3, 
                   color='red', label='B-field')
    print('Quivering B-field...')    

#E-Field simulation
if (Efield & ~stop) == True:
    Eq = ax.quiver(X, Y, Z, Ex, Ey, Ez, length=0.005, normalize=True, alpha=0.6, arrow_length_ratio=0.3,
                   color='blue', label='E-field') 
    print('Quivering E-field...')
```

- 繪圖參數。
```
#Visulization
m = str(m); n = str(n)
plt.title('Vector-Field Visualization - '+ (mode+m+n) +' mode')
plt.axis('equal')

ax.set_xlim3d(0, H )
ax.set_ylim3d(0, WD)
ax.set_zlim3d(Range)
ax.set_xlabel('X axis(m)')
ax.set_ylabel('Y axis(m)')
ax.set_zlabel('Z axis(m)')

#plot waveguide L:z-dir; WD:y-dir; H:x-dir
#if ~stop == True:
op=1    
ax.plot([0,0,H,H,0],[0,WD,WD,0,0],[0,0,0,0,0],'--k',alpha=op)
ax.plot([0,0,H,H,0],[0,WD,WD,0,0],[L,L,L,L,L],'--k',alpha=op)
ax.plot([0,0,0,0,0],[0,WD,WD,0,0],[0,0,L,L,0],'--k',alpha=op)
ax.plot([H,H,H,H,H],[0,WD,WD,0,0],[0,0,L,L,0],'--k',alpha=op, label='WaveGuide')
#ax.axis('off')

ax.legend(loc='center right')
#ax.legend(loc='best')
plt.show()
```
- 程式模擬完成。
```
print('End of Simulation')
```
