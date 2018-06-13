此資料夾為各階段數值解的檔案，最後607磁場、電場數值解為最後的測試結果，但結果始終與標準解有所差異。

1.5/21磁場單一平面數值解：
    
    在網路找到求解二階偏微分方程式，藉由參考其作法，重新定義基本參數，在了解其運算過程後，依照波導的形式，重新寫入新的邊界條件，
    並加入一階偏微分求得磁場x.y方向、電場x.y方向，算出TE模式下，z方向剖面的數值解。

2.5/31磁場、電場數值解:

    在先前的數值解解出來後，我們將其數值作圖觀察，發現到趨勢大小與預期並不相符，在重新分析波導的解後，發現到我們並沒有考量到解本
    身是為複數的形式，在521磁場單一平面數值解的大小，是為純考量實數部分的運算，因此，改變最初的矩陣定義，從原先的實數0的矩陣，改
    為複數0的矩陣運算求解。並在檢查程式的過程中，發現到邊界條件的輸入錯誤，重新輸入正確的邊界條件。

3.6/07磁場、電場數值解:

    改變x.y方向解方程式的方法，使用不同函數解出方程式。在先前的531磁場、電場數值解中所求得磁場x.y方向、電場x.y方向，依然與標準解
    有差異。以磁場來說，在TE10模式中其y方向的數值解應為0，但在我們的的數值解當中，其數值並不為0，但在檢查531磁場、電場數值解過後
    ，推測是在求Bx.By.Ex.Ey的偏微分方式是有問題的，因此改為使用函數gradient來重新運算，但測試結果依然無法求得接近的結果。
    最終在多次的檢查程式碼後，我們依然認為在最前段求解Bz的二次偏微分方程式求解，應該是沒有問題的，但最後在偏微分各點數值時出了問題
    ，可能與其運算定義上還要再多加釐清。


I.環境設定

    import numpy as np
    import numpy.linalg as lg
    import math
 
 II.參數設定(SI-unit)
    1.波導管性質
   
      a = 5*10**-2
      b = 5*10**-2
      c=5*10**-1
   
   2.入射波物理性質
    
     C=3*10**8
     fq=5*10**9
     w=2*math.pi*fq
     kz=np.sqrt(w*w/C/C-(math.pi*math.pi*(m*m/a/a+n*n/b/b)))
     B0=1
     complexcoe=w*w/C/C-kz*kz
   
   3.選擇波導管種類
    
     mode：TE/TM；mn：mode調整(m、n為非負數之整數)
     m=1
     n=0
     
 III.計算數值解
   1.定義函數  

     def pde_waveguide_TE(f,g,BFx,Domain,Mx,My,Mz,c,tol,MaxIter):
  
   2.設定參數
     
     Bx=[]
     By=[]
     Bz=[]
     Ex=[]
     Ey=[]
     Ez=[]
     z=np.linspace(0,c,Mz)

     bx0, bxf, by0, byf = BFx[0], BFx[1], BFx[2], BFx[3]
     x0, xf, y0, yf = Domain[0], Domain[1], Domain[2], Domain[3]
     
     hx = float(xf-x0)/Mx
     hy = float(yf-y0)/My

     x = [x0+j*hx for j in range(1,Mx)]
     y = [y0+i*hy for i in range(1,My)]

     Mx, My = Mx+1, My+1
     
   3.設定複數
       
     for intz in z:
       
       #set complex number to calculate 
        u = np.zeros([My,Mx])*1j
        F = np.zeros([My,Mx])*1j
        G = np.zeros([My,Mx])*1j
        u0 = np.zeros([My,Mx])*1j
        ux = np.zeros([My,Mx])*1j
        uy = np.zeros([My,Mx])*1j
        ez = np.zeros([My,Mx])*1j
        ex = np.zeros([My,Mx])*1j
        ey = np.zeros([My,Mx])*1j
   
   4.設定邊界條件
   
     j = 1
     for xj in x:
        u[0,j], u[My-1,j] = by0(xj), byf(xj)
        j+=1
     i = 1
     for yi in y:
        u[i,0], u[i,Mx-1] = bx0(yi), bxf(yi)
        i+=1
   
   5.設定z方向係數
       
     #adding the z direction coefficient 
     u = u*(np.cos((kz*intz))+1j*np.sin(kz*intz))
        
        
   6.平衡邊界值
        
      sum_of_bv = sum(u[0,:])+sum(u[My-1,:])+sum(u[1:My-1,0])+sum(u[1:My-1,Mx-1])
      u[1:My-1,1:Mx-1] = float(sum_of_bv)/(2*(Mx+My-2))
        
   7.設定f(xj,yi) & g(xj,yi)
        
     for i in range(1,My-1):
         for j in range(1,Mx-1):
             F[i,j], G[i,j] = f(x[j-1],y[i-1]), g(x[j-1],y[i-1])   
   8.定義變數
   
     dx2, dy2 = hx**2, hy**2
     dxy2 = 2*(dx2+dy2)
     rx, ry = dx2/dxy2, dy2/dxy
     rxy = rx*dy2   
        
   9.解u(Bz)數值解     
   
     for itr in range(MaxIter):
         for i in range(1,My-1):
             for j in range(1,Mx-1):
                 u[i,j] = ry*(u[i,j+1]+u[i,j-1])+rx*(u[i+1,j]+u[i-1,j])+rxy*(G[i,j]*u[i,j]-F[i,j])
             
         Err = abs(u-u0)

            
         if (itr>1) & (Err.max()<tol):
                break
            
         u0=u        
        
   10.用差分法計算Bx,By,Ex,Ey                  
        
      for i in range(Mx-1):
          for j in range(My-1): 
              ux[i,j] = (u[i,j+1]-u[i,j])/hx
              uy[i,j] = (u[i+1,j]-u[i,j])/hy
        
        
      for i in range(Mx-1):
          for j in range(My-1): 
              ex[i,j] = (u[i+1,j]-u[i,j])/hy
              ey[i,j] = (u[i,j+1]-u[i,j])/hx        
    
   11.算出完整的Bx,By,Ex,Ey   
        
      ux = ux*kz*1j/(w*w/C/C-kz*kz)
      uy = uy*kz*1j/(w*w/C/C-kz*kz)
      ex = ex*w*1j/(w*w/C/C-kz*kz)
      ey = ey*(-w)*1j/(w*w/C/C-kz*kz)        
        
   12.將計算結果取實數部分  
        
      u = np.real(u[1:My-1,1:Mx-1])
      ux = np.real(ux[1:My-1,1:Mx-1])
      uy = np.real(uy[1:My-1,1:Mx-1])
      ex = np.real(ex[1:My-1,1:Mx-1])
      ey = np.real(ey[1:My-1,1:Mx-1])
      ez = np.real(ez[1:My-1,1:Mx-1])        
    
   13.將每單一平面數值解整合
   
      Bx.append(ux)
      By.append(uy)
      Bz.append(u)
      Ex.append(ex)
      Ey.append(ey)
      Ez.append(ez)   

   14.定義函數結束
   
      return Bx,By,Bz,Ex,Ey,Ez   
-----------------------------------------------------------------------------------------------------------------------------

IV.輸出結果

   1.定義邊界

     BF=[lambda y:B0*np.cos(m*math.pi*0/a)*np.cos(n*math.pi*y/b)
     ,lambda y:B0*np.cos(m*math.pi*a/a)*np.cos(n*math.pi*y/b)
     ,lambda x:B0*np.cos(m*math.pi*x/a)*np.cos(n*math.pi*0/b)
     ,lambda x:B0*np.cos(m*math.pi*x/a)*np.cos(n*math.pi*b/b)]  
     Domain = [0,a,0,b]
     
   2.計算完整數值解輸出
     
     Bx,By,Bz,Ex,Ey,Ez = pde_waveguide_TE(lambda x,y:0.0, lambda x,y:complexcoe, BF,Domain,numberx+1,numbery+1,numberz,c,10**(-8),1000)
