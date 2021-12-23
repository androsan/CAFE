#!/usr/bin/env python

""" Program za nukleacijo in rast zrn po FE analizi v programu ABAQUS, avtor: Andraž Kocjan
     Inštitut za kovinske materiale in tehnologije (IMT), Lepi pot 11, 1000 Ljubljana
     Februar, 2019 """

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random, math, time
plt.ion()

# Rezultat FE analize v Abaqus-u
PATH = 'C:/sm-2018-w64-0-3/WORK/SLM_2D_Source_Files/post processing database/0002/'
file = 'SLM_2D_Source_0002.npy'
Matrix_3D  =  np.load(PATH+file)

matrix_2D = Matrix_3D[:,8,12:28,:]
 
#  N - širina, M - dolžina domene v številu celic, ki se ujema z FE mrežo

N = 16 #matrix_2D.shape[1]
M = matrix_2D.shape[2]

time_frames = matrix_2D.shape[0]
imported_time_matrix = matrix_2D.reshape(time_frames, N, M)

nuk_num=1
old_grid = np.zeros((nuk_num, N, M))
new_grid = np.zeros((nuk_num, N, M))
total_matrix=np.zeros((N,M))

""" -----------------------------------------------------------------------------------------------------------------------------"""
ns_max = 5e10      # število novih nukleusov pri optimalni dT [m^-2]
dTs_max0 = 2       # optimalna dT - kjer nastane največ nukleusov [K]
dTs_sigma = 0.5    # standardna deviacija Gaussove porazdelitve števila nukleusov v odvisnosti od dT [K]

nv_max = 5e14      # parametri nukleacije v tekočem analogija s parametri za površino
dTv_max0 = 2
dTv_sigma = 0.5 

def nukleacija_povrsina(pp):   # heterogena nukleacija
   Ns_Nas =  ns_max/(math.sqrt(2*3.1428))*math.e**(-(pp-dTs_max0)**2/(2*dTs_sigma**2))
   return Ns_Nas

def nukleacija_volumen(vv):    # homogena nukleacija
   Nv_Nav =  nv_max/(math.sqrt(2*3.1428))*math.e**(-(vv-dTv_max0)**2/(2*dTv_sigma**2))
   #Nv_Nav=0
   return Nv_Nav

nukleusi_pov = np.vectorize(nukleacija_povrsina)
nukleusi_vol = np.vectorize(nukleacija_volumen)
""" ---------------------------------------------------------------------------------------------------------------------------------"""
""" OZNAKE STANJA CELIC """

st_prah = 0          # 0 - prah
st_talina = 10       # 10 - talina
st_haz = 5           # 5 - HAZ (heat affected zone)
st_sled = 3          # 3 - sled  (sintered track)
""" ---------------------------------------------------- """

Tm_Celsius = 1507                            # temperatura tališča v st. Celzija
Tmelt=  Tm_Celsius + 273                    # temperatura tališča [K]
Thaz = 1.2 * Tmelt                          # mejna temperatura HAZ-a [K]


def nakljucje(msm):
   for x in range(N):
      for y in range(M):
         if msm[x][y]==0:
            msm[x][y]=random.randint(0,10)     
   return msm


P=np.zeros((N,M))
def phase(u, Tm):
   global P
   plt.clf()
   for x,i in enumerate(u):
      for y,j in enumerate(i):
         if j>=Tm:
            P[x][y]= st_talina  # talina (melt_pool)
         elif j<Tm and j>=Thaz:
            P[x][y]= st_haz   # HAZ
         else:
            P[x][y]= st_prah  # prah ALI h.a.z. ALI sled
   return P


taula=np.zeros((N,M))
def taljenje(u, Tm):
   global taula
   plt.clf()
   for x,i in enumerate(u):
      for y,j in enumerate(i):
         if j>=Tm:
            taula[x][y]= 0  # talina
         else:
            taula[x][y]= -1  # prah ALI h.a.z. ALI sled
   return taula


grain_counter=0
def show_nuclei(interface, bulk, timesnap):
      global old_grid, new_grid, old_fresh, novi_nukleusi, grain_counter
      novi_nukleusi=[]
      live= nakljucje(taula)
      time_step = 1
      TEMP_now = imported_time_matrix[timesnap]
      """~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ dodano v 2020 ~~"""
      try:
         TEMP_next = imported_time_matrix[timesnap + time_step]
      except IndexError:
         pass
      """~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
      for i in range(N):
            for j in range(M):
               if taula[i][j]==0 and total_matrix[i][j]==0 and (live[i][j]/1 < interface[i][j] or live[i][j]/1 < bulk[i][j])and TEMP_next[i][j] < TEMP_now[i][j]:              
                  novi_nukleusi.append('nuk')
                  grain_counter +=1
                  #print(grain_counter)
                  old_fresh= np.zeros((1, N,M)).astype(int)
                  old_fresh[0][i][j]=1
                  old_grid=np.concatenate([old_grid, old_fresh])
      
      #print('novi nukleusi:   ',len(novi_nukleusi),'     grain_counter:  ',grain_counter,'    old_grid.shape:   ',old_grid.shape)
                  
      new_fresh=np.zeros((len(novi_nukleusi),N,M)).astype(int)
      new_grid=np.concatenate([new_grid,new_fresh])
      return old_grid
      


def Growth(count):
   global old_grid, new_grid
   
   for i in range(N):
        for j in range(M):
                if old_grid[count][i][j]==1 and taula[i][j]==0 and total_matrix[i][j]<2:
                    new_grid[count][i][j]=1
                    try:
                        if new_grid[count][i][j+1]==0:
                            new_grid[count][i][j+1]= 1
                    except IndexError:
                        pass
                    try:
                        if new_grid[count][i][abs(j-1)]==0:
                            new_grid[count][i][abs(j-1)]= 1
                    except IndexError :
                        pass
                    try:
                        if new_grid[count][i+1][j]==0:
                            new_grid[count][i+1][j]=1
                    except IndexError :
                        pass
                    try:
                        if new_grid[count][abs(i-1)][j]==0:
                            new_grid[count][abs(i-1)][j]=1
                    except IndexError :
                        pass
                
   old_grid[count]=new_grid[count].copy()
   return old_grid[count]
 

""" ------------------------------------------------------------------------------------------------------------"""				
def show_vg(dT_field):               # it shows the interface growth rate (vg)vs. undercooling temperatures (dT_field)at given time step
   vg=2.03e-4 * dT_field**2 - 0.544e-4 * dT_field
   plt.clf()
   plt.imshow(vg, extent=(np.amin(xim), np.amax(xim), np.amin(yim), np.amax(yim)),  cmap=plt.get_cmap('rainbow'), interpolation='bessel')
   #plt.plot(dT_field[1], vg[1])
   #print('dT= ',dT_field[1][6],'vg= ',vg[1][6])
   return vg

def get_growth_length(dT_field, cas):
   vg=2.03e-4 * dT_field**2 - 0.544e-4 * dT_field
   length= vg*cas
   return length
""" -----------------------------------------------------------------------------------------------------------"""

def vsota(nuk_rate):
   global total_matrix
   plt.clf()
   """ Commenting cancels the grain-growth and only formation of nuclei is observed"""
   for i in range(nuk_rate):
      Growth(i)
   
   total_matrix=sum(old_grid[i] for i in range(nuk_rate))
   plt.imshow(total_matrix)


tt=-1
def animate(i):
   global tt
   tt +=1
   plt.clf()
   try:
      T=imported_time_matrix[tt]
      taljenje(T, Tmelt)
      plt.imshow(taula)
      """
      dTt = T - Tmelt       # termično podhlajanje, kjer je T iz Abaqus-a uvožena matrika dejanskih temperatur posamezne celice
      nuk_matrix_pov=nukleusi_pov(dTt)                     
      nuk_matrix_vol=nukleusi_vol(dTt)
      show_nuclei(nuk_matrix_pov, nuk_matrix_vol, tt)
      """
      #show_vg(dTt)
      #dolzine=get_growth_length(dTt,1)
      #print(dolzine[0])
      """
      nuks=len(old_grid)
      vsota(nuks)
      """
      #for i in range(N):
         #for j in range(M):
          #  if old_grid[i][j]==1:
               #print(dolzine[i][j], tt)
   except IndexError:
      print('Done..')
      ani.event_source.stop()
      
# Odkomentiraj spodnje 4 vrstice, če želiš pognati animacija_sledi
"""
time_frames = time_frames
temp=np.zeros((2,N,M)) # talina
haz=np.zeros((2,N,M))  #  talina +h.a.z. (heat affected zone)
out=np.zeros((2,N,M))   #  talina + sled
"""

tt=0
def animacija_sledi(i):  
    global tt
    tt+=1
    #print('---------------- frame number:  ',tt,' -----------------') 
    if tt<time_frames:
        plt.clf()
        T=imported_time_matrix[tt]
        P=phase(T,Tmelt)
        
        if tt==1:
            temp[0]= P
            #print('TEMP0= ',temp[0],'TEMP1= ',temp[1])
            for j in range(len(temp[0])):
                for i in range(len(temp[0][j])):
                    if temp[0][j][i]==st_talina and temp[1][j][i]==st_prah:
                        out[0][j][i]=st_talina
                    else:
                        out[0][j][i]=temp[1][j][i]
            plt.imshow(temp[0])
        
        else:
            temp[1]= P
            #print('TEMP0= ',temp[0],'TEMP1= ',temp[1])
            for j in range(len(temp[0])):
                for i in range(len(temp[0][j])):
                    if (temp[0][j][i]==st_haz and temp[1][j][i]==st_prah):
                        out[1][j][i]=st_sled
                    else:
                        out[1][j][i]=temp[1][j][i]
                    if out[0][j][i]==st_sled and out [1][j][i]==st_prah:
                        out[1][j][i]=st_sled

            #print('OUT0= ',out[0], 'OUT1= ',out[1]); print('-------------------------------------------------'); print()

            temp[0]=temp[1].copy()
            out[0]=out[1].copy()

            plt.imshow(out[1])
            cmap = plt.cm.jet
            plt.imshow(out[1], interpolation='bessel', cmap=cmap) #  cmap='jet'
         
    else:
        ani.event_source.stop()
        print('konec..')


#fig, ax = plt.subplots()
#ani= animation.FuncAnimation(fig, animate, interval=500)  # simulate new_grid data coming in..
#plt.show()



""" Analogija brez animacije - samo računanje ---------------------------------- """

def vsota_brez_animacije(nuk_rate):
   global total_matrix
   for i in range(nuk_rate):
      Growth(i)
   total_matrix=sum(old_grid[i] for i in range(nuk_rate))
   return total_matrix

tt=-1 
def brez_animacije(cas):
   global tt
   start_time = time.time()
   while tt < cas-1:
      tt+=1
      T=imported_time_matrix[tt]
      taljenje(T, Tmelt)
      #phase(T, Tmelt)
      dTt = T - Tmelt       # termično podhlajanje, kjer je T uvožena matrika dejanskih temperatur posamezne celice
      nuk_matrix_pov=nukleusi_pov(dTt)                     
      nuk_matrix_vol=nukleusi_vol(dTt)
      show_nuclei(nuk_matrix_pov, nuk_matrix_vol, tt)
      #show_vg(dTt)
      #dolzine=get_growth_length(dTt,1)
      #print(dolzine[0])
      nuks=len(old_grid)
      total_matrix=vsota_brez_animacije(nuks)
      #for i in range(N):
         #for j in range(M):
          #  if old_grid[i][j]==1:
               #print(dolzine[i][j], tt)
   np.save('total_matrix_2D_ABAQUS.npy', total_matrix)
   elapsed_time = time.time() - start_time
   print('Python time:   ',round(elapsed_time, 3),'  sec.')
      
brez_animacije(time_frames)

def c():
   rezultat=np.load('total_matrix_2D_ABAQUS.npy')
   plt.imshow(rezultat)
   plt.show()

#c()

""" --------------------------------------------------------------------------------------------- """

# DIAGRAM PORAZDELITVE VELIKOSTI ZRN..
velikosti_zrn=[]
def count(x):      #  sešteje vse 1 v matriki posameznega zrna (npr. x=old[4])in določi njegovo velikost..
   q=[]
   for i in x:
      d=list(i)
      m=d.count(1)
      q.append(m)
   return(sum(q))
for i in old_grid:
	velikosti_zrn.append(count(i))
del velikosti_zrn[0]
print(' Nastalo je',len(velikosti_zrn),' zrn.')
y= np.array(velikosti_zrn)
def make_histogram_x_axis(seq, step):
    interval=[0, step]
    z=step
    for i in seq:
        if i < z:
            pass
        else:
            z += step
            interval.append(z)
    return interval

x=make_histogram_x_axis(y,5)
#(n, bins, patches)= plt.hist(y, bins=x)
#y=[100*i/sum(n)for i in n]
plt.hist(y, bins=x,edgecolor='white')
plt.title("Porazdelitev velikosti zrn")
plt.xlabel("Velikost (a.u.)")
plt.ylabel("Število zrn posamezne frakcije")
plt.show()

# *********************************************************************************

""" Creating a .mp4 video of simulation results -----------------------------"""

import matplotlib
#matplotlib.use("Agg")
from matplotlib.animation import FFMpegWriter

def animacija_sledi_film(roll):  
    #print('---------------- frame number:  ',tt,' -----------------') 
    
     #plt.clf()
     T=imported_time_matrix[roll]
     P=phase(T,Tmelt)
     
     if roll==1:
         temp[0]= P
         #print('TEMP0= ',temp[0],'TEMP1= ',temp[1])
         for j in range(len(temp[0])):
             for i in range(len(temp[0][j])):
                 if temp[0][j][i]==st_talina and temp[1][j][i]==st_prah:
                     out[0][j][i]=st_talina
                 else:
                     out[0][j][i]=temp[1][j][i]
         #plt.imshow(temp[0])
         return temp[0]
     
     else:
         temp[1]= P
         #print('TEMP0= ',temp[0],'TEMP1= ',temp[1])
         for j in range(len(temp[0])):
             for i in range(len(temp[0][j])):
                 if (temp[0][j][i]==st_haz and temp[1][j][i]==st_prah):
                     out[1][j][i]=st_sled
                 else:
                     out[1][j][i]=temp[1][j][i]
                 if out[0][j][i]==st_sled and out [1][j][i]==st_prah:
                     out[1][j][i]=st_sled

         #print('OUT0= ',out[0], 'OUT1= ',out[1]); print('-------------------------------------------------'); print()

         temp[0]=temp[1].copy()
         out[0]=out[1].copy()

         #plt.imshow(out[1])
         #cmap = plt.cm.jet
         #plt.imshow(out[1], interpolation='bessel', cmap=cmap) #  cmap='jet'
         return out[1]

"""
from tkinter import *
root=Tk()
gc=StringVar()
lab1=Label(root, textvariable=gc, font=('Arial', '50')).pack()
gc.set(0)
"""


def animacija_rasti_zrn(roll):
   plt.clf()
   try:
      T=imported_time_matrix[roll]
      taljenje(T, Tmelt)
      dTt = T - Tmelt
      nuk_matrix_pov=nukleusi_pov(dTt)                     
      nuk_matrix_vol=nukleusi_vol(dTt)
      show_nuclei(nuk_matrix_pov, nuk_matrix_vol, roll)
   except IndexError:
      pass
   nuks=len(old_grid)
   #print(nuks)
   #gc.set(nuks)
   vsota(nuks)
   plt.title('št. novo nastalih nukleusov v izbranem časovnem intervalu:         '+str(len(novi_nukleusi))+'\n     števec zrn:         '+str(grain_counter)
             +'\n    old_grid.shape:   '+str(old_grid.shape), loc='right', fontsize=20)
   
   #return total_matrix
   

#fig = plt.figure(figsize = (20,5))
fig, ax = plt.subplots()

def make_results_movie(hitrost, okvir_start, okvir_stop, resolucija):    # hitrost [frames per second - fps],  okvirjev [stevilo okvirjev zajetih v film], resolucija [dpi]
    metadata = dict(title='Simulacija rasti zrn pri 3D tisku', artist='made by Andro', comment='Movie support!')
    writer = FFMpegWriter(fps=hitrost, metadata=metadata)
    
    with writer.saving(fig, "simulacija_rasti.mp4", resolucija):
        for i in range(okvir_start, okvir_stop - 1):
            plt.clf()
            #************************************************ TALINA *********************************
            fig.add_subplot()


            plt.imshow(animacija_rasti_zrn(i))#, interpolation='bessel', cmap=plt.cm.jet)
            #animate(i)
      
            plt.title('T(tališča)=  '+str(Tm_Celsius)+' 'u'\u2103')     #   T(h.a.z.)=  830 'u'\u2103')
            writer.grab_frame()

fps = 50  # slik na sekundo
start_frame = 500   # zaporedna številka začetnega okvirja
stop_frame = 600   # zaporedna številka zadnjega okvirja
dpi = 100  # resolucija
"""
go=time.time()
make_results_movie(fps, start_frame, stop_frame, dpi)
over=time.time()
print('The animation in Python code was made in ',round(over-go, 3),' seconds..')
"""
#ani= animation.FuncAnimation(fig, animacija_rasti_zrn, interval=100)  # simulate new_grid data coming in..

#animacija_rasti_zrn(400)
