#!/usr/bin/env python

""" Program za nukleacijo in rast zrn po FE analizi v programu SALOME, avtor: Andraž Kocjan
     Inštitut za kovinske materiale in tehnologije (IMT), Lepi pot 11, 1000 Ljubljana
     Junij, 2020 """

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random, math, time
plt.ion()

# Rezultat FE analize v Abaqus-u
PATH = 'C:/sm-2018-w64-0-3/WORK/kode za Salome/Cellular Automata/'
file = 'salome_0.npy'         # oznaka simulacije: rp0_ss800_f5R8c4f8b23N30_7.5W_t1_fin
zero_time_matrix  =  np.load(PATH+file)

#matrix_2D = Matrix_3D[:,8,:,:]

time_frames = 3
Z =                  zero_time_matrix.shape[0]
X =                  zero_time_matrix.shape[1]
Y =                  zero_time_matrix.shape[2]


#imported_time_matrix = Matrix_3D.reshape(time_frames, N, M)

nuk_num=1
old_grid = np.zeros((nuk_num, Z, X, Y))
new_grid = np.zeros((nuk_num, Z, X, Y))
total_matrix=np.zeros((Z,X,Y))

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

st_prah = 0           # 0 - prah
st_talina = 10       # 10 - talina
st_haz = 5           # 5 - HAZ (heat affected zone)
st_sled = 3          # 3 - sled  (sintered track)
""" ---------------------------------------------------- """

Tm_Celsius = 900                        # temperatura tališča v st. Celzija
Tmelt=  Tm_Celsius + 273              # temperatura tališča [K]
Thaz = 1.2 * Tmelt                          # mejna temperatura HAZ-a [K]


def nakljucje(msm):
   for z in range(Z):
      for x in range(X):
         for y in range(Y):
            if msm[z][x][y]==0:
               msm[z][x][y]=random.randint(0,10)     
   return msm

""" ~~~~~~~~~~~~~~~~~~~~~~~~ za prikaz staljene cone ~~~~~~~~~~~~~"""
P=np.zeros((Z,X,Y))
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


taula=np.zeros((Z,X,Y))
def taljenje(u, Tm):
   global taula
   plt.clf()
   for z,k in enumerate(u):
      for x,i in enumerate(k):
         for y,j in enumerate(i):
            if j>=Tm:
               taula[z][x][y]= 0  # talina
            else:
               taula[z][x][y]= -1  # prah ALI h.a.z. ALI sled
   return taula


grain_counter=0
def show_nuclei(interface, bulk, timesnap):
      global old_grid, new_grid, old_fresh, novi_nukleusi, grain_counter
      novi_nukleusi=[]
      live= nakljucje(taula)
      time_step = 1
      #TEMP_now = Matrix_3D[timesnap]
      TEMP_now = np.load('salome_'+str(timesnap)+'.npy')
      """~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ dodano v 2020 ~~"""
      try:
         #TEMP_next = Matrix_3D[timesnap + time_step]
         TEMP_next = np.load('salome_'+str(timesnap+time_step)+'.npy')
         print('TEMP_next.shape: ', TEMP_next.shape)
      except IndexError:
         raise IndexError('zmanjkalo je časovnih korakov...')
      """~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
      for k in range(Z):
         for i in range(X):
            for j in range(Y):
               if taula[k][i][j]==0 and total_matrix[k][i][j]==0 and (live[k][i][j]/1 < interface[k][i][j] or live[k][i][j]/1 < bulk[k][i][j])and TEMP_next[k][i][j] < TEMP_now[k][i][j]:              
                  novi_nukleusi.append('nuk')
                  grain_counter +=1
                  #print(grain_counter)
                  old_fresh= np.zeros((1, Z,X,Y)).astype(int)
                  old_fresh[0][k][i][j]=1
                  old_grid=np.concatenate([old_grid, old_fresh])
      
      #print('novi nukleusi:   ',len(novi_nukleusi),'     grain_counter:  ',grain_counter,'    old_grid.shape:   ',old_grid.shape)
                  
      new_fresh=np.zeros((len(novi_nukleusi),Z,X,Y)).astype(int)
      new_grid=np.concatenate([new_grid,new_fresh])
      return old_grid


""" Definicija smeri kristalizacije
            
                   dlup    up    drup
                          \   |    /
                           \  |   /
                left -----  0,0  ----- right
                           /   |   \
                          /    |     \
                dldown  down  drdown

Legenda:

left ..... horizontal-left ,  right .... horizontal-right
up ..... vertical-up  ,  down .... vertical-down

primer za diagonale:   dlup = d-l-up ===> d.... diagonal, l .... left, up ... up
"""
vertical=1
horizontal=1
diagonal_right=1
diagonal_left=1

z_00=1

orient_matrix = [vertical, horizontal, diagonal_right, diagonal_left, z_00]

def Growth(count, om):
   global old_grid, new_grid
   
   for k in range(Z):
      for i in range(X):
           for j in range(Y):
               if old_grid[count][k][i][j]==1 and taula[k][i][j]==0 and total_matrix[k][i][j]<2:     # pogoj za zacetek rasti - obstoj nukleusa (vrednost celice = 1)
                 new_grid[count][k][i][j]=1
                 # smer: --------------------------------------- up ------
                 try:
                     if new_grid[count][k][i][j+om[0]]==0:
                         new_grid[count][k][i][j+om[0]]= 1
                 except IndexError:
                     pass
                 # smer: ------------------------------------ down -----
                 try:
                     if new_grid[count][k][i][abs(j-om[0])]==0:
                         new_grid[count][k][i][abs(j-om[0])]= 1
                 except IndexError :
                     pass
                 # smer: -------------------------------------- left ------
                 try:
                     if new_grid[count][k][i-om[1]][j]==0:
                         new_grid[count][k][i-om[1]][j]=1
                 except IndexError :
                     pass
                 # smer: ------------------------------------- right ------
                 try:
                     if new_grid[count][k][i+om[1]][j]==0:
                         new_grid[count][k][i+om[1]][j]=1
                 except IndexError :
                     pass
                  # smer: -------------------------------------  Z-axis UP (3D)------
                 try:
                     if new_grid[count][k+om[4]][i][j]==0:
                         new_grid[count][k+om[4]][i][j]=1
                 except IndexError :
                     pass
                  # smer: -------------------------------------  Z-axis DOWN (3D)------
                 try:
                     if new_grid[count][k-om[4]][i][j]==0:
                         new_grid[count][k-om[4]][i][j]=1
                 except IndexError :
                     pass
                 
                 """
                 # diagonals
                 # smer: ------------------------------------ dlup ------
                 try:
                     if new_grid[count][abs(i-om[3])][j-om[3]]==0:
                         new_grid[count][abs(i-om[3])][j-om[3]]=1
                 except IndexError :
                     pass
                     
                 # smer: ------------------------------------ drup ------
                 try:
                     if new_grid[count][abs(i-om[2])][j+om[2]]==0:
                         new_grid[count][abs(i-om[2])][j+om[2]]=1
                 except IndexError :
                     pass
                 # smer: --------------------------------- dldown ------
                 try:
                     if new_grid[count][i+om[2]][abs(j-om[2])]==0:
                         new_grid[count][i+om[2]][abs(j-om[2])]=1
                 except IndexError :
                     pass
                 # smer: --------------------------------- drdown ------
                 try:
                     if new_grid[count][i+om[3]][j+om[3]]==0:
                         new_grid[count][i+om[3]][j+om[3]]=1
                 except IndexError :
                     pass
                 """
                 
   old_grid[count]=new_grid[count].copy()
   return old_grid[count]


"""
def Growth_2D(count):
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
 """

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

def vsota(nuk_rate, om):
   global total_matrix
   plt.clf()
   """ Commenting cancels the grain-growth and only formation of nuclei is observed"""
   for i in range(nuk_rate):
      Growth(i, om)
   
   total_matrix=sum(old_grid[i] for i in range(nuk_rate))
   plt.imshow(total_matrix)


tt=-1
def animate(i):
   global tt
   tt +=1
   plt.clf()
   try:
      #T=Matrix_3D[tt]
      T = np.load('salome_'+str(tt)+'.npy')
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
        #T=Matrix_3D[tt]
        T=np.load('salome_'+str(tt)+'.npy')
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
#ani= animation.FuncAnimation(fig, animate, interval=100)  # simulate new_grid data coming in..
#plt.show()



""" Analogija brez animacije - samo računanje ---------------------------------- """

def vsota_brez_animacije(nuk_rate, om):
   global total_matrix
   for i in range(nuk_rate):
      Growth(i, om)
   total_matrix=sum(old_grid[i] for i in range(nuk_rate))
   return total_matrix

tt=-1 
def brez_animacije(cas, om):
   global tt
   start_time = time.time()
   while tt < cas-1:
      tt+=1
      #T=Matrix_3D[tt]
      T = np.load('salome_'+str(tt)+'.npy')
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
      total_matrix=vsota_brez_animacije(nuks, om)
      #for i in range(N):
         #for j in range(M):
          #  if old_grid[i][j]==1:
               #print(dolzine[i][j], tt)
   np.save('total_matrix_3D_SALOME.npy', total_matrix)
   elapsed_time = time.time() - start_time
   print('Python time:   ',round(elapsed_time, 3),'  sec.')
      
brez_animacije(time_frames, orient_matrix)

def c():
   rezultat=np.load('total_matrix_2D_ABAQUS.npy')
   plt.imshow(rezultat)
   plt.show()

#c()

""" --------------------------------------------------------------------------------------------- """
"""
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
"""
# *********************************************************************************

""" Creating a .mp4 video of simulation results -----------------------------"""

import matplotlib
#matplotlib.use("Agg")
#from matplotlib.animation import FFMpegWriter

def animacija_sledi_film(roll):  
    #print('---------------- frame number:  ',tt,' -----------------') 
    
     #plt.clf()
     #T=Matrix_3D[roll]
     T = np.load('salome_'+str(roll)+'.npy')
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


def animacija_rasti_zrn(roll, om):
   plt.clf()
   try:
      #T=Matrix_3D[roll]
      T = np.load('salome_'+str(roll)+'.npy')
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
   vsota(nuks, om)
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
#ani= animation.FuncAnimation(fig, animacija_rasti_zrn(None, orient_matrix), interval=100)  # simulate new_grid data coming in..

#animacija_rasti_zrn(400)
