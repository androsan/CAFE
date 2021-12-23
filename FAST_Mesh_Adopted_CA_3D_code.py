''' FAST code for Mesh_Adopted 3D CA, August 2020 '''
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import multi_dot
import random, math, time
plt.ion()

# FE analysis result from SALOME
PATH = 'C:/sm-2018-w64-0-3/WORK/kode za Salome/Cellular Automata/mesh_dependency_experiment/'

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

mapa =                       'CA_'                                                    #  Folder name
flashies=                     'flashies'

result =                       mapa + ' result .npy'                                 #   NPY file with results                           

time_range =             (0, 15)                                                 #  Range of time frames;  unit: - [-]
#time_step =                 7*5e-05                         # 0.0015625 for CA_4 unit msec.             #  FE time step; unit: SECOND [s]
time_step =  0.056    # 0.00385

cell =                           5e-05 / 3                                                  #   FE cell size; unit:  METERS [m]

Tmelt_Celsius =       1650             #     2217                                         #   Melting point; unit:  Deg.  CELSIUS  [deg. C]
dTliquidus =                                   65                                          #   Difference between Absolute Liquidus Line and Melting point; unit: Deg. CELSIUS [deg. C]


reduce_domain_size = False

orientacije = np.array(['001', '00_1', '010', '0_10', '100', '_100', #])    #   1st shell (6-nearest neighbours)---> grain-growth directions with respect to domain axes('ZXY')
                     '011', ])#'01_1', '0_11', '0_1_1', '101', '110', '10_1', '1_10', '_101', '_110', '_10_1', '_1_10',  #])    #   2nd shell (12-nearest neighbours)directions..                                            
                     #'111', '11_1', '1_11', '1_1_1', '_111', '_11_1', '_1_11', '_1_1_1',   ])      #   3rd shell (8-nearest neighbours)directions..

np.random.seed(75689803)                                               #   Random Seed

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
Tmelt= Tmelt_Celsius + 273                                                  #   Melting point; unit:  KELVIN [K]
L_1st_shell =   cell                                                              #   Distance to six nearest neighbours;  unit:  METERS [m]
L_2nd_shell =  math.sqrt(2)* cell                                       #   Distance to twelve second-nearest neighbours;  unit:  METERS [m]
L_3rd_shell =  math.sqrt(3)*cell                                         #   Distance to eight third-nearest neighbours;  unit:  METERS [m]

def m(x):
    mic = np.load(PATH+mapa+'/flashies/flashy_snap_'+str(x)+'.npy')
    #plt.title(u'\u03B1 = '+a)
    plt.imshow(mic[0])



''' domain limits ----> reduction of FE domain to size a bit larger than the molten track, where the crystallization occurs .. '''

def Domain_Size_Reduction(domain, threshold):
   global z_min, z_max, x_min, x_max, y_min, y_max
   ax2=np.where(np.any(domain > threshold, axis=2))
   ax1=np.where(np.any(domain > threshold, axis=1))
   z_min = np.min(ax2[0]);  z_max = np.max(ax2[0])+1
   x_min = np.min(ax2[1]);  x_max = np.max(ax2[1])+1
   y_min = np.min(ax1[1]);  y_max = np.max(ax1[1])+1
   
   reduced_domain = domain[z_min:z_max, x_min:x_max, y_min:y_max]
   return reduced_domain



if reduce_domain_size:
   limits = [ [], [], [], [], [], [] ]
   for i in time_range:
      Domain_Size_Reduction(np.load(PATH+mapa+'/'+'salome_'+str(i)+'.npy'), Tmelt)
      limits[0].append(z_min); limits[1].append(z_max)
      limits[2].append(x_min); limits[3].append(x_max)
      limits[4].append(y_min); limits[5].append(y_max)

   z_min, z_max, x_min, x_max, y_min, y_max = min(limits[0]), max(limits[1]), min(limits[2]), max(limits[3]), min(limits[4]), max(limits[5])

else:

   domena = np.full((1,201,201), Tmelt+5); np.save(PATH+mapa+'/'+'salome_'+str(time_range[0])+'.npy', domena)
   #domena = np.load(PATH+mapa+'/'+'salome_'+str(time_range[0])+'.npy')
   z_min =                           0
   z_max =                          domena.shape[0]

   x_min =                           0
   x_max =                          domena.shape[1]

   y_min =                           0
   y_max =                          domena.shape[2]

''' ........................................... domain constraints - FINAL and Saving as .Npy file ..............................................................'''
#z_min =    0
#z_max =   23
try:
   np. save(PATH+mapa+'/domain_constraints.npy', np.array([z_min, z_max, x_min, x_max, y_min, y_max]))
except FileNotFoundError:
   os.mkdir(PATH+mapa)
   np. save(PATH+mapa+'/domain_constraints.npy', np.array([z_min, z_max, x_min, x_max, y_min, y_max]))
''' ....................................................................................................................................................................................'''


Z =                  z_max - z_min                           #  Number of cells along domain Z-axis;  unit: - [-]
X =                  x_max - x_min                           #  Number of cells along domain X-axis;  unit: - [-]
Y =                  y_max - y_min                           #  Number of cells along domain Y-axis;  unit: - [-]

old_grid = np.zeros((Z, X, Y))                        
new_grid = np.zeros((Z, X, Y))
rs_old = np.zeros((Z, X, Y))
rs_new = np.zeros((Z, X, Y))


""" ---------------------------------------------------------  *  F  U  N  C  T  I  O  N  S *  ----------------------------------------------------------------------------------------------------------"""
''' parameters of heterogeneous nucleation'''
ns_max =                      5e10                                                #  Number of new nuclei at optimal dT; unit: per SQUARE METER  [m^-2]
dTs_max0 =                       2                                                #   Optimal undercooling dT, where the formation of nuclei is the highest; unit:  KELVIN  [K]
dTs_sigma =                    0.5                                                #   Standard deviation of Gaussian nuclei distribution [K]

''' parameters of homogeneous nucleation '''
nv_max =                      5e14                                             
dTv_max0 =                       2
dTv_sigma =                    0.5 

def nukleacija_povrsina(pp):                                             # heterogeneous nucleation
   Ns_Nas =  ns_max/(math.sqrt(2*math.pi))*math.e**(-(pp-dTs_max0)**2/(2*dTs_sigma**2))
   #Ns_Nas = 0
   return Ns_Nas

def nukleacija_volumen(vv):                                              # homogeneous nucleation
   #Nv_Nav =  nv_max/(math.sqrt(2*math.pi))*math.e**(-(vv-dTv_max0)**2/(2*dTv_sigma**2))
   Nv_Nav=0
   return Nv_Nav

nukleusi_pov = np.vectorize(nukleacija_povrsina)
nukleusi_vol = np.vectorize(nukleacija_volumen)

def nakljucje(msm):
   rand = np.random.random_sample(msm.shape)
   rand[taula==-1]=-1
               
   return rand

def taljenje(u, Tm):
   global taula
   plt.clf()
   taula=np.zeros((Z,X,Y))
   taula[u<Tm] = -1
   return taula

def liquidus(u, liquidus_temp):
   global likvid
   likvid=np.zeros((Z,X,Y))
   likvid[u>liquidus_temp]=-1
   return likvid

# ****************************************** Grain orientation and Cube Rotation Matrices ******************************************

def Rotate_the_Cube_XY(alpha, beta):
    global RmX, RmY
    cube = np.array([
                              [1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1]
                            ])
    RmY = np.array([ # rotation matrix around Y-axis
        [math.cos(math.radians(alpha)), 0, math.sin(math.radians(alpha))],
        [0, 1, 0],
        [-math.sin(math.radians(alpha)), 0, math.cos(math.radians(alpha))],
        ])
    RmX = np.array([ # rotation matrix around X-axis
        [1, 0, 0],
        [0, math.cos(math.radians(beta)), math.sin(math.radians(beta))],
        [0, - math.sin(math.radians(beta)), math.cos(math.radians(beta))],
        ])
    oix = multi_dot([cube, RmY, RmX])
    return oix

def Rotate_the_Cube_Z(xy_cube, gamma):
    RmZ = np.array([# rotation matrix around Z-axis
        [math.cos(math.radians(gamma)), math.sin(math.radians(gamma)), 0],
        [- math.sin(math.radians(gamma)), math.cos(math.radians(gamma)), 0],
        [0, 0, 1],
        ])
    oiz = np.dot(xy_cube, RmZ)
    oiz_switched_rows_and_columns = oiz[np.ix_([2,0,1], [2,0,1])]
    return oiz_switched_rows_and_columns

def get_color(x,y):
    u = 1-x
    v = x-y
    w = y
    rgb = np.array([u,v,w])
    RGB = 255*rgb/np.max(rgb)
    return RGB.astype('int')

def random_xy_tilt():
    x = np.random.uniform(0,1)
    y = np.random.uniform(0,x)
    return x,y
   
''' ======================================================= NUCLEATION ======================================================= '''
asc = {}
grain_counter=0
def show_nuclei(interface, bulk, timesnap):
      global old_grid, new_grid, grain_counter, asc, live
      
      live= nakljucje(taula)
      time_step = 1
      ß = 0.1                              # larger ß will produce more grains, ß ranges from 0 to 10000 and beyond ..

      
      TEMP_now = np.load(PATH+mapa+'/salome_'+str(timesnap)+'.npy')[z_min:z_max, x_min:x_max, y_min:y_max]
      try:
         TEMP_next = np.load(PATH+mapa+'/salome_'+str(timesnap+time_step)+'.npy')[z_min:z_max, x_min:x_max, y_min:y_max]
      except IndexError:
         raise IndexError('Please, correct the time range!')
      
      for k in range(Z):
         for i in range(X):
            for j in range(Y):
               #if taula[k][i][j]==0 and old_grid[k][i][j]==0 and (live[k][i][j]<interface[k][i][j] or live[k][i][j]<bulk[k][i][j])and TEMP_next[k][i][j] > TEMP_now[k][i][j]:              
               if taula[k][i][j]==0 and old_grid[k][i][j]==0 and ((live[k][i][j]<interface[k][i][j] and interface[k][i][j]<ß) or (live[k][i][j]<bulk[k][i][j] and bulk[k][i][j]<ß))and TEMP_next[k][i][j] < TEMP_now[k][i][j]:              
                  grain_counter +=1
                  old_grid[k][i][j]=grain_counter
                  """ generation of random grain orientation """
                  x,y = random_xy_tilt()
                  rgb = get_color(x,y)
                  alfa  = math.degrees(math.atan(x))
                  beta = math.degrees(math.atan(y))- 9.7356103173*y
                  cub_xy = Rotate_the_Cube_XY(alfa, beta)
                  gama = math.degrees(math.atan(np.random.uniform(0,1)))
                  oi = Rotate_the_Cube_Z(cub_xy, gama)
                  asc[grain_counter] ={'oi': oi, 'alfa':alfa, 'beta':beta, 'gama':gama, 'rgb':rgb, 'coords':(timesnap,k,i,j), 'temp': TEMP_now[k,i,j]-273, }    # ALL data about nuclei
                  
      #condition = (taula==0) & (old_grid==0) & ((live/1 < interface) | (live/1 < bulk))& (TEMP_next < TEMP_now)  # concept of fast code
      #grain_counter+=1
      #np.save('jedra_dinamicno.npy', old_grid)
      
      return old_grid
				
def show_vg(kos):               # it shows the interface growth rate (vg)vs. undercooling temperatures (dT_field)at given time step
   T = np.load(PATH+mapa+'/salome_'+str(kos)+'.npy')[z_min:z_max, x_min:x_max, y_min:y_max]
   dT_field = T - Tmelt                                     # undecooling [K]
   vg=2.03e-4 * dT_field**2 - 0.544e-4 * dT_field
   plt.clf()
   plt.imshow(vg,  cmap=plt.get_cmap('rainbow'), interpolation='bessel')
   return vg

def get_growth_length(dT_field, cas):
   vg = 2.03e-4 * dT_field**2 - 0.544e-4 * dT_field
   length= vg*cas
   return length

def get_growth_speed(dT_field):
    vg = 2.03e-4 * dT_field**2 - 0.544e-4 * dT_field
    return vg

def merging(prva,druga):
    tot = np.where(druga==0, prva, druga)
    return tot

def grain_activity_check(seznam):
   pogoji=np.array([])
   for x in seznam:
      gac=np.array([x[0].size,x[1].size,x[2].size,])
      pog=np.all(gac==0)
      np.append(pogoji, pog)
   inactive = np.all(pogoji==True)
   return inactive 

   

def Growth_FAST_code(count, direction):                                                      #   One step of grain-growth for given grain in given direction
   global old_grid, new_grid, total, wij, vij, lij, rs_old, rs_new

   lij = vij*rs_new
   
   ''' 1st shell -------------------> six nearest neighbours '''
   if direction == '001':
      # ----> direction:  001    (1/6,      1/26)
      dij = np.array([0,0,1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      A001 = np.where((new_grid==0)& (taula==0)& (likvid==0)& (wij*lij >= L_1st_shell) &(np.pad(new_grid,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==count), count, new_grid) 
      rs001 = np.where((rs_new == 0) & (np.pad(rs_new,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==time_step), time_step, rs_new)  
      #np.save(mapa+'/meshes/001.npy', A001); del A001
      return A001, rs001
      '''
      try:
         total = merging(total, A001); del A001
      except NameError:
         total = merging(new_grid, A001); del A001
      '''

   elif direction == '00_1':
      # ----> direction:  00_1    (2/6,      2/26)
      dij = np.array([0,0,-1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      A00_1 = np.where((new_grid==0)& (taula==0)& (likvid==0)& (wij*lij >= L_1st_shell)&(np.pad(new_grid,((0,0),(0,0), (0,1)), 'constant')[:,:,1:] ==count), count, new_grid) 
      rs00_1 = np.where((new_grid == 0) & (np.pad(rs_new,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==time_step), time_step, rs_new)  
      #np.save(mapa+'/meshes/00_1.npy', A00_1); del A00_1
      return A00_1, rs00_1
      '''
      try:
         total = merging(total, A00_1); del A00_1
      except NameError:
         total = merging(new_grid, A00_1); del A00_1
      '''

   elif direction == '010':
      # ----> direction:  010    (3/6,      3/26)
      dij = np.array([0,1,0])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      A010 = np.where((new_grid==0)& (taula==0)& (likvid==0)& (wij*lij >= L_1st_shell)&(np.pad(new_grid,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:] ==count), count, new_grid) 
      rs010 = np.where((new_grid == 0) & (np.pad(rs_new,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==time_step), time_step, rs_new)  
      #np.save(mapa+'/meshes/010.npy', A010); del A010
      return A010, rs010
      '''
      try:
         total = merging(total, A010); del A010
      except NameError:
         total = merging(new_grid, A010); del A010
      '''

   elif direction == '0_10':
      # ----> direction:  0_10    (4/6,      4/26)
      dij = np.array([0,-1,0])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      A0_10 = np.where((new_grid==0)& (taula==0)& (likvid==0)& (wij*lij >= L_1st_shell)&(np.pad(new_grid,((0,0),(0,1), (0,0)), 'constant')[:,1:,:] ==count), count, new_grid) 
      rs0_10 = np.where((new_grid == 0) & (np.pad(rs_new,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==time_step), time_step, rs_new)  
      #np.save(mapa+'/meshes/0_10.npy', A0_10); del A0_10
      return A0_10, rs0_10
      '''
      try:
         total = merging(total, A0_10); del A0_10
      except NameError:
         total = merging(new_grid, A0_10); del A0_10
      '''

   elif direction == '100':
      # ----> direction:  100    (5/6,      5/26)
      dij = np.array([1,0,0])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      A100 = np.where((new_grid==0)& (taula==0)& (likvid==0)& (wij*lij >= L_1st_shell)&(np.pad(new_grid,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:] ==count), count, new_grid) 
      rs100 = np.where((new_grid == 0) & (np.pad(rs_new,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==time_step), time_step, rs_new)  
      #np.save(mapa+'/meshes/100.npy', A100); del A100
      return A100, rs100
      '''
      try:
         total = merging(total, A100); del A100
      except NameError:
         total = merging(new_grid, A100); del A100
      '''

   elif direction == '_100':
      # ----> direction:  _100    (6/6,      6/26)
      dij = np.array([-1,0,0])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      A_100 = np.where((new_grid==0)& (taula==0)& (likvid==0)& (wij*lij >= L_1st_shell)&(np.pad(new_grid,((0,1),(0,0), (0,0)), 'constant')[1:,:,:] ==count), count, new_grid) 
      rs_100 = np.where((new_grid == 0) & (np.pad(rs_new,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==time_step), time_step, rs_new)  
      #np.save(mapa+'/meshes/_100.npy', A_100); del A_100
      return A_100, rs_100
      '''
      try:
         total = merging(total, A_100); del A_100
      except NameError:
         total = merging(new_grid, A_100); del A_100
      '''

      ''' ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, END of 1st shell ~//~ '''

   
      ''' 2nd shell --------------> twelve second-nearest neighbours '''
   elif direction == '011':
      # ----> direction:  011    (1/12,      7/26)
      dij = np.array([0,1,1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      B011 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_2nd_shell) &(np.pad(new_grid,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==count), count, new_grid)
      rs011 = np.where((new_grid == 0) & (np.pad(rs_new,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==time_step), time_step, rs_new)
      #np.save(mapa+'/meshes/011.npy', B011); del B011
      return B011, rs011
      '''
      try:
         total = merging(total, B011); del B011
      except NameError:
         total = merging(new_grid, B011); del B011
      '''

   elif direction == '01_1':
      # ----> direction:  01_1    (2/12,      8/26)
      dij = np.array([0,1,-1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      B01_1 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_2nd_shell) &(np.pad(new_grid,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==count), count, new_grid)
      #np.save(mapa+'/meshes/01_1.npy', B01_1); del B01_1
      return B01_1
      '''
      try:
         total = merging(total, B01_1); del B01_1
      except NameError:
         total = merging(new_grid, B01_1); del B01_1
      '''

   elif direction == '0_11':
      # ----> direction:  0_11    (3/12,      9/26)
      dij = np.array([0,-1,1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      B0_11 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_2nd_shell) &(np.pad(new_grid,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==count), count, new_grid)
      #np.save(mapa+'/meshes/0_11.npy', B0_11); del B0_11
      return B0_11
      '''
      try:
         total = merging(total, B0_11); del B0_11
      except NameError:
         total = merging(new_grid, B0_11); del B0_11
      '''

   elif direction == '0_1_1':
      # ----> direction:  0_1_1    (4/12,      10/26)
      dij = np.array([0,-1,-1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      B0_1_1 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_2nd_shell) &(np.pad(new_grid,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==count), count, new_grid)
      #np.save(mapa+'/meshes/0_1_1.npy', B0_1_1); del B0_1_1
      return B0_1_1
      '''
      try:
         total = merging(total, B0_1_1); del B0_1_1
      except NameError:
         total = merging(new_grid, B0_1_1); del B0_1_1
      '''

   elif direction == '101':
      # ----> direction:  101    (5/12,      11/26)
      dij = np.array([1,0,1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      B101 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_2nd_shell) &(np.pad(new_grid,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==count), count, new_grid)
      #np.save(mapa+'/meshes/101.npy', B101); del B101
      return B101
      '''
      try:
         total = merging(total, B101); del B101
      except NameError:
         total = merging(new_grid, B101); del B101
      '''

   elif direction == '10_1':
      # ----> direction:  10_1    (6/12,      12/26)
      dij = np.array([1,0,-1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      B10_1 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_2nd_shell) &(np.pad(new_grid,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==count), count, new_grid)
      #np.save(mapa+'/meshes/10_1.npy', B10_1); del B10_1
      return B10_1
      '''
      try:
         total = merging(total, B10_1); del B10_1
      except NameError:
         total = merging(new_grid, B10_1); del B10_1
      '''

   elif direction == '110':
      # ----> direction:  110    (7/12,      13/26)
      dij = np.array([1,1,0])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      B110 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_2nd_shell) &(np.pad(new_grid,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==count), count, new_grid)
      #np.save(mapa+'/meshes/110.npy', B110); del B110
      return B110
      '''
      try:
         total = merging(total, B110); del B110
      except NameError:
         total = merging(new_grid, B110); del B110
      '''

   elif direction == '1_10':
      # ----> direction:  1_10    (8/12,      14/26)
      dij = np.array([1,-1,0])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      B1_10 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_2nd_shell) &(np.pad(new_grid,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==count), count, new_grid)
      #np.save(mapa+'/meshes/1_10.npy', B1_10); del B1_10
      return B1_10
      '''
      try:
         total = merging(total, B1_10); del B1_10
      except NameError:
         total = merging(new_grid, B1_10); del B1_10
      '''
         
   elif direction == '_101':
      # ----> direction:  _101    (9/12,      15/26)
      dij = np.array([-1,0,1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      B_101 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_2nd_shell) &(np.pad(new_grid,((0,1),(0,0), (1,0)), 'constant')[1:,:,:-1]==count), count, new_grid)
      #np.save(mapa+'/meshes/_101.npy', B_101); del B_101
      return B_101
      '''
      try:
         total = merging(total, B_101); del B_101
      except NameError:
         total = merging(new_grid, B_101); del B_101
      '''

   elif direction == '_10_1':
      # ----> direction:  _10_1    (10/12,      16/26)
      dij = np.array([-1,0,-1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      B_10_1 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_2nd_shell) &(np.pad(new_grid,((0,1),(0,0), (0,1)), 'constant')[1:,:,1:]==count), count, new_grid)
      #np.save(mapa+'/meshes/_10_1.npy', B_10_1); del B_10_1
      return B_10_1
      '''
      try:
         total = merging(total, B_10_1); del B_10_1
      except NameError:
         total = merging(new_grid, B_10_1); del B_10_1
      '''

   elif direction == '_110':
      # ----> direction:  _110    (11/12,      17/26)
      dij = np.array([-1,1,0])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      B_110 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_2nd_shell) &(np.pad(new_grid,((0,1),(1,0), (0,0)), 'constant')[1:,:-1,:]==count), count, new_grid)
      #np.save(mapa+'/meshes/_110.npy', B_110); del B_110
      return B_110
      '''
      try:
         total = merging(total, B_110); del B_110
      except NameError:
         total = merging(new_grid, B_110); del B_110
      '''

   elif direction == '_1_10':
      # ----> direction:  _1_10    (12/12,      18/26)
      dij = np.array([-1,-1,0])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      B_1_10 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_2nd_shell) &(np.pad(new_grid,((0,1),(0,1), (0,0)), 'constant')[1:,1:,:]==count), count, new_grid)
      #np.save(mapa+'/meshes/_1_10.npy', B_1_10); del B_1_10
      return B_1_10
      '''
      try:
         total = merging(total, B_1_10); del B_1_10
      except NameError:
         total = merging(new_grid, B_1_10); del B_1_10
      '''
      ''' ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, END of 2nd shell ~//~ '''


      ''' 3rd shell --------------> eight third-nearest neighbours '''
   elif direction == '111':
      # ----> direction: 111    (1/8,      19/26)
      dij = np.array([1,1,1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      C111 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_3rd_shell) &(np.pad(new_grid,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==count), count, new_grid)
      #np.save(mapa+'/meshes/111.npy', C111); del C111
      return C111
      '''
      try:
         total = merging(total, C111); del C111
      except NameError:
         total = merging(new_grid, C111); del C111
      '''

   elif direction == '11_1':
      # ----> direction: 11_1    (2/8,      20/26)
      dij = np.array([1,1,-1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      C11_1 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_3rd_shell) &(np.pad(new_grid,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==count), count, new_grid)
      #np.save(mapa+'/meshes/11_1.npy', C11_1); del C11_1
      return C11_1
      '''
      try:
         total = merging(total, C11_1); del C11_1
      except NameError:
         total = merging(new_grid, C11_1); del C11_1
      '''

   elif direction == '1_11':
      # ----> direction: 1_11    (3/8,      21/26)
      dij = np.array([1,-1,1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      C1_11 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_3rd_shell) &(np.pad(new_grid,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==count), count, new_grid)
      #np.save(mapa+'/meshes/1_11.npy', C1_11); del C1_11
      return C1_11
      '''
      try:
         total = merging(total, C1_11); del C1_11
      except NameError:
         total = merging(new_grid, C1_11); del C1_11
      '''

   elif direction == '1_1_1':
      # ----> direction: 1_1_1    (4/8,      22/26)
      dij = np.array([1,-1,-1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      C1_1_1 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_3rd_shell) &(np.pad(new_grid,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==count), count, new_grid)
      #np.save(mapa+'/meshes/1_1_1.npy', C1_1_1); del C1_1_1
      return C1_1_1
      '''
      try:
         total = merging(total, C1_1_1); del C1_1_1
      except NameError:
         total = merging(new_grid, C1_1_1); del C1_1_1
      '''

   elif direction == '_111':
      # ----> direction: _111    (5/8,      23/26)
      dij = np.array([-1,1,1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      C_111 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_3rd_shell) &(np.pad(new_grid,((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==count), count, new_grid)
      #np.save(mapa+'/meshes/_111.npy', C_111); del C_111
      return C_111
      '''
      try:
         total = merging(total, C_111); del C_111
      except NameError:
         total = merging(new_grid, C_111); del C_111
      '''

   elif direction == '_11_1':
      # ----> direction: _11_1    (6/8,      24/26)
      dij = np.array([-1,1,-1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      C_11_1 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_3rd_shell) &(np.pad(new_grid,((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==count), count, new_grid)
      #np.save(mapa+'/meshes/_11_1.npy', C_11_1); del C_11_1
      return C_11_1
      '''
      try:
         total = merging(total, C_11_1); del C_11_1
      except NameError:
         total = merging(new_grid, C_11_1); del C_11_1
      '''

   elif direction == '_1_11':
      # ----> direction: _1_11    (7/8,      25/26)
      dij = np.array([-1,-1,1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      C_1_11 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_3rd_shell) &(np.pad(new_grid,((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==count), count, new_grid)
      #np.save(mapa+'/meshes/_1_11.npy', C_1_11); del C_1_11
      return C_1_11
      '''
      try:
         total = merging(total, C_1_11); del C_1_11
      except NameError:
         total = merging(new_grid, C_1_11); del C_1_11
      '''

   elif direction == '_1_1_1':
      # ----> direction: _1_1_1    (8/8,      26/26)
      dij = np.array([-1,-1,-1])
      wij = np.max(np.linalg.norm(asc[count]['oi'] * dij, axis=1))
      C_1_1_1 = np.where((new_grid==0)& (taula==0)& (wij*lij >= L_3rd_shell) &(np.pad(new_grid,((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==count), count, new_grid)
      #np.save(mapa+'/meshes/_1_1_1.npy', C_1_1_1); del C_1_1_1
      return C_1_1_1
      '''
      try:
         total = merging(total, C_1_1_1); del C_1_1_1
      except NameError:
         total = merging(new_grid, C_1_1_1); del C_1_1_1
      ''' 
   ''' ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, END of 3rd shell ~//~ '''


   #new_grid = total.copy()         
   #old_grid=new_grid.copy()
   
   moves = [A001, A010, A00_1, A0_10, A100, A_100, ]
   grain_inactive = grain_activity_check(moves)
   if grain_inactive==True:
      del asc[count]
      print(f'Grain number {count} is no longer active..')
   

def random_selection(x):
    global item
    if len(x)> 0:
        item = random.choice(x)
        x=np.delete(x,np.where(x==item))
    else:
        item = x[0]
        x=np.delete(x,np.where(x==item))
    return x

""" ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""

avtomatska_nukleacija = False
rast_zrn = True

mn = {1:{'ß':(0, 50, 50), 'Ł': (0,0,0)},
           2: {'ß':(0, 150, 150), 'Ł': (0,0,45)},                    #  data of manually created nuclei, 'ß' are the (Z,X,Y) coordinates, 'Ł' are tilting parameters (x,y,gama)
           #3:{'ß':(1, 10, 20), 'Ł': (0,0,0)},
           #4: {'ß':(1, 0, 40), 'Ł': (1,1,20)},
      
         }

if not avtomatska_nukleacija:
   '''||||||||||||||||||||||||||| nucleation - manual |||||||||||||||||||||||||| '''
   for i in mn:
      old_grid[mn[i]['ß'][0],mn[i]['ß'][1],mn[i]['ß'][2]]=i; new_grid[mn[i]['ß'][0],mn[i]['ß'][1],mn[i]['ß'][2]]=i
      rs_old[new_grid != 0] += time_step; rs_new[new_grid != 0] += time_step
      #x,y = random_xy_tilt()
      x,y = mn[i]['Ł'][0], mn[i]['Ł'][1] 
      rgb = get_color(x,y)
      alfa  = math.degrees(math.atan(x))
      beta = math.degrees(math.atan(y))- 9.7356103173*y
      cub_xy = Rotate_the_Cube_XY(alfa, beta)
      #gama = math.degrees(math.atan(np.random.uniform(0,1)))
      gama = mn[i]['Ł'][2]
      oi = Rotate_the_Cube_Z(cub_xy, gama)
      asc[i] ={'oi': oi, 'alfa':alfa, 'beta':beta, 'gama':gama, 'rgb':rgb,}
   ''' ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| '''


T = domena
S = ([0,0,1], [0,1,0], [0,0,-1], [0,-1,0], [1,0,0], [-1,0,0], )
selected_grains = [1]

def Lij(arr, zrno):
    f1stdij = np.array(arr)
    f1stwij = round(np.max(np.linalg.norm(asc[zrno]['oi'] * f1stdij, axis=1)), 3)
    f1stlij = round(f1stwij*lij[old_grid==zrno][0]*1e06, 1);   share_f1stlij = round(100*f1stwij*lij[old_grid==zrno][0]/cell, 2)
    return f1stwij, f1stlij, share_f1stlij

def Mesh_Adopted_Growth(cas, real_step, smeri, growth, auto_nucleation):
   global new_grid, old_grid, vij, ascar, wij, rs_old, rs_new, tt
   #rs=real_step
   #rs=0.00385
   #rs=0.004725
   
   for tt in range(cas[0], cas[1]-1):
      step_time_start = time.time()
      #T = np.load(PATH+mapa+'/salome_'+str(tt)+'.npy')[z_min:z_max, x_min:x_max, y_min:y_max]
      
      taljenje(T, Tmelt)                              # condition for melting  ----> taula matrix of cell values;  if value -1 then solid (or powder)elif value 0 then liquid
      liquidus(T,Tmelt+dTliquidus)              # absolute liquidus line
      
      dTt = T - Tmelt                                     # undecooling [K]

      if auto_nucleation:
         nuk_matrix_pov=nukleusi_pov(dTt)                                              
         nuk_matrix_vol=nukleusi_vol(dTt)
         show_nuclei(nuk_matrix_pov, nuk_matrix_vol, tt)           #   formation of new nuclei 
      
      new_grid = old_grid.copy()
      if growth:
         #rs[new_grid != 0] += real_step
         #lij = get_growth_length(dTt, rs)#; rs+=real_step
         vij = get_growth_speed(dTt)
         
         u=0
         ascar = np.array(list(asc.keys()))
         for grain in ascar:
            ascar = random_selection(ascar)
            random_grain = item
            directions = smeri.copy()
            seznam_premikov = []
            for s in directions:
               directions = random_selection(directions)
               random_direction = item
               q=Growth_FAST_code(random_grain, random_direction)
               seznam_premikov.append(q)

            #*********** MERGING the MESHES *************#
            '''
            mesh = 0
            total = merging(np.load(mapa+'/meshes/'+smeri[mesh]+'.npy'), np.load(mapa+'/meshes/'+smeri[mesh+1]+'.npy'))
            for mesh in range(2, len(smeri)):
               total = merging(total, np.load(mapa+'/meshes/'+smeri[mesh]+'.npy'))
            '''
            total = merging(seznam_premikov[0][0], seznam_premikov[1][0])
            total_rs = merging(seznam_premikov[0][1], seznam_premikov[1][1])
            for mesh in range(2, len(smeri)):
               total = merging(total, seznam_premikov[mesh][0])
               total_rs = merging(total_rs, seznam_premikov[mesh][1])
            
            new_grid = total.copy()
            rs_new = total_rs.copy()
            old_grid = new_grid.copy()
            rs_old = rs_new.copy()
            #rs_old[new_grid != 0] += time_step; rs_new[new_grid != 0] += time_step
      #rs_old[rs_old != 0] += time_step; rs_new[rs_new != 0] += time_step
      
      
         
      flashy_snap = np.zeros((old_grid.shape[0], old_grid.shape[1], old_grid.shape[2], 3)).astype('int')
      for zrno in asc:
         flashy_snap[old_grid==zrno]=asc[zrno]['rgb']
      try:
         np.save(PATH+mapa+'/'+flashies+'/flashy_snap_'+str(tt)+'.npy', flashy_snap)
      except FileNotFoundError:
         os.mkdir(PATH+mapa+'/'+flashies)
         np.save(PATH+mapa+'/'+flashies+'/flashy_snap_'+str(tt)+'.npy', flashy_snap)

      step_cpu = round(time.time()- step_time_start, 1)
      print('Časovni korak: ',tt,'............ čas računanja: ',step_cpu,' sekund........... št. zrn: ',grain_counter)
      try:
         H=[]
         for sg in selected_grains:
             lij_database=[]
             for s in S:
                lij_database.append(Lij(s, sg))
             H.append((sg, lij_database))

      except (IndexError, KeyError):
          f1stlij = 'N.A.'; share_f1stlij='N.A.'
          print('we have index error')
      except NameError:
          print('wij does not exist')
          f1stlij = 'N.A.'; share_f1stlij='N.A.'
          
      #print('Real Step (rs): ',round(rs[0][0], 6), ' sekund.................')
      for h in H:
          print('Velikost zrna #',h[0],':   ',np.count_nonzero(old_grid==h[0]),'       Orientacija=(',u'\u03B1=',asc[h[0]]['alfa'],u'\u00B0,   ', u'\u03B2=',asc[h[0]]['beta'],u'\u00B0,   ', u'\u03B3=',asc[h[0]]['gama'],u'\u00B0)')
          print('orientation weight coefficient (Wij),  growth length grain (Lij=Wij*vg*rs),  Lij share of cell-size (%)')
          for sm, hh in enumerate(h[1]):
              word=['      <'+str(S[sm])[1:-1]+'>','    |',"{:.3f}".format(hh[0]),'     |', str(hh[1])+u' \u00B5m           |', str(hh[2])+' % of cell-size']
              print('{:>16}{:>10}  {:>15}{:>12}  {:>28} {:>30}'.format(word[0], word[1], word[2], word[3], word[4], word[5]))

          print()
      print()

   #with open(PATH+mapa+'/nuclei_data.json', 'w') as nuks:              # Writing data of ALL nuclei as .json file
    #  json.dump(asc, nuks)

   elapsed_time = time.time() - start_time
   print(185*'*')
   print(50*' ','DONE!!  Total computing time =  ',round(elapsed_time, 3),'  seconds.')
   print(grain_counter,' grains were formed!'); print(185*'*'); print()
   #np.save(mapa+'/'+result, old_grid)
   
         
start_time = time.time()
Mesh_Adopted_Growth(time_range, time_step, orientacije, rast_zrn, avtomatska_nukleacija)

#alfas=[asc[i]['alfa'] for i in asc]
#betas=[asc[i]['beta'] for i in asc]


      
