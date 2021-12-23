"""
Program za nukleacijo in rast zrn po 'Cellular Automata Mesh Dependency' metodi v 3D,
avtor:  dr. Andraž Kocjan
          Inštitut za kovinske materiale in tehnologije
          Oktober 2020
"""
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import multi_dot
import random, math, time
plt.ion()

""" =====================================  F  U  N  C  T  I  O  N  S ================================================================"""
# **** Domain Size Reduction ****
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

# **** Nucleation functions ****
def nukleacija_povrsina(pp):                                                  # heterogeneous nucleation (at the liquid-solid interface ---> melt-pool border)
    ''' parameters of heterogeneous nucleation'''
    ns_max =                      5e10                                                #  Number of new nuclei at optimal dT; unit: per SQUARE METER  [m^-2]
    dTs_max0 =                       2                                                 #   Optimal undercooling dT, where the formation of nuclei is the highest; unit:  KELVIN  [K]
    dTs_sigma =                    0.5                                                 #   Standard deviation of Gaussian nuclei distribution [K]
    Ns_Nas =  ns_max/(math.sqrt(2*math.pi))*math.e**(-(pp-dTs_max0)**2/(2*dTs_sigma**2))
    #Ns_Nas = 0
    return Ns_Nas

def nukleacija_volumen(vv):                                                   # homogeneous nucleation
    ''' parameters of homogeneous nucleation '''
    nv_max =                      5e14                                             
    dTv_max0 =                       2
    dTv_sigma =                    0.5 
    Nv_Nav =  nv_max/(math.sqrt(2*math.pi))*math.e**(-(vv-dTv_max0)**2/(2*dTv_sigma**2))
    Nv_Nav=0
    return Nv_Nav

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

# **** Growth functions ****
def liquidus(u, liquidus_temp):
   global likvid
   likvid=np.zeros((Z,X,Y))
   likvid[u>liquidus_temp]=-1
   return likvid

def growth_speed(dT_field):
    vg = 2.03e-4 * dT_field**2 - 0.544e-4 * dT_field
    return vg


# **** Grain orientation and Cube Rotation Matrices ****

def Rotate_the_Cube_XY(alpha, beta):              # rotation of cube around Y- axis (alpha angle) and  X-axis (beta angle)
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

def Rotate_the_Cube_Z(xy_cube, gamma):                   # rotation of cube around Z- axis (gamma angle)
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

def random_xy_tilt(posibilities):
    #x = np.random.uniform(0,1)
    #y = np.random.uniform(0,1)   
    xran = np.random.randint(0, posibilities+1)
    x = xran / posibilities
    y = np.random.randint(0, xran+1)/posibilities
    return x,y

def merging(prva,druga):
    tot = np.where(druga==0, prva, druga)
    return tot

def time_counter(neg_index):
   Negatives[neg_index] += dt

def random_selection(x):
    global item
    if len(x)> 0:
        item = random.choice(x)
        x=np.delete(x,np.where(x==item))
    else:
        item = x[0]
        x=np.delete(x,np.where(x==item))
    return x

def m(x):
    mic = np.load(PATH+mapa+'/'+flashies+'/flashy_snap_'+str(x)+'.npy')
    #plt.title(u'\u03B1 = '+a)
    plt.imshow(mic[13])

def TT(x):
   temp = np.load(PATH+mapa+'/salome_'+str(x)+'.npy')[z_min:z_max, x_min:x_max, y_min:y_max]
   return temp

def W(O, L):
    O = np.array(O)
    cos_Z = np.dot(O[0], L) / (np.linalg.norm(O[0]* np.linalg.norm(L)))
    cos_X = np.dot(O[1], L) / (np.linalg.norm(O[1] * np.linalg.norm(L)))
    cos_Y = np.dot(O[2], L) / (np.linalg.norm(O[2] * np.linalg.norm(L)))
    return np.max(np.absolute(np.array([cos_Z, cos_X, cos_Y])))

def Dissipate_Distances(eps):
    R_min  = 0.5 - epsilon
    R_max = 1-R_min
    RAN = np.random.uniform(R_min, R_max, (Z,X,Y))
    return RAN


""" ====================================  F  I  L  E  S    &   F  O  L  D  E  R  S ========================================================"""
# Path to FE analysis results from Salome Meca 2018

#case =       'mesh_dependency_experiment'      # special case (for development purpose)
#case =     'cooling_plate_1650'
#case =     '3D_ss800_f5R8c4f8b23d1N30_10W_t5'
#case =     '3D_ss800_f5R8c4f8b23d2N30_10W_t5'
#case =       'rp1_ss800_f5R8c4f8b23N30_7.5W_t5'
#case = 'salome_4D'
#case =   'SLM_2D_Source'
#PATH = 'C:/sm-2018-w64-0-3/WORK/kode za Salome/Cellular Automata/'+case+'/'


case =      'SLM_2D_Source'
subcase = '0002'

#case =       'cooling_plate'
#subcase =  '0001'


PATH = 'C:/sm-2018-w64-0-3/WORK/'+case+'_Files/post processing database/'+subcase+'/'

#mapa    =                       'INTER_40_8  Z[7-8], X[18-24], Y[12-97], 1500°C'           #  Folder name (should be created manually, if it doesn't exist)
mapa    =                       'INTER  time=40, space=8  Z[7-8], X[15-27], Y[12-97], 1500°C, N=12' 


flashies =                       'flashies_001'                                                              #  Subfolder with time-snap 3D matrices, i.e. flashies

'''
z_min =         7
z_max =        8
Z =                   z_max - z_min

x_min =        
x_max =
X=                  x_max - x_min

y_min =
y_max =
Y=                 y_max - y_min
'''


""" =====================================  I  N  P  U  T    D  A  T  A  ======================================================================="""
# ````` DOMAIN Input variables `````

avtomatska_nukleacija   =     False
forced_nucleation =  False                          #  It counts only if avtomatska_nukleacija == True

reduce_domain_size     =     False
delete_inactive_grains_ID = False

save_flashy_as_RGB = False


cell=1              # for mesh dependency development (MDD)
#dt=cell/4
dt=cell

#dt =               (0.03/8)/ 14              #  for cooling_plate   (CA_edge_680_8)
#dt =               (0.03/8)/36.5            #  for cooling_plate    (CA_edge_48_8)
#dt =                  0.03/8                      # for cooling_plate     (CA_1_8)
#cell =                5e-05/8                    # for cooling plate
#dt =                   cell

#cell =               5e-06/4                     # for SLM                                                      # Size of cubic cell (here considered as voxel), a.u.
#dt   =               5e-06/4                     #  for SLM (3D source)                                 # Size of time-step, a.u.

#cell =              5e-06/4                      # for SLM (2D source)
#dt =                6.25e-6/4                   #  for SLM (2D source)

#cell =                 5e-06/8                    # case: SLM_2D_Source
#dt =                   cell/3                       # case: SLM_2D_Source


'''**************************'''
START_step =         0                       # Starting time step number
END_step     =         90                    # Ending time step number
'''**************************'''


smeri = np.array(['001', '00_1', '010', '0_10',  ]) #'100', '_100',                                                                                                #   1st shell (6-nearest neighbours)---> grain-growth directions with respect to domain axes('ZXY')
                           #'011', '01_1', '0_11', '0_1_1', #])#'101', '110', '10_1', '1_10', '_101', '_110', '_10_1', '_1_10',  #])         #   2nd shell (12-nearest neighbours)directions..                                            
                           #'111', '11_1', '1_11', '1_1_1', '_111', '_11_1', '_1_11', '_1_1_1',   ])                                            #   3rd shell (8-nearest neighbours)directions..
                           #'012', '01_2', '0_12', '0_1_2',    '021', '02_1', '0_21', '0_2_1',  #])                                                           #   Extra neighbours  (22.5 degrees direction)
                           #'002', '00_2', '020', '0_20', #])
                           #'022', '02_2', '0_22', '0_2_2', ])

# `````` Material properties ``````

''' SS 316 at melting point, i.e. 1650 Kelvin '''
Lambda = 35
Rho = 7284
Cp = 678
TherDiff = Lambda / (Rho*Cp)
DTime = cell**2 / TherDiff

#Tmelt_Celsius =                          1643             #     for cooling_plate_1650               #   Melting point; unit:  Deg.  CELSIUS  [deg. C]
Tmelt_Celsius =                            1007             #     for SLM
 
#dTliquidus =                 30        # for cooling_plate                                                       #   Difference between Absolute Liquidus Line and Melting point; unit: Deg. CELSIUS [deg. C]
dTliquidus =                   80        # for SLM ; salome_4D


# ````` Supplemental variables, functions and calculated parameters `````
Negatives = {-1:0}; asc = {}; grain_counter = 0; S = {}; np.random.seed(75689803)
Tmelt= Tmelt_Celsius + 273         #   Melting point; unit:  KELVIN [K]
rp = 100                                       #  Number of possible random alfa, beta, gama choices for grain orientation randomization

''' ..........................NEIGHBOURHOOD  DISTANCES   (fixed or randomized)......................... '''

'''********************'''
epsilon =   0
'''********************'''

if epsilon == 0:
   random_distances = False
elif epsilon > 0 and epsilon <= 0.49:
   random_distances = True


def Dsr_1st(r):
   if random_distances:
      cell_1st = cell * (0.5 + 1.158*r)                              #"""------------------1st shell ::: [001], [010], [100] ::: 6 neighbours ------------------"""
   else:
      cell_1st=cell
   return cell_1st


def Dsr_2nd(r):
   if random_distances:
      cell_2nd=cell*math.sqrt(2)*(0.5 + 1.081*r)           #''' ------------------2nd shell ::: [011] ::: 4 neighbours ------------------ '''
   else:
      cell_2nd=cell*math.sqrt(2)
   return cell_2nd

def Dsr_2nd_II(r):
   if random_distances:
      cell_2nd_II = None    # to be constructed                  #''' ------------------2nd shell_II ::: [101], [110] ::: 8 neighbours ------------------ '''
   else:
      cell_2nd_II = cell*math.sqrt(2)
   return cell_2nd_II


def Dsr_3rd(r):
   if random_distances:
      cell_3rd= None         # !!! to be constructed              #''' ------------------3rd shell ::: [111] ::: 8 neighbours ------------------ '''
   else:
      cell_3rd=cell*math.sqrt(3)                                      
   return cell_3rd


def Dsr_4th(r):
   if random_distances:
      cell_4th=cell*(1.581 + 1.376*r)                                #''' ------------------4th shell ::: [012], [021] ::: 8 neighbours ------------------ '''
   else:
      cell_4th=cell*math.sqrt(5)
   return cell_4th

def Dsr_5th(r):
   if random_distances:
      cell_5th=cell*(1.5 + 1.098*r)                                    #''' ------------------5th shell ::: [002], [020], [200] ::: 6 neighbours ------------------ '''
   else:
      cell_5th=cell*2
   return cell_5th

def Dsr_6th(r):
   if random_distances:
      cell_6th=cell*(2.12132 + 1.44939*r)                          #''' ------------------6th shell ::: [022] ::: 4 neighbours ------------------ '''
   else:
      cell_6th=cell*math.sqrt(8)
   return cell_6th

   

# ````` Automatic determination of system size due to length, width and depth of melt-pool at starting and final time step `````

if avtomatska_nukleacija:
   if reduce_domain_size:
      limits = [ [], [], [], [], [], [] ]
      for i in (START_step, END_step):
         Domain_Size_Reduction(np.load(PATH+mapa+'/'+'salome_'+str(i)+'.npy'), Tmelt)
         limits[0].append(z_min); limits[1].append(z_max)
         limits[2].append(x_min); limits[3].append(x_max)
         limits[4].append(y_min); limits[5].append(y_max)

      z_min, z_max, x_min, x_max, y_min, y_max = min(limits[0]), max(limits[1]), min(limits[2]), max(limits[3]), min(limits[4]), max(limits[5])
      z_min = 0
      z_max = 1

   else:
      domena = np.load(PATH+mapa+'/'+'salome_'+str(START_step)+'.npy')
      z_min =                           7
      z_max =                          8#     domena.shape[0]

      x_min =                           0
      x_max =                                 domena.shape[1]

      y_min =                           0
      y_max =                                 domena.shape[2]

   ''' ........................................... domain constraints - FINAL and Saving as .Npy file ..............................................................'''
   #z_min=11
   #z_max=12
   try:
      np. save(PATH+mapa+'/domain_constraints.npy', np.array([z_min, z_max, x_min, x_max, y_min, y_max]))
   except FileNotFoundError:
      os.mkdir(PATH+mapa)
      np. save(PATH+mapa+'/domain_constraints.npy', np.array([z_min, z_max, x_min, x_max, y_min, y_max]))
   ''' ....................................................................................................................................................................................'''

   Z =                  z_max - z_min                           #  Number of cells along domain Z-axis;  unit: - [-]
   X =                  x_max - x_min                           #  Number of cells along domain X-axis;  unit: - [-]
   Y =                  y_max - y_min                           #  Number of cells along domain Y-axis;  unit: - [-]


if avtomatska_nukleacija:
   faza = np.zeros((Z,X,Y))                               # fazna matrika
   cas = np.zeros((Z,X,Y))                                # časovna matrika 
   vg = np.zeros((Z,X,Y))                                  # matrika hitrosti rasti
   NP = np.vectorize(nukleacija_povrsina)
   NV = np.vectorize(nukleacija_volumen)


''' ======================================= N  U  C  L  E  A  T  I  O  N ============================================================== '''
#  ~~~~ Manual Nucleation (MN)~~~~

if not avtomatska_nukleacija:
    '''||||||||||||||||||||||||||| nucleation - manual |||||||||||||||||||||||||| '''
    Z,X,Y = 1, 201, 201                                              # Size of domain in terms of cells in Z,X,Y directions, respectively, for testing and code development
    faza = np.zeros((Z,X,Y))                               # fazna matrika
    cas = np.zeros((Z,X,Y))                                # časovna matrika
    vg = np.zeros((Z,X,Y))                                  # matrika hitrosti rasti
    vg =   1                                                               # Value of homogeneous 'growing velocity' field, a.u., for testing and code development


    M = { 1: {'ß':(0, 100, 100), 'Ł': (0,0,0)},
               #2: {'ß':(0, 51, 51), 'Ł': (0,0,0)},
               #3: {'ß':(0, 45, 45), 'Ł': (0.5,0,22.5)},                    #  data of manually created nuclei, 'ß' are the (Z,X,Y) coordinates, 'Ł' are tilting parameters (x,y,gama)
               #4: {'ß':(0, 150, 150), 'Ł': (0,0,45*0.5)},
               #5: {'ß':(0, 125, 125), 'Ł': (0,0,45*0.75)},
               #6: {'ß':(0, 150, 150), 'Ł': (0,0,45*0.875)},
               #7: {'ß':(0, 250, 250), 'Ł': (0,0,45)},

               #8: {'ß':(0, 3, 3), 'Ł': (0,0,0)},
             }

    for i in M:
        faza[M[i]['ß'][0],M[i]['ß'][1],M[i]['ß'][2]]=i                               # define nucleus ID in faza matrix
        x,y = M[i]['Ł'][0], M[i]['Ł'][1]
        rgb = get_color(x,y)
        alfa  = math.degrees(math.atan(x))- 9.7356103173*x*y
        beta = math.degrees(math.atan(y))- 9.7356103173*x*y
        cub_xy = Rotate_the_Cube_XY(alfa, beta)
        gama = M[i]['Ł'][2]
        oi = Rotate_the_Cube_Z(cub_xy, gama)
        asc[i] ={'oi': oi, 'alfa':alfa, 'beta':beta, 'gama':gama, 'rgb':rgb,}

    grain_ID = np.array(list(asc.keys()))
    cas[np.isin(faza, grain_ID, invert=False)] = -1
    taula=0;likvid=0; grain_counter=len(grain_ID)

   
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Growth (time FOR-LOOP, grains FOR-LOOP, directions FOR-LOOP)..
start_time = time.time(); inactive_grains=[]

if random_distances:
       R=Dissipate_Distances(epsilon)
else:
       R=1

fn =     1                                                               # Additional weight (factor) of orientation weight (W) for first neighbours (fn)
sn =    1                                                               # Additional weight (factor) of orientation weight (W) for second neighbours (sn)
en =    1  #1/math.sqrt(2)                                   # Additional weight (factor) of orientation weight (W) for extra neighbours (en)
en2 =  1
en3 =  1

IG = {}        #  IG (saves time-steps when indivudual grain didn't grow IN ANY DIRECTION, IG stands for Inactive Grains)
AG = {}       #  AG (saves time-steps when indivudual grain did grow IN ANY DIRECTION, AG stands for Active Grains)

for i in range(START_step, END_step):
    step_time_start = time.time()
    # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   1.   N  U  C  L  E  A  T  I  O  N   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
     
    if avtomatska_nukleacija:
        ''' avtomatska nukleacija '''
        T = np.load(PATH+mapa+'/salome_'+str(i)+'.npy')[z_min:z_max, x_min:x_max, y_min:y_max]        # temperaturno polje izbranega časovnega koraka
        taljenje(T, Tmelt)                              # condition for melting  ----> taula matrix of cell values;  if value -1 then solid (or powder)elif value 0 then liquid
        liquidus(T,Tmelt+dTliquidus)              # absolute liquidus line
      
        dTt = T - Tmelt                                     # undecooling [K]
        
        interface = NP(dTt)                                              
        bulk = NV(dTt)
        live= nakljucje(taula)
        time_shift = 1

        #ß = 1000              # for SLM without mesh dependency                                               # larger ß will produce more grains, ß ranges from 0 to 10000 and beyond ..
        ß = 0.1                  # case: SLM_2D_Source(cell/3)

        #ß = 0.001     # for cooling_plate

        try:
            T_next = np.load(PATH+mapa+'/salome_'+str(i+time_shift)+'.npy')[z_min:z_max, x_min:x_max, y_min:y_max]
        except IndexError:
            raise IndexError('Please, correct the time range!')

        '''............. the following code can be written with numpy broadcasting to avoid outrageously slow for-loops.. .....................................'''
        new_grains_ID = []
        if not forced_nucleation:
           for k in range(Z):
            for ii in range(X):
               for j in range(Y):          
                  if taula[k][ii][j]==0 and faza[k][ii][j]==0 and ((live[k][ii][j]<interface[k][ii][j] and interface[k][ii][j]<ß) or (live[k][ii][j]<bulk[k][ii][j] and bulk[k][ii][j]<ß))and T_next[k][ii][j] < T[k][ii][j]:              
                     grain_counter +=1
                     new_grains_ID.append(grain_counter)
                     IG[grain_counter]=[]
                     AG[grain_counter]=[]
                     faza[k][ii][j]=grain_counter
                     """ generation of random grain orientation """
                     x,y = random_xy_tilt(rp)
                     rgb = get_color(x,y)
                     alfa  = math.degrees(math.atan(x))- 9.7356103173*x*y
                     beta = math.degrees(math.atan(y))- 9.7356103173*x*y
                     cub_xy = Rotate_the_Cube_XY(alfa, beta)
                     gama = math.degrees(math.atan(np.random.randint(0, rp+1)/rp))
                     #gama = math.degrees(math.atan(np.random.uniform(0,1)))
                     oi = Rotate_the_Cube_Z(cub_xy, gama)
                     asc[grain_counter] ={'oi': oi, 'alfa':alfa, 'beta':beta, 'gama':gama, 'rgb':rgb, 'coords':(i,k,ii,j), 'temp': T[k,ii,j]-273, }    # ALL data about nuclei
        
        if i==START_step and forced_nucleation:
           '''forced nucleation - at initial time step'''
           grain_counter +=1
           new_grains_ID.append(grain_counter)
           IG[grain_counter]=[]
           AG[grain_counter]=[]
           k=       0      # Z
           ii=     40      # X
           j=      14      # Y
           faza[k][ii][j]=grain_counter
           """ manual grain orientation setup """
           x =      0          
           y =      0            
           #x,y = random_xy_tilt(rp)
           rgb = get_color(x,y)
           #alfa  = math.degrees(math.atan(x))- 9.7356103173*x*y
           #beta = math.degrees(math.atan(y))- 9.7356103173*x*y
           cub_xy = Rotate_the_Cube_XY(alfa, beta)
           gama = 45    #math.degrees(math.atan(np.random.randint(0, rp+1)/rp))
           oi = Rotate_the_Cube_Z(cub_xy, gama)
           asc[grain_counter] ={'oi': oi, 'alfa':alfa, 'beta':beta, 'gama':gama, 'rgb':rgb, 'coords':(i,k,ii,j), 'temp': T[k,ii,j]-273, }
           
        grain_ID = np.array(list(asc.keys()))
        cas[np.isin(faza, new_grains_ID, invert=False)] = -1          # vrednost časovne matrike vseh novih nukleusov je -1
        vg = growth_speed(dTt)
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 

    # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   2.   T  I  M  E   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    for negind in Negatives:
        time_counter(negind)

    for s in smeri[:]:
        for negind in Negatives:
            ''' ----------------------------------------- 1st shell -------------------> 6 nearest neighbours '''
            if s == '001':
                if negind == -1:
                    cas001 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==negind), Negatives[negind], cas)
                    S[s]=cas001; del cas001
                else:
                    cas001 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==negind), Negatives[negind], S[s])
                    S[s]=cas001; del cas001

            elif s == '00_1':
                if negind == -1:
                    cas00_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==negind), Negatives[negind], cas)
                    S[s]=cas00_1; del cas00_1
                else:
                    cas00_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==negind), Negatives[negind], S[s])
                    S[s]=cas00_1; del cas00_1

            elif s == '010':
                if negind == -1:
                    cas010 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==negind), Negatives[negind], cas)
                    S[s]=cas010; del cas010
                else:
                    cas010 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==negind), Negatives[negind], S[s])
                    S[s]=cas010; del cas010

            elif s == '0_10':
                if negind == -1:
                    cas0_10 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==negind), Negatives[negind], cas)
                    S[s]=cas0_10; del cas0_10
                else:
                    cas0_10 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==negind), Negatives[negind], S[s])
                    S[s]=cas0_10; del cas0_10

            elif s == '100':
                if negind == -1:
                    cas100 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==negind), Negatives[negind], cas)
                    S[s]=cas100; del cas100
                else:
                    cas100 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==negind), Negatives[negind], S[s])
                    S[s]=cas100; del cas100

            elif s == '_100':
                if negind == -1:
                    cas_100 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==negind), Negatives[negind], cas)
                    S[s]=cas_100; del cas_100
                else:
                    cas_100 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==negind), Negatives[negind], S[s])
                    S[s]=cas_100; del cas_100
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

                    ''' ----------------------------------------- 2nd shell -------------------> 12 nearest neighbours '''

            elif s == '011':
                if negind == -1:
                    cas011 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==negind), Negatives[negind], cas)
                    S[s]=cas011; del cas011
                else:
                    cas011 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==negind), Negatives[negind], S[s])
                    S[s]=cas011; del cas011

            elif s == '01_1':
                if negind == -1:
                    cas01_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==negind), Negatives[negind], cas)
                    S[s]=cas01_1; del cas01_1
                else:
                    cas01_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==negind), Negatives[negind], S[s])
                    S[s]=cas01_1; del cas01_1

            elif s == '0_11':
                if negind == -1:
                    cas0_11 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==negind), Negatives[negind], cas)
                    S[s]=cas0_11; del cas0_11
                else:
                    cas0_11 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==negind), Negatives[negind], S[s])
                    S[s]=cas0_11; del cas0_11

            elif s == '0_1_1':
                if negind == -1:
                    cas0_1_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==negind), Negatives[negind], cas)
                    S[s]=cas0_1_1; del cas0_1_1
                else:
                    cas0_1_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==negind), Negatives[negind], S[s])
                    S[s]=cas0_1_1; del cas0_1_1

            elif s == '101':
                if negind == -1:
                    cas101 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==negind), Negatives[negind], cas)
                    S[s]=cas101; del cas101
                else:
                    cas101 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==negind), Negatives[negind], S[s])
                    S[s]=cas101; del cas101

            elif s == '10_1':
                if negind == -1:
                    cas10_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==negind), Negatives[negind], cas)
                    S[s]=cas10_1; del cas10_1
                else:
                    cas10_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==negind), Negatives[negind], S[s])
                    S[s]=cas10_1; del cas10_1

            elif s == '110':
                if negind == -1:
                    cas110 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==negind), Negatives[negind], cas)
                    S[s]=cas110; del cas110
                else:
                    cas110 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==negind), Negatives[negind], S[s])
                    S[s]=cas110; del cas110

            elif s == '1_10':
                if negind == -1:
                    cas1_10 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==negind), Negatives[negind], cas)
                    S[s]=cas1_10; del cas1_10
                else:
                    cas1_10 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==negind), Negatives[negind], S[s])
                    S[s]=cas1_10; del cas1_10

            elif s == '_101':
                if negind == -1:
                    cas_101 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,0), (1,0)), 'constant')[1:,:,:-1]==negind), Negatives[negind], cas)
                    S[s]=cas_101; del cas_101
                else:
                    cas_101 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,0), (1,0)), 'constant')[1:,:,:-1]==negind), Negatives[negind], S[s])
                    S[s]=cas_101; del cas_101

            elif s == '_10_1':
                if negind == -1:
                    cas_10_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,0), (0,1)), 'constant')[1:,:,1:]==negind), Negatives[negind], cas)
                    S[s]=cas_10_1; del cas_10_1
                else:
                    cas_10_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,0), (0,1)), 'constant')[1:,:,1:]==negind), Negatives[negind], S[s])
                    S[s]=cas_10_1; del cas_10_1

            elif s == '_110':
                if negind == -1:
                    cas_110 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(1,0), (0,0)), 'constant')[1:,:-1,:]==negind), Negatives[negind], cas)
                    S[s]=cas_110; del cas_110
                else:
                    cas_110 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(1,0), (0,0)), 'constant')[1:,:-1,:]==negind), Negatives[negind], S[s])
                    S[s]=cas_110; del cas_110

            elif s == '_1_10':
                if negind == -1:
                    cas_1_10 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,1), (0,0)), 'constant')[1:,1:,:]==negind), Negatives[negind], cas)
                    S[s]=cas_1_10; del cas_1_10
                else:
                    cas_1_10 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,1), (0,0)), 'constant')[1:,1:,:]==negind), Negatives[negind], S[s])
                    S[s]=cas_1_10; del cas_1_10
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

                    ''' ----------------------------------------- 3rd shell -------------------> 8 nearest neighbours '''
            elif s == '111':
                if negind == -1:
                    cas111 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==negind), Negatives[negind], cas)
                    S[s]=cas111; del cas111
                else:
                    cas111 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==negind), Negatives[negind], S[s])
                    S[s]=cas111; del cas111

            elif s == '11_1':
                if negind == -1:
                    cas11_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==negind), Negatives[negind], cas)
                    S[s]=cas11_1; del cas11_1
                else:
                    cas11_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==negind), Negatives[negind], S[s])
                    S[s]=cas11_1; del cas11_1

            elif s == '1_11':
                if negind == -1:
                    cas1_11 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==negind), Negatives[negind], cas)
                    S[s]=cas1_11; del cas1_11
                else:
                    cas1_11 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==negind), Negatives[negind], S[s])
                    S[s]=cas1_11; del cas1_11

            elif s == '1_1_1':
                if negind == -1:
                    cas1_1_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==negind), Negatives[negind], cas)
                    S[s]=cas1_1_1; del cas1_1_1
                else:
                    cas1_1_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==negind), Negatives[negind], S[s])
                    S[s]=cas1_1_1; del cas1_1_1

            elif s == '_111':
                if negind == -1:
                    cas_111 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==negind), Negatives[negind], cas)
                    S[s]=cas_111; del cas_111
                else:
                    cas_111 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==negind), Negatives[negind], S[s])
                    S[s]=cas_111; del cas_111

            elif s == '_11_1':
                if negind == -1:
                    cas_11_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==negind), Negatives[negind], cas)
                    S[s]=cas_11_1; del cas_11_1
                else:
                    cas_11_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==negind), Negatives[negind], S[s])
                    S[s]=cas_11_1; del cas_11_1

            elif s == '_1_11':
                if negind == -1:
                    cas_1_11 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==negind), Negatives[negind], cas)
                    S[s]=cas_1_11; del cas_1_11
                else:
                    cas_1_11 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==negind), Negatives[negind], S[s])
                    S[s]=cas_1_11; del cas_1_11

            elif s == '_1_1_1':
                if negind == -1:
                    cas_1_1_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==negind), Negatives[negind], cas)
                    S[s]=cas_1_1_1; del cas_1_1_1
                else:
                    cas_1_1_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==negind), Negatives[negind], S[s])
                    S[s]=cas_1_1_1; del cas_1_1_1
                ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

                ''' ----------------------------------------- extra neighbours -------------------> gettin' rid of mesh dependency :)'''
            elif s == '012':
               if negind == -1:
                    cas012 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==negind), Negatives[negind], cas)
                    S[s]=cas012; del cas012
               else:
                    cas012 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==negind), Negatives[negind], S[s])
                    S[s]=cas012; del cas012

            elif s == '01_2':
               if negind == -1:
                    cas01_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==negind), Negatives[negind], cas)
                    S[s]=cas01_2; del cas01_2
               else:
                    cas01_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==negind), Negatives[negind], S[s])
                    S[s]=cas01_2; del cas01_2

            elif s == '0_12':
               if negind == -1:
                    cas0_12 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==negind), Negatives[negind], cas)
                    S[s]=cas0_12; del cas0_12
               else:
                    cas0_12 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==negind), Negatives[negind], S[s])
                    S[s]=cas0_12; del cas0_12

            elif s == '0_1_2':
               if negind == -1:
                    cas0_1_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==negind), Negatives[negind], cas)
                    S[s]=cas0_1_2; del cas0_1_2
               else:
                    cas0_1_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==negind), Negatives[negind], S[s])
                    S[s]=cas0_1_2; del cas0_1_2

            elif s == '021':
               if negind == -1:
                    cas021 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==negind), Negatives[negind], cas)
                    S[s]=cas021; del cas021
               else:
                    cas021 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==negind), Negatives[negind], S[s])
                    S[s]=cas021; del cas021

            elif s == '02_1':
               if negind == -1:
                    cas02_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==negind), Negatives[negind], cas)
                    S[s]=cas02_1; del cas02_1
               else:
                    cas02_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==negind), Negatives[negind], S[s])
                    S[s]=cas02_1; del cas02_1

            elif s == '0_21':
               if negind == -1:
                    cas0_21 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==negind), Negatives[negind], cas)
                    S[s]=cas0_21; del cas0_21
               else:
                    cas0_21 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==negind), Negatives[negind], S[s])
                    S[s]=cas0_21; del cas0_21

            elif s == '0_2_1':
               if negind == -1:
                    cas0_2_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==negind), Negatives[negind], cas)
                    S[s]=cas0_2_1; del cas0_2_1
               else:
                    cas0_2_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==negind), Negatives[negind], S[s])
                    S[s]=cas0_2_1; del cas0_2_1

               '''~~~~~~~~~~~~~~~~~~~~~'''
               '''........................................ another set of extra neighbours .................. '''

            elif s == '002':
                if negind == -1:
                    cas002 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==negind), Negatives[negind], cas)
                    S[s]=cas002; del cas002
                else:
                    cas002 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==negind), Negatives[negind], S[s])
                    S[s]=cas002; del cas002

            elif s == '00_2':
                if negind == -1:
                    cas00_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==negind), Negatives[negind], cas)
                    S[s]=cas00_2; del cas00_2
                else:
                    cas00_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==negind), Negatives[negind], S[s])
                    S[s]=cas00_2; del cas00_2

            elif s == '020':
                if negind == -1:
                    cas020 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==negind), Negatives[negind], cas)
                    S[s]=cas020; del cas020
                else:
                    cas020 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==negind), Negatives[negind], S[s])
                    S[s]=cas020; del cas020

            elif s == '0_20':
                if negind == -1:
                    cas0_20 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==negind), Negatives[negind], cas)
                    S[s]=cas0_20; del cas0_20
                else:
                    cas0_20 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==negind), Negatives[negind], S[s])
                    S[s]=cas0_20; del cas0_20

            elif s == '022':
               if negind == -1:
                    cas022 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==negind), Negatives[negind], cas)
                    S[s]=cas022; del cas022
               else:
                    cas022 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==negind), Negatives[negind], S[s])
                    S[s]=cas022; del cas022

            elif s == '02_2':
               if negind == -1:
                    cas02_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==negind), Negatives[negind], cas)
                    S[s]=cas02_2; del cas02_2
               else:
                    cas02_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==negind), Negatives[negind], S[s])
                    S[s]=cas02_2; del cas02_2

            elif s == '0_22':
               if negind == -1:
                    cas0_22 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==negind), Negatives[negind], cas)
                    S[s]=cas0_22; del cas0_22
               else:
                    cas0_22 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==negind), Negatives[negind], S[s])
                    S[s]=cas0_22; del cas0_22

            elif s == '0_2_2':
               if negind == -1:
                    cas0_2_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==negind), Negatives[negind], cas)
                    S[s]=cas0_2_2; del cas0_2_2
               else:
                    cas0_2_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==negind), Negatives[negind], S[s])
                    S[s]=cas0_2_2; del cas0_2_2

          
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 
    




    # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   3.   P  H  A  S   E   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    F = faza.copy()
    for g in grain_ID:
           grain_ID = random_selection(grain_ID)
           grain = item
           directions = smeri.copy()
           seznam_premikov = []
           for d in directions:
               directions = random_selection(directions)
               s = item
               
               ''' ----------------------------------------- 1st shell -------------------> 6 nearest neighbours '''
               if s == '001':
                   dij = np.array([0,0,1])
                   wij = fn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza001 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain), grain, faza)
                   seznam_premikov.append(faza001); del faza001

               elif s == '00_1':
                   dij = np.array([0,0,-1])
                   wij = fn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza00_1 = np.where((faza==0)& (taula==0)& (likvid==0)& (lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain), grain, faza)
                   seznam_premikov.append(faza00_1); del faza00_1

               elif s == '010':
                   dij = np.array([0,1,0])
                   wij = fn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza010 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain), grain, faza)
                   seznam_premikov.append(faza010); del faza010

               elif s == '0_10':
                   dij = np.array([0,-1,0])
                   wij = fn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza0_10 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain), grain, faza)
                   seznam_premikov.append(faza0_10); del faza0_10

               elif s == '100':
                   dij = np.array([1,0,0])
                   wij = fn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza100 = np.where((faza==0)& (taula==0)& (likvid==0)& (lij>=Dsr_1st(R) )&(np.pad(faza,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain), grain, faza)
                   seznam_premikov.append(faza100); del faza100

               elif s == '_100':
                   dij = np.array([-1,0,0])
                   wij = fn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza_100 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_1st(R) )&(np.pad(faza,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain), grain, faza)
                   seznam_premikov.append(faza_100); del faza_100
                   ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

                   ''' ----------------------------------------- 2nd shell -------------------> 12 nearest neighbours '''
               elif s == '011':
                   dij = np.array([0,1,1])
                   wij = W(asc[grain]['oi'], dij)
                   lij = sn*wij*vg*S[s]
                   faza011 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R) )&(np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain), grain, faza)
                   seznam_premikov.append(faza011); del faza011

               elif s == '01_1':
                   dij = np.array([0,1,-1])
                   wij = sn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza01_1 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R) )&(np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain), grain, faza)
                   seznam_premikov.append(faza01_1); del faza01_1

               elif s == '0_11':
                   dij = np.array([0,-1,1])
                   wij = sn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza0_11 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R) )&(np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain), grain, faza)
                   seznam_premikov.append(faza0_11); del faza0_11

               elif s == '0_1_1':
                   dij = np.array([0,-1,-1])
                   wij = sn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza0_1_1 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R) )&(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain), grain, faza)
                   seznam_premikov.append(faza0_1_1); del faza0_1_1

               elif s == '101':
                   dij = np.array([1,0,1])
                   wij = sn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza101 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==grain), grain, faza)
                   seznam_premikov.append(faza101); del faza101

               elif s == '10_1':
                   dij = np.array([1,0,-1])
                   wij = sn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza10_1 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==grain), grain, faza)
                   seznam_premikov.append(faza10_1); del faza10_1

               elif s == '110':
                   dij = np.array([1,1,0])
                   wij = sn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza110 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==grain), grain, faza)
                   seznam_premikov.append(faza110); del faza110

               elif s == '1_10':
                   dij = np.array([1,-1,0])
                   wij = sn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza1_10 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain), grain, faza)
                   seznam_premikov.append(faza1_10); del faza1_10

               elif s == '_101':
                   dij = np.array([-1,0,1])
                   wij = sn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza_101 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((0,1),(0,0), (1,0)), 'constant')[1:,:,:-1]==grain), grain, faza)
                   seznam_premikov.append(faza_101); del faza_101

               elif s == '_10_1':
                   dij = np.array([-1,0,-1])
                   wij = sn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza_10_1 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((0,1),(0,0), (0,1)), 'constant')[1:,:,1:]==grain), grain, faza)
                   seznam_premikov.append(faza_10_1); del faza_10_1

               elif s == '_110':
                   dij = np.array([-1,1,0])
                   wij = sn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza_110 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((0,1),(1,0), (0,0)), 'constant')[1:,:-1,:]==grain), grain, faza)
                   seznam_premikov.append(faza_110); del faza_110

               elif s == '_1_10':
                   dij = np.array([-1,-1,0])
                   wij = sn*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza_1_10 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((0,1),(0,1), (0,0)), 'constant')[1:,1:,:]==grain), grain, faza)
                   seznam_premikov.append(faza_1_10); del faza_1_10
                   ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''          

                   ''' ----------------------------------------- 3rd shell -------------------> 8 nearest neighbours '''
               elif s == '111':
                   dij = np.array([1,1,1])
                   wij = W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza111 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain), grain, faza)
                   seznam_premikov.append(faza111); del faza111

               elif s == '11_1':
                   dij = np.array([1,1,-1])
                   wij = W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza11_1 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain), grain, faza)
                   seznam_premikov.append(faza11_1); del faza11_1

               elif s == '1_11':
                   dij = np.array([1,-1,1])
                   wij = W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza1_11 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain), grain, faza)
                   seznam_premikov.append(faza1_11); del faza1_11

               elif s == '1_1_1':
                   dij = np.array([1,-1,-1])
                   wij = W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza1_1_1 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain), grain, faza)
                   seznam_premikov.append(faza1_1_1); del faza1_1_1

               elif s == '_111':
                   dij = np.array([-1,1,1])
                   wij = W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza_111 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==grain), grain, faza)
                   seznam_premikov.append(faza_111); del faza_111

               elif s == '_11_1':
                   dij = np.array([-1,1,-1])
                   wij = W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza_11_1 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==grain), grain, faza)
                   seznam_premikov.append(faza_11_1); del faza_11_1

               elif s == '_1_11':
                   dij = np.array([-1,-1,1])
                   wij = W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza_1_11 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==grain), grain, faza)
                   seznam_premikov.append(faza_1_11); del faza_1_11

               elif s == '_1_1_1':
                   dij = np.array([-1,-1,-1])
                   wij = W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza_1_1_1 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==grain), grain, faza)
                   seznam_premikov.append(faza_1_1_1); del faza_1_1_1
                   ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

                   ''' ----------------------------------------- extra neighbours -------------------> gettin' rid of mesh dependency :)'''
               elif s == '012':
                   dij = np.array([0,1,2])
                   wij = en*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza012 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )
                                                   &(np.pad(faza,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)
                                                   & ((np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)                     # 001
                                                          | (np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)                # 011
                                                          |  (np.pad(faza,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)                  # 002
                                                          |  (np.pad(faza,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)               # 021
                                                          |  (np.pad(faza,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==grain))            # 022
                                                    , grain, faza) 
                   seznam_premikov.append(faza012); del faza012

               elif s == '01_2':
                   dij = np.array([0,1,-2])
                   wij = en*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza01_2 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )
                                                      &(np.pad(faza,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)
                                                      & ((np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)                     # 00_1
                                                             | (np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)                # 01_1
                                                             |  (np.pad(faza,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)                  # 00_2
                                                             | (np.pad(faza,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)                # 02_1
                                                             |  (np.pad(faza,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==grain))            # 02_2
                                                     , grain, faza) 
                   seznam_premikov.append(faza01_2); del faza01_2

               elif s == '0_12':
                   dij = np.array([0,-1,2])
                   wij = en*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza0_12 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )
                                                      &(np.pad(faza,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)
                                                      & ((np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)                    # 001
                                                             | (np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)                # 0_11
                                                             |  (np.pad(faza,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)                 # 002
                                                             | (np.pad(faza,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)                # 0_21
                                                             |  (np.pad(faza,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==grain))            # 0_22
                                                      , grain, faza)
                   seznam_premikov.append(faza0_12); del faza0_12

               elif s == '0_1_2':
                   dij = np.array([0,-1,-2])
                   wij = en*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza0_1_2 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )
                                                         &(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)
                                                         & ((np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)                  # 00_1
                                                                | (np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)              # 0_1_1
                                                                |  (np.pad(faza,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)               # 00_2
                                                                | (np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)              # 0_2_1
                                                                |  (np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain))          # 0_2_2
                                                         , grain, faza) 
                   seznam_premikov.append(faza0_1_2); del faza0_1_2

               elif s == '021':
                   dij = np.array([0,2,1])
                   wij = en*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza021 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )
                                                   &(np.pad(faza,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)
                                                   & ((np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)                      # 010
                                                          | (np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)                # 011
                                                          |  (np.pad(faza,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)                  # 020
                                                          | (np.pad(faza,((0,0),(3,0), (0,0)), 'constant')[:,:-3,:]==grain)                   # 030
                                                          |  (np.pad(faza,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain))            # 012
                                                   , grain, faza)
                   seznam_premikov.append(faza021); del faza021

               elif s == '02_1':
                   dij = np.array([0,2,-1])
                   wij = en*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza02_1 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )
                                                      &(np.pad(faza,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)
                                                      & ((np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)                  # 010
                                                             | (np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)              # 01_1
                                                             |  (np.pad(faza,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)               # 020
                                                             | (np.pad(faza,((0,0),(3,0), (0,0)), 'constant')[:,:-3,:]==grain)                # 030
                                                             |  (np.pad(faza,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain))          # 01_2
                                                      , grain, faza) 
                   seznam_premikov.append(faza02_1); del faza02_1

               elif s == '0_21':
                   dij = np.array([0,-2,1])
                   wij = en*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza0_21 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )
                                                      &(np.pad(faza,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)
                                                      & ((np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)                   # 0_10
                                                             | (np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)             # 0_11
                                                             |  (np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)               # 0_20
                                                             | (np.pad(faza,((0,0),(0,3), (0,0)), 'constant')[:,3:,:]==grain)                # 0_30
                                                             |  (np.pad(faza,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain))         # 0_12
                                                       , grain, faza)
                   seznam_premikov.append(faza0_21); del faza0_21

               elif s == '0_2_1':
                   dij = np.array([0,-2,-1])
                   wij = en*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza0_2_1 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )
                                                         &(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)
                                                         & ((np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)               # 0_10
                                                                | (np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)          # 0_1_1
                                                                |  (np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)           # 0_20
                                                                | (np.pad(faza,((0,0),(0,3), (0,0)), 'constant')[:,3:,:]==grain)            # 0_30
                                                               |  (np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain))       # 0_1_2
                                                         , grain, faza) 
                   seznam_premikov.append(faza0_2_1); del faza0_2_1

                   '''~~~~~~~~~~~~~~~~~~~~~'''
                   '''........................................ another set of extra neighbours   ::: Dsr_5th(R) ::: .................. '''

               elif s == '002':
                   dij = np.array([0,0,2])
                   wij = en2*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza002 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )
                                                   &(np.pad(faza,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)
                                                   & ((np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)                  # 001
                                                          | (np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)             # 011
                                                          |  (np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)             # 0_11
                                                          | (np.pad(faza,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)             # 012           
                                                          |  (np.pad(faza,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain))          # 0_12
                                                , grain, faza) 
                   seznam_premikov.append(faza002); del faza002

               elif s == '00_2':
                   dij = np.array([0,0,-2])
                   wij = en2*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza00_2 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )
                                                      &(np.pad(faza,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)
                                                      & ((np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)                # 00_1
                                                             | (np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)          # 01_1
                                                             |  (np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)          # 0_1_1
                                                             | (np.pad(faza,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)          # 01_2
                                                             |  (np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain))       # 0_1_2
                                                     , grain, faza)
                   seznam_premikov.append(faza00_2); del faza00_2

               elif s == '020':
                   dij = np.array([0,2,0])
                   wij = en2*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza020 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )
                                                   &(np.pad(faza,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)
                                                   & ((np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)                  # 010 
                                                          | (np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)             # 011
                                                          |  (np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)             # 01_1
                                                          | (np.pad(faza,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)             # 021
                                                          |  (np.pad(faza,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain))          # 02_1  
                                                , grain, faza)
                   seznam_premikov.append(faza020); del faza020

               elif s == '0_20':
                   dij = np.array([0,-2,0])
                   wij = en2*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza0_20 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )
                                                      &(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)
                                                      & ((np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)               # 0_10
                                                             | (np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)         # 0_11
                                                             |  (np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)         # 0_1_1
                                                             | (np.pad(faza,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)         # 0_21
                                                             |  (np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain))      # 0_2_1
                                                      , grain, faza)     
                   seznam_premikov.append(faza0_20); del faza0_20

                   '''~~~~~~~~~~~~~~~~~~~~~'''
                   '''........................................ another set of extra neighbours   ::: Dsr_6th(R) ::: .................. '''

               elif s == '022':
                   dij = np.array([0,2,2])
                   wij = en3*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza022 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )
                                                   &(np.pad(faza,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==grain)
                                                   & ((np.pad(faza,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)             # 012
                                                          | (np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)          # 011
                                                          |  (np.pad(faza,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)         # 021
                                                          | (np.pad(faza,((0,0),(1,0), (3,0)), 'constant')[:,:-1,:-3]==grain)          # 013
                                                          |  (np.pad(faza,((0,0),(3,0), (1,0)), 'constant')[:,:-3,:-1]==grain))      # 031
                                                   , grain, faza) 
                   seznam_premikov.append(faza022); del faza022

               elif s == '02_2':
                   dij = np.array([0,2,-2])
                   wij = en3*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza02_2 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )
                                                      &(np.pad(faza,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==grain)
                                                      & ((np.pad(faza,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)           # 01_2
                                                             | (np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)        # 01_1
                                                             |  (np.pad(faza,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)       # 02_1
                                                             | (np.pad(faza,((0,0),(1,0), (0,3)), 'constant')[:,:-1,3:]==grain)        # 01_3
                                                             |  (np.pad(faza,((0,0),(3,0), (0,1)), 'constant')[:,:-3,1:]==grain))    # 03_1
                                                       , grain, faza) 
                   seznam_premikov.append(faza02_2); del faza02_2

               elif s == '0_22':
                   dij = np.array([0,-2,2])
                   wij = en3*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza0_22 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )
                                                      &(np.pad(faza,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==grain)
                                                      & ((np.pad(faza,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)           # 0_12
                                                             | (np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)        # 0_11
                                                             |  (np.pad(faza,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)       # 0_21
                                                             | (np.pad(faza,((0,0),(0,1), (3,0)), 'constant')[:,1:,:-3]==grain)        # 0_13
                                                             |  (np.pad(faza,((0,0),(0,3), (1,0)), 'constant')[:,3:,:-1]==grain))    # 0_31 
                                                      , grain, faza) 
                   seznam_premikov.append(faza0_22); del faza0_22

               elif s == '0_2_2':
                   dij = np.array([0,-2,-2])
                   wij = en3*W(asc[grain]['oi'], dij)
                   lij = wij*vg*S[s]
                   faza0_2_2 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )
                                                         &(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)
                                                         & ((np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)         # 0_1_2
                                                                | (np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)      # 0_1_1
                                                                |  (np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)     # 0_2_1
                                                                | (np.pad(faza,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)      # 0_1_3
                                                                |  (np.pad(faza,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain))  # 0_3_1
                                                         , grain, faza)
                   seznam_premikov.append(faza0_2_2); del faza0_2_2



           total=merging(seznam_premikov[0], seznam_premikov[1])
           #total=seznam_premikov[0]
           times = np.array(list(S.values()))
           cas_total = merging(times[0], times[1])
           #cas_total = times[0]
           if len(smeri)> 2:
               for move in range(2, len(smeri)):
                   total=merging(total, seznam_premikov[move])
                   cas_total = merging(cas_total, times[move])
              
           faza = total.copy()
           cas = cas_total.copy()

           grain_size = np.count_nonzero(faza == grain)
           asc[grain]['grain size'] = grain_size

           # ----------------------------------------- this function selects inactive grains (i.e. grains, which no longer grow - either domain limits were reached and/or has reached another grain)
           if delete_inactive_grains_ID:
              active=np.isin(faza-F, grain)
              
              if np.all(active==False):
                 inactive_grains.append(grain)
                 IG[grain].append(i)
              else:
                 AG[grain].append(i)
           
           # -----------------------------------------
        
    grain_ID = np.array(list(asc.keys()))
    if delete_inactive_grains_ID:
       grain_ID = np.array(list(set(grain_ID)-set(inactive_grains)))    # erasing inactive grains from grain_ID register
    
    if not np.all(faza==F):
        negind -=1
        cas_total[np.isin(faza-F, grain_ID)]= negind
        Negatives[negind]=0

    try:
        cas  = cas_total.copy()
    except NameError:
        pass

    rgb_snap = np.zeros((faza.shape[0], faza.shape[1], faza.shape[2], 3)).astype('int')
    for zrno in asc:
        rgb_snap[faza==zrno]=asc[zrno]['rgb']
    try:
        if save_flashy_as_RGB:
           np.save(PATH+mapa+'/'+flashies+'/flashy_snap_'+str(i)+'.npy', rgb_snap)
        else:
           np.save(PATH+mapa+'/'+flashies+'/flashy_snap_'+str(i)+'.npy', faza)
    except FileNotFoundError:
        os.mkdir(PATH+mapa+'/'+flashies)
        if save_flashy_as_RGB:
           np.save(PATH+mapa+'/'+flashies+'/flashy_snap_'+str(i)+'.npy', rgb_snap)
        else:
           np.save(PATH+mapa+'/'+flashies+'/flashy_snap_'+str(i)+'.npy', faza)
        

    step_cpu = round(time.time()- step_time_start, 1)
    #print('Časovni korak: ',i,'............ čas računanja: ',step_cpu,' sekund........... št. aktivnih zrn: ',str(grain_ID.shape[0]),'.........št. vseh zrn: ',grain_counter)

    if not bool(grain_ID.shape[0]) and grain_counter:
       break
    

process_time = time.time() - start_time
print(150*'*')
print(50*' ','DONE!!  Total computing time =  ',round(process_time, 3),'  seconds.')
print(grain_counter,' grains were formed!'); print(150*'*'); print()

with open(PATH+mapa+'/nuclei_data.json', 'w') as nuks:              # Writing data of ALL nuclei as .json file, but values must be list NOT np.array !!!
   asc_list =asc.copy()
   for nuk in asc:
      asc_list[nuk]['oi']=asc[nuk]['oi'].tolist()
      asc_list[nuk]['rgb']=asc[nuk]['rgb'].tolist()
   json.dump(asc_list, nuks)


# ......................................... PRIKAZ in IZPIS REZULTATOV ..................................................

#plt.imshow(rgb_snap[0]); plt.figure()
plt.imshow(faza[0])

def Figure_Save():
   plt.gca().set_axis_off()
   plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
   plt.margins(0,0)
   plt.gca().xaxis.set_major_locator(plt.NullLocator())
   plt.gca().yaxis.set_major_locator(plt.NullLocator())

   #delo_doma_path = 'C:/Users/akocjan/Desktop/delo od doma II. val/komunikacija s Tijanom/ISO_MD/cell-6/multi grains/'
   delo_doma_path = 'C:/Users/akocjan/Desktop/delo od doma II. val/komunikacija s Tijanom/real_SLM_MD/cell-3/'

   #ime_slike = 'ISO_MD_alfa_'+str(alfa)+'_beta_'+str(beta)+'_gama_'+str(gama)  +',  fn= '+str(fn)+', sn= '+str(sn)+', en= '+str(en)+', size= '+str(asc[1]['grain size'])+' px'
   ime_slike = 'SLM_MD alfa='+str(round(alfa, 3))+', beta='+str(round(beta, 3))+', gama='+str(round(gama, 3))  +', eps= '+str(epsilon)+', size= '+str(asc[1]['grain size'])+' px'

   plt.savefig(delo_doma_path+ime_slike+'.png', bbox_inches='tight', pad_inches = 0)

#Figure_Save()

''' ---------------------------------------------------------------------------------- the end -------------------------------------------------------------------------------------------------------------- '''























# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Growth (time FOR-LOOP, grains FOR-LOOP, single direction)..
"""
for i in range(0,4):
    for negind in Negatives:
        time_counter(negind)
        if negind == -1:
            cas001 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==negind), Negatives[negind], cas)
        else:
            cas001 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas001,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==negind), Negatives[negind], cas001)

    '''............................................. tole naj bo funkcija rasti kot npr. def Growth_FAST_code(count, direction):  ..........................................................'''
    lij = wij*vg*cas001
    grc=0
    for grain in grain_ID:
        if grc ==0:                  
            faza001 = np.where((faza==0)&(lij>=cell)&(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain), grain, faza)
        else:
            faza001 = np.where((faza==0)&(lij>=cell)&(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain), grain, faza001)
        grc+=1
    '''----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'''
    
    if not np.all(faza001==faza):
        negind -=1
        cas001[np.isin(faza001-faza, grain_ID)]= negind
        Negatives[negind]=0
    
    faza = faza001.copy()
    cas  = cas001.copy()
"""

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#__________________________________________Growth (step by step, single grain, single direction)..
"""
#negind = -1
N = list(Negatives.keys())
time_counter(negind)

''' 1st time step ................................................................................................................................................................................................................. '''
zeit_1=0; zeit_1+=dt
cas001 = np.where((faza != grain_ID) & (np.pad(cas,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==-1), zeit_1, cas)
lij = wij*vg*cas001
faza001 = np.where((faza==0)&(lij>=cell)&(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain_ID), grain_ID, faza) 
if not np.all(faza001==faza):
    negind -=1
    cas001[faza001-faza == grain_ID]= negind
#------------------------------------------------>>>
faza = faza001.copy()
cas  = cas001.copy()
print(50*'.','STEP #1 .... zeit_1 = ',zeit_1)
print('cas'); print(cas)
print('faza'); print(faza)  
#------------------------------------------------>>> ..........................................................................................................................................................................


''' 2nd time step .................................................................................................................................................................................................................'''
zeit_1+=dt
zeit_2=0; zeit_2+=dt
cas001 = np.where((faza != grain_ID) & (np.pad(cas,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==-1), zeit_1, cas)
cas001 = np.where((faza != grain_ID) & (np.pad(cas,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==-2), zeit_2, cas001)

#cas001_1 = np.where((faza != grain_ID) & (np.pad(cas,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==-1), zeit_1, cas)
#cas001_2 = np.where((faza != grain_ID) & (np.pad(cas,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==-2), zeit_2, cas)
#cas001=merging(cas001_2, cas001_1)

#lij_1 = wij*vg*cas001_1
#lij_2 = wij*vg*cas001_2

lij = wij*vg*cas001

#faza001_1 = np.where((faza==0)&(lij_1>=cell)&(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain_ID), grain_ID, faza) 
#faza001_2 = np.where((faza==0)&(lij_2>=cell)&(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain_ID), grain_ID, faza) 

faza001 = np.where((faza==0)&(lij>=cell)&(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain_ID), grain_ID, faza) 

#faza001 = merging(faza001_1, faza001_2)

if not np.all(faza001==faza):
    negind -=1
    cas001[faza001-faza == grain_ID]= negind

#cas001[faza001-faza == grain_ID]= -2

#------------------------------------------------>>>
faza = faza001.copy()
cas  = cas001.copy()
print(50*'.','STEP #2 .... zeit_1 = ',zeit_1,',  zeit_2 = ',zeit_2)
print('cas'); print(cas)
print('faza'); print(faza)
#------------------------------------------------>>> ........................................................................................................................................................................


''' 3rd time step .................................................................................................................................................................................................................'''
zeit_1+=dt
zeit_2+=dt
zeit_3=0; zeit_3+=dt

cas001 = np.where((faza != grain_ID) & (np.pad(cas,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==-1), zeit_1, cas)
cas001 = np.where((faza != grain_ID) & (np.pad(cas,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==-2), zeit_2, cas001)
cas001 = np.where((faza != grain_ID) & (np.pad(cas,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==-3), zeit_3, cas001)

lij = wij*vg*cas001

faza001 = np.where((faza==0)&(lij>=cell)&(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain_ID), grain_ID, faza) 

if not np.all(faza001==faza):
    negind -=1
    cas001[faza001-faza == grain_ID]= negind

#cas001[faza001-faza == grain_ID]= -3

#------------------------------------------------>>>
faza = faza001.copy()
cas  = cas001.copy()
 #------------------------------------------------>>> ........................................................................................................................................................................

# vpelji slovar Negatives v step-by-step
"""











