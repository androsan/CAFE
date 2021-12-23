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
np.random.seed(90983)

""" =====================================  F  U  N  C  T  I  O  N  S ================================================================"""
def f(n):
    f=np.load('C:/sm-2018-w64-0-3/WORK/SLM_2D_Source_Files/post processing database/0002/'+
              'INTER  time=1, space=8  Z[0-9], X[15-27], Y[12-97], 1500°C, N=12/2D 1st order Moore, real field/flashies_RGB/flashy_RGB_{}.npy'.format(n))
    return f

# MPI executor functions
import MPI_executor_functions as mef

# Loading temperature field (unit: KELVIN)in corresponding space (YP) and time (TR) related folder
def Load_NPY(yp, tr, x):
    npy = np.load(PATH+mapa+yp+'/'+tr+time_factor_folder+'/salome_'+str(x)+'.npy')
    return npy

def Stack_2 (yp1st, yp2nd, tm, cif):
    yp_1st =  Load_NPY(yp1st, tm, cif)
    yp_2nd = Load_NPY(yp2nd, tm, cif)
    yp_stacked = np.dstack((yp_1st,yp_2nd ))
    return yp_stacked

def Make_Y_Partitions(Ymin, Ymax, slices):
    dDY = (Ymax-Ymin)/slices
    a=[Ymin+i*dDY for i in range(slices)]
    b=[Ymin+i*dDY+1 for i in range(1,slices+1)]
    YP={}
    for num,(i,j) in enumerate(zip(a,b)):
        YP['YP'+str(num)] = [int(i),int(j)]
    return YP

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
def nukleacija_povrsina(pp):                                                   # heterogeneous nucleation (at the liquid-solid interface ---> melt-pool border)
    ''' parameters of heterogeneous nucleation'''
    ns_max =                      5e10                                                 #  Number of new nuclei at optimal dT; unit: per SQUARE METER  [m^-2]
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
   taula[u<(Tm)] = -1
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
    xran1 = np.random.randint(0, posibilities+1)
    xran2 = np.random.randint(0, posibilities+1)
    xran3 = np.random.randint(0, posibilities+1)
    xran = (xran1+xran2+xran3)/3  
    x = xran / posibilities
    y = np.random.randint(0, xran+1)/posibilities
    return x,y

def merging(prva,druga):
    tot = np.where(druga==0, prva, druga)
    return tot

def time_counter(neg_index):
        global Negatives
        Negatives[neg_index] += dt
        
def random_selection(x, erase):
    global item
    if len(x)!=0:
        item = random.choice(x)
        if erase:
            x=np.delete(x,np.where(x==item))
    else:
        pass
    return x

def Ordered_Grains_n_Directions(g_ID, smeri_arrays):
    selekcija=[]
    for gid in g_ID:
        for sa in smeri_arrays:
            for sm in sa:
                selekcija.append((gid, sm))
    return selekcija

def Random_Grains_Ordered_Directions(g_ID, smer):
    selekcija=[]
    while True:
        if any(g_ID):
            g_ID = random_selection(g_ID, True)
            grain = item
            for j in smer:
                selekcija.append((grain,j))
        else:
            break     
    return selekcija

def Random_Grains_n_Directions(g_ID, smer):
    GD = {}
    for gid in g_ID:
       GD[gid] = smer.copy()
    selekcija=[]
    while True:
        if any(g_ID):
            g_ID = random_selection(g_ID, False)
            grain = item
            counter = 0
        else:
            break
        while counter==0:
            if any(GD[grain]):
                GD[grain] = random_selection(GD[grain], True)
                s = item
                selekcija.append((grain,s))
                counter+=1
                break
            else:
                g_ID=g_ID[g_ID!=grain]
                break
        continue
    return selekcija

def Random_Grains_Directions_Segmented(g_ID, smeri_arrays):
    SM=[]
    for sm in smeri_arrays:
        SM+=Random_Grains_n_Directions(g_ID, sm)
    return SM

def TT(x):
   temp = np.load(PATH+mapa+'salome_'+str(x)+'.npy')[z_min:z_max, x_min:x_max, y_min:y_max]
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

def Save_KickOff():
    kickoff_folder = 'kickoff_data/'
    if not os.path.isdir(PATH+mapa+track+kickoff_folder):
       os.mkdir(PATH+mapa+track+kickoff_folder)

    np.save(PATH+mapa+track+kickoff_folder+'faza_kickoff.npy', faza)
    np.save(PATH+mapa+track+kickoff_folder+'cas_kickoff.npy', cas)
    np.save(PATH+mapa+track+kickoff_folder+'rgb_snap_kickoff.npy', rgb_snap)
    try:
        np.save(PATH+mapa+cuts_RGB+'cut_RGB_'+str(cut_count)+'.npy', rgb_snap[:,:,:cutoff_limit,:])
        np.save(PATH+mapa+cuts_RGB+'cut_RGB_'+str(cut_count+1)+'.npy', rgb_snap[:,:,cutoff_limit:,:])
        np.save(PATH+mapa+cuts_faza+'cut_faza_'+str(cut_count)+'.npy', faza[:,:,:cutoff_limit])
        np.save(PATH+mapa+cuts_faza+'cut_faza_'+str(cut_count+1)+'.npy', faza[:,:,cutoff_limit:])
    except FileNotFoundError:
        os.mkdir(PATH+mapa+cuts_RGB)
        os.mkdir(PATH+mapa+cuts_faza)
        np.save(PATH+mapa+cuts_RGB+'cut_RGB_'+str(cut_count)+'.npy', rgb_snap[:,:,:cutoff_limit,:])
        np.save(PATH+mapa+cuts_RGB+'cut_RGB_'+str(cut_count+1)+'.npy', rgb_snap[:,:,cutoff_limit:,:])
        np.save(PATH+mapa+cuts_faza+'cut_faza_'+str(cut_count)+'.npy', faza[:,:,:cutoff_limit])
        np.save(PATH+mapa+cuts_faza+'cut_faza_'+str(cut_count+1)+'.npy', faza[:,:,cutoff_limit:])
    with open(PATH+mapa+track+kickoff_folder+'nuclei_kickoff.json', 'w') as nuks:              # Writing data of ALL nuclei as .json file, but values must be list NOT np.array !!!
        asc_list = asc.copy()
        for nuk in asc:
            try:
                asc_list[nuk]['oi']=asc[nuk]['oi'].tolist()
                asc_list[nuk]['rgb']=asc[nuk]['rgb'].tolist()
            except AttributeError:
                asc_list[nuk]['oi']=asc[nuk]['oi']
                asc_list[nuk]['rgb']=asc[nuk]['rgb']
        json.dump(asc_list, nuks)
    with open(PATH+mapa+track+kickoff_folder+'negatives_kickoff.json', 'w') as negs:
        json.dump(Negatives, negs)
    with open(PATH+mapa+track+kickoff_folder+'S_kickoff.json', 'w') as S_kick:
        S_list = S.copy()
        for _s_ in S:
            try:
                S_list[_s_]=S[_s_].tolist()
            except AttributeError:
                S_list[_s_]=S[_s_]
        json.dump(S_list, S_kick)

    np.save(PATH+mapa+track+kickoff_folder+'grain_ID_kickoff.npy', grain_ID_)
    np.save(PATH+mapa+track+kickoff_folder+'FF_kickoff.npy', FF)
    np.save(PATH+mapa+track+kickoff_folder+'inactive_grains_kickoff.npy', inactive_grains)
    with open(PATH+mapa+track+kickoff_folder+'AG_kickoff.json', 'w') as ag:
        json.dump(AG, ag)
    with open(PATH+mapa+track+kickoff_folder+'IG_kickoff.json', 'w') as ig:
        json.dump(IG, ig)

    counts={}; counts['tm_count']=tm_count ; counts['yp_count']=yp_count ; counts['cut_count']=cut_count; counts['start_step_intermediate'] = i
    counts['h_intermediate']=h-1
    with open(PATH+mapa+track+kickoff_folder+'counts_kickoff.json', 'w') as cnt:
        json.dump(counts, cnt)


def Load_KickOff():
    kickoff_folder = 'kickoff_data/'
    
    faza=np.load(PATH+mapa+track+kickoff_folder+'faza_kickoff.npy')
    rgb_snap=np.load(PATH+mapa+track+kickoff_folder+'rgb_snap_kickoff.npy')
    cas= np.load(PATH+mapa+track+kickoff_folder+'cas_kickoff.npy')
    with open(PATH+mapa+track+kickoff_folder+'nuclei_kickoff.json', 'r') as nuks:              # Writing data of ALL nuclei as .json file, but values must be list NOT np.array !!!
        asc=json.load(nuks)
        asc ={int(k):v for k,v in asc.items()}
        for nuk in asc:
            asc[nuk]['oi']=np.array(asc[nuk]['oi'])
            asc[nuk]['rgb']=np.array(asc[nuk]['rgb'])

    grain_counter = len(list(asc.keys()))

    grain_ID_ = np.load(PATH+mapa+track+kickoff_folder+'grain_ID_kickoff.npy')
    FF = list(np.load(PATH+mapa+track+kickoff_folder+'FF_kickoff.npy'))
    inactive_grains = np.load(PATH+mapa+track+kickoff_folder+'inactive_grains_kickoff.npy').tolist()
           
    with open(PATH+mapa+track+kickoff_folder+'negatives_kickoff.json', 'r') as negs:
        Negatives=json.load(negs)
        Negatives ={int(k):v for k,v in Negatives.items()}
    with open(PATH+mapa+track+kickoff_folder+'S_kickoff.json', 'r') as S_kick:
        S=json.load(S_kick)
        for _s_ in S:
            S[_s_]=np.array(S[_s_])
    with open(PATH+mapa+track+kickoff_folder+'AG_kickoff.json', 'r') as ag:
        AG=json.load(ag)
        AG ={int(k):v for k,v in AG.items()}
    with open(PATH+mapa+track+kickoff_folder+'IG_kickoff.json', 'r') as ig:
        IG=json.load(ig)
        IG ={int(k):v for k,v in IG.items()}
    with open(PATH+mapa+track+kickoff_folder+'counts_kickoff.json', 'r') as cnt:
        counts = json.load(cnt)

    return faza, rgb_snap, cas, asc, Negatives, grain_counter, S, AG, IG, FF, inactive_grains, counts['tm_count'], counts['yp_count'], counts['cut_count'], counts['start_step_intermediate'], grain_ID_, counts['h_intermediate']

def Selection_Mechanism(zrno_ID, smeri_podatkovna, pick_one):
    selection_mechanisms ={
            1: 'Ordered_Grains_n_Directions(zrno_ID, smeri_podatkovna)',                                           # fully ORDERED ........OK
            2: 'Ordered_Grains_n_Directions(zrno_ID, np.flip(smeri_podatkovna))',                            # fully ORDERED, flipped

            3: 'Random_Grains_Ordered_Directions(zrno_ID,  smeri)',                                                    # random selection of grains, ordered selection of directions, from low to high order ......... not ok :( 
            4: 'Random_Grains_Ordered_Directions(zrno_ID, np.flip( smeri_podatkovna.flatten()))', # random selection of grains, ordered selection of directions, from high to low order ......... OK
            
            5: 'Random_Grains_Directions_Segmented(zrno_ID, smeri_podatkovna)',                             # random selection of grains, ordered selection of segments of random directions within, from low to high order ......... not ok :(
            6: 'Random_Grains_Directions_Segmented(zrno_ID, np.flip(smeri_podatkovna))',              # random selection of grains, ordered selection of segments {[001,010,011],  [012,021],  [002,020],  [022]}of random directions within, from high to low order ......... OK

            7: 'Random_Grains_n_Directions(zrno_ID, smeri_podatkovna.flatten())', }                       # fully RANDOM ....... not OK

    selekcija = eval(selection_mechanisms[pick_one])
    return selekcija

""" ====================================  F  I  L  E  S    &   F  O  L  D  E  R  S ========================================================"""
# Path to FE analysis results from Salome Meca 2018

case =           'SLM_2D_Source'   ;  subcase = '0002'
PATH =         'C:/sm-2018-w64-0-3/WORK/'+case+'_Files/post processing database/'+subcase+'/'

#mapa    =      'INTER  time=40, space=8  Z[7-8], X[15-27], Y[12-97], 1500°C, N=12/'
mapa    =       'INTER  time=1, space=8  Z[0-9], X[15-27], Y[12-97], 1500°C, N=12/'


time_factor_folder = '/time_factor_24/time_factor_3'

# Folders for systematic mesh dependency study
tracks_database = {
    'iso': ['2D 1st order Moore, iso field/', '2D 2nd order Moore, iso field/', '3D 1st order Moore, iso field/', '3D 2nd order Moore, iso field/'],
    'real': ['2D 1st order Moore, real field/', '2D 2nd order Moore, real field/', '3D 1st order Moore, real field/', '3D 2nd order Moore, real field/'],
}
track =                                                                               tracks_database['real'][0]

if not os.path.isdir(PATH+mapa+track):
   os.mkdir(PATH+mapa+track)

flashies_RGB =       track+'flashies_RGB/'                          #  Subfolder with time-snap 3D matrices, i.e. flashies
flashies_faza =        track+'flashies_faza/'

cuts_RGB =            track+'cuts_RGB/'                               #  Subfolder with cut 3D matrices, i.e. cuts
cuts_faza =             track+'cuts_faza/'

''' ~~~~~~~~~~~~~~~~ Long_track_CA  ::: Domain Constraints ::: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
z_min =          71
z_max =         72
Z =                  z_max - z_min

X=                  96                     #  [27-15] =  12 cells before space interpolation, and 12*8 = 96 cells after space interpolation
Y=                 128                   #  N=12 ; length of two YPs (pair), each 120 cells long: YP0 [12,27] --->>> 27-12 = 15*8 = 120 + 8
#Y =                32                   #   N=80 ;  

''' ~~~~~~~~~~~~~~~~ Making of  ::: Time (tm_list):::  and  ::: Space (yp_list):::  Partitions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
#y_min =           16       # N = 80 
y_min =             12       # N = 12
y_max =            97

N =                    12         # Number of equally sized segments along Y-axis of 4D matrix

YP = Make_Y_Partitions(y_min, y_max-1, N)
yp_list = [i+'  '+str(YP[i]).replace(' ', '')for i in YP]                                  #  ::: Creation of Y partitions names ::: list of strings  (yp_list)
tm_list = ['TR{0}  [{0},{1}]'.format(i,i+1) for i in range (0,16)]         #  ::: Creation of time ranges names ::: list of strings  (tm_list)
''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''


""" =====================================  I  N  P  U  T    D  A  T  A  ======================================================================="""
# ````` DOMAIN Input variables `````

GDSM = 'new'        #  choose  'old'  or   'new'   , Grain and Directions Selecting Mechanism (GDSM)- see explanation below:
'''
GDSM explanation:

GDSM OLD =    Picks random grain and performs growth in ALL directions by merging function. Here a list
                        seznam premikov has the length of number of given directions.

GDSM NEW =   Picks random grain and performs growth in ONE direction, which can be choosed either in
                        order or randomly. Then, matrix faza is refreshed and the process is repeated for the rest of
                        the directions. '''

pick_selection_mechanism =         7            # pick 7 for fully random grains and directions;  5 for segmented directions

from_beginning =      True                           # True to start from beginning, False to continue from previously saved simulation (KickOff)
save_kickoff =          True                            # if True it saves kickoff parameters (runs Save_KickOff())
save_last_cut_figure = True

avtomatska_nukleacija   =              False
delete_inactive_grains_ID =         False     ;  FF_length = 30     #  Inactive grains are deleted every FF_length time step

save_flashy_as_RGB =               False
save_flashy_as_faza =                False

save_cut_as_RGB =                    True
save_cut_as_faza =                     True

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~ these two parameters in combination with vg (growth velocity)determine CA --->
~~~~~~~~~~~~~~~~~~~~~~~~~~~~ if dt is less than cell, then we get rid of mesh dependency, but then the
~~~~~~~~~~~~~~~~~~~~~~~~~~~~ interpolation time factor should be way larger, which would make CA simulation time very large!'''
FEM_cell_size =                 5e-06
space_factor =                     8

cell=FEM_cell_size/space_factor          # case: SLM_2D_Source

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
FEM_scanning_speed =      800                 # FEM laser scanning speed; unit: MILIMETERS per SECOND [mm/s]
FEM_time_step = FEM_cell_size * 1000 / FEM_scanning_speed           #  time_step to define loads in FEM; unit: SECOND [s]

FEM_time_factor =               5
extra_time_factor =               8

time_factor =     extra_time_factor * FEM_time_factor
'''............................................................................................'''
dt = FEM_time_step / extra_time_factor   #*0.99
'''............................................................................................'''

# Filter Negatives by keys
critical_negind = - 100
negatives_thresh = 500
# Filter Negatives by values
dt_thresh = 500*dt

''' Frequency of new nuclei formed, ß ranges from 0 to 1, high to low, respectively.'''
#ß =   0.995  # low nuclei concentration
ß =   0.95
#ß =     0.97

'''****************************'''
START_step =          0                               # Starting time step number
END_step     =          100                         # Ending time step number

h =                      0                                    #  Relative time step counter

yp_count =         0
tm_count =         1                                    #int((START_step+1)/(time_factor*space_factor))

'''****************************'''
time_shift =   1

cut_count =   0
Cut_Off_Percent =   50                    #  Percent of Stack_2 domain lenght at which this domain should be cut off  ---> the left part is saved (.npy and/or .png)the right goes on to CA and so on and on and on..
cutoff_limit =  64

'''******************* Number of Possible Cubic Unit Cell Random Orientations'''
rp =         45      # Number of possible random alfa, beta, gama choices for grain orientation randomization
'''******************** NEIGHBOURHOOD  DISTANCES   (fixed or randomized)'''
epsilon =   0      # (it's a float between zero and 0.49)
'''********************'''
if epsilon == 0:
   random_distances = False
elif epsilon > 0 and epsilon <= 0.49:
   random_distances = True
R=Dissipate_Distances(epsilon)if random_distances else 1

   
'''======================================================================================================================================================================'''
'''                                                     Moore Neighbourhood - Crystallographic Orientations, Directions & GROUPS                                                                                                                                                    '''
'''======================================================================================================================================================================'''
'''........................................................................................................................................... I. order Moore neighbourhood '''
G1  =  np.array(['001', '00_1', '010', '0_10']) 
G2  =  np.array(['100']) 
G3  =  np.array(['_100']) 
G4  =  np.array(['011', '01_1', '0_11', '0_1_1'])
G5  =  np.array(['101', '10_1', '110', '1_10'])
G6  =  np.array(['_101', '_10_1', '_110', '_1_10'])
G7  =  np.array(['111', '11_1', '1_11', '1_1_1'])
G8  =  np.array(['_111', '_11_1', '_1_11', '_1_1_1'])

'''.......................................................................................................................................... II. order Moore neighbourhood '''
G9  =  np.array(['012', '01_2', '0_12', '0_1_2', '021', '02_1', '0_21', '0_2_1'])
G10 =  np.array(['102', '10_2', '120', '1_20'])
G11 =  np.array(['_102', '_10_2', '_120', '_1_20'])
G12 =  np.array(['201', '20_1', '210', '2_10'])
G13 =  np.array(['_201', '_20_1', '_210', '_2_10'])
G14 = np.array(['002', '00_2', '020', '0_20'])
G15 = np.array(['200'])
G16 = np.array(['_200'])
G17 = np.array(['022', '02_2', '0_22', '0_2_2'])
G18 = np.array(['202', '20_2', '220', '2_20'])
G19 = np.array(['_202', '_20_2', '_220', '_2_20'])
G20 = np.array(['112', '11_2', '1_12', '1_1_2',    '121', '12_1', '1_21', '1_2_1'])
G21 = np.array(['_112', '_11_2', '_1_12', '_1_1_2',    '_121', '_12_1', '_1_21', '_1_2_1'])
G22 = np.array(['211', '21_1', '2_11', '2_1_1'])
G23 = np.array(['_211', '_21_1', '_2_11', '_2_1_1'])
G24 = np.array(['122', '12_2', '1_22', '1_2_2'])
G25 = np.array(['_122', '_12_2', '_1_22', '_1_2_2'])
G26 = np.array(['212', '21_2', '2_12', '2_1_2',    '221', '22_1', '2_21', '2_2_1'])
G27 = np.array(['_212', '_21_2', '_2_12', '_2_1_2',    '_221', '_22_1', '_2_21', '_2_2_1'])
G28 = np.array(['222', '22_2', '2_22', '2_2_2'])
G29 = np.array(['_222', '_22_2', '_2_22', '_2_2_2'])


''' . . . . . . . . . . . . . . . . . . . . . STRUCTURING (subarrays)groups (G1,..)within smeri_database to form segments . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .'''

smeri_database = np.array([np.concatenate((G1,G4))])                                                                             # 1 segment --- 2D , I. order Moore
#smeri_database = np.array([np.concatenate((G1,G4, G9, G14, G17))])                                                       # 1 segment --- 2D , II. order Moore
#smeri_database = np.array([np.concatenate([G1, G4]), np.concatenate([G9, G14, G17]),])                            # 2 segments --- 2D , II. order Moore


#smeri_database = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,9)]))])          # 1 segment --- 3D , I. order Moore
#smeri_database = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,30)]))])        # 1 segment --- 3D , II. order Moore

#smeri_database = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,30)]))])
#                            .reshape((4,31))                                                                                                              # 4 segments (equal size)--- 3D , II. order Moore

#smeri_database = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,9)])),
#                                           np.concatenate(tuple([eval('G{}'.format(i)) for i in range(9,30)])),])       # 2 segments --- 3D , II. order Moore

#smeri_database = np.array([np.concatenate((G1,G2,G3)),np.concatenate((G4,G5,G6)),np.concatenate((G7,G8)),  # 9 segments --- 3D , II. order Moore
#        np.concatenate((G9,G10,G11,G12,G13)), np.concatenate((G14,G15,G16)), np.concatenate((G17,G18,G19)),
#        np.concatenate((G20,G21,G22,G23)), np.concatenate((G24,G25,G26,G27)), np.concatenate((G28,G29)), ])

''' . . . . . . . . . . . . . . . . . . . . >>> smeri_database.flatten() to get ::: 'smeri' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .'''
smeri = smeri_database.flatten()
#smeri = np.array([np.concatenate((G1,G4, G9, G14, G17))]).flatten()

#smeri = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,30)]))]).flatten()

#smeri_database = np.array([np.array(['001']), np.array(['00_1'])])
#smeri = smeri_database.flatten()



# `````` Material properties ``````

''' SS 316 at melting point, i.e. 1650 Kelvin '''
#Lambda = 35
#Rho = 7284
#Cp = 678
#TherDiff = Lambda / (Rho*Cp)
#DTime = cell**2 / TherDiff

''' *********************** Critical Temperatures ******************************************************************************************************************************** '''
Tmelt_Celsius =                             1507               #   Melting point; unit:  deg. Celsius [deg. C]  
Tmelt= Tmelt_Celsius + 273                                  #   Melting point; unit:  KELVIN [K]

''' Absolute liquidus temperature '''
dTliquidus =    50          
''' ********************************************************************************************************************************************************************************* '''



#  ~~~~ AUTO Nucleation (AN)~~~~
if avtomatska_nukleacija:
   vg = np.zeros((Z,X,Y))        # matrika hitrosti rasti
   NP = np.vectorize(nukleacija_povrsina)
   NV = np.vectorize(nukleacija_volumen)
   if from_beginning:
      faza = np.zeros((Z,X,Y))  # fazna matrika
      cas = np.zeros((Z,X,Y))   # časovna matrika
      Negatives = {-1:0}; asc = {}; grain_counter = 0; S = {}; inactive_grains=[]
      IG = {}     #  IG (saves time-steps when indivudual grain didn't grow IN ANY DIRECTION, IG stands for Inactive Grains)
      AG = {}   #  AG (saves time-steps when indivudual grain did grow IN ANY DIRECTION, AG stands for Active Grains)
      FF = []
      
   elif not from_beginning:
      faza, rgb_snap, cas, asc, Negatives, grain_counter, S, AG, IG, FF, inactive_grains, tm_count, yp_count, cut_count, START_step, grain_ID_, h  = Load_KickOff()
      #if negind<=critical_negind:
            #Negatives ={key:val for key, val in Negatives.items() if key < (negind+negatives_thresh)}
      Negatives ={key:val for key, val in Negatives.items() if val <  dt_thresh}

#yps_test=Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], START_step)
#cutoff_limit = int(yps_test.shape[2]*Cut_Off_Percent/100)



#  ~~~~ Manual Nucleation (MN)~~~~
if not avtomatska_nukleacija:
    '''||||||||||||||||||||||||||| nucleation - manual |||||||||||||||||||||||||| '''
    Z,X,Y = 1, 100, 100                                              # Size of domain in terms of cells in Z,X,Y directions, respectively, for testing and code development
    faza = np.zeros((Z,X,Y))                               # fazna matrika
    cas = np.zeros((Z,X,Y))                                # časovna matrika
    vg = np.zeros((Z,X,Y))                                  # matrika hitrosti rasti
    vg =   1                                                               # Value of homogeneous 'growing velocity' field, a.u., for testing and code development
    cell =  1                                                              # for mesh dependency development (MDD)
    dt = cell/8
    T=1; T_next=0
    Negatives = {-1:0}; asc = {}; grain_counter = 0; S = {}; inactive_grains=[]; FF = []
    critical_negind = - 2000
    negatives_thresh = 15

    M = { 1: {'ß':(0, 50, 50), 'Ł': (0,0,0)},
               #2: {'ß':(11, 11, 11), 'Ł': (1,0,0)},
               #3: {'ß':(9, 9, 9), 'Ł': (1,1,0)},                    #  data of manually created nuclei, 'ß' are the (Z,X,Y) coordinates, 'Ł' are tilting parameters (x,y,gama)
               #4: {'ß':(0, 100, 100), 'Ł': (0,0,45*0.5)},
               #5: {'ß':(0, 125, 125), 'Ł': (0,0,45*0.75)},
               #6: {'ß':(0, 150, 150), 'Ł': (0,0,45*0.875)},
               #7: {'ß':(0, 175, 175), 'Ł': (0,0,45)},

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

    Selection = Selection_Mechanism(grain_ID, smeri_database, pick_selection_mechanism)
    
    cas[np.isin(faza, grain_ID, invert=False)] = -1
    taula=0;likvid=0; grain_counter=len(grain_ID)



# Writing the log..
with open(PATH+mapa+track+'Logfile.txt', 'w')as cuttxt:
    cuttxt.write(100*'*'+'\n'+'z_min = '+str(z_min)+' ,  z_max = '+str(z_max)+' ,  Z = '+str(Z)+' ,\n'+
                 'X = '+str(X)+' ,  y_min = '+str(y_min)+' ,  y_max = '+str(y_max)+' ,  Y = '+str(Y)+' ,\n'+
                 'N = '+str(N)+' ,\n\n'+

                 'FEM_time_step = '+str(FEM_time_step)+' sec.'+' ,\n'+
                 'dt = '+str(dt)+' sec. ,  FEM_scanning_speed = '+str(FEM_scanning_speed)+' mm/s'+' ,\n'+
                 'FEM_cell_size = '+str(1e6*FEM_cell_size)+u'\u00B5m ,  cell = '+str(1e6*cell)+u'\u00B5m'+' ,\n'+
                 'space_factor = '+str(space_factor)+' ,  FEM_time_factor = '+str(FEM_time_factor)+' ,  extra_time_factor = '+str(extra_time_factor)+' ,\n\n'+

                 'START_step = '+str(START_step)+' ,\n\n'+

                 'smeri = '+str(list(smeri))+' ,\n\n'+

                 'ß = '+str(ß)+' ,\n\n'+

                 'Tmelt_Celsius = '+str(Tmelt_Celsius)+u'\N{DEGREE SIGN}C ,  dTliquidus = '+str(dTliquidus)+u'\N{DEGREE SIGN}C'+' ,\n'+
                 'delete_inactive_grains_ID = '+str(delete_inactive_grains_ID)+'\n\n'+
                 100*'*'+'\n\n')



'''  M P I  $$$$$$$$$$$$$$$$$$$  MPI  tryin' to make the code run on several cores using multiprocessing module  MPI  $$$$$$$$$$$$  M P I  '''

import multiprocessing
import concurrent.futures


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ... STARTING the CA ..
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Growth (time FOR-LOOP, grains FOR-LOOP, directions FOR-LOOP)~~~~~~~~~~~~~~~~~~~~~~~~~
start_time = time.time()

for i in range(START_step, END_step+1):
    step_time_start = time.time(); h+=1; print(i)
    # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   1.   N  U  C  L  E  A  T  I  O  N   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    
    if avtomatska_nukleacija:
        ''' avtomatska nukleacija '''

        T = Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], i)[z_min:z_max]
        taljenje(T, Tmelt)                              # condition for melting  ----> taula matrix of cell values;  if value -1 then solid (or powder)elif value 0 then liquid
        
        if np.all(taula[:,:,:cutoff_limit]==-1):   # CUT condition and consequences
            print(); print(50*'*',' CUT! ',50*'*')
            print('Cut_Off_Percent = ',Cut_Off_Percent,' %    , cutoff limit = ',cutoff_limit)
            cut_text = 'Time step number  '+str(i)+',  real time: '+str(round(1000*dt*h, 3))+' msec.'
            print(cut_text)
            print(106*'*'); print()

            with open(PATH+mapa+track+'cut_data.txt', 'a')as cuttxt:
                cuttxt.write(cut_text+'\n')

            if save_cut_as_RGB:
                try:
                    np.save(PATH+mapa+cuts_RGB+'cut_RGB_'+str(cut_count)+'.npy', rgb_snap[:,:,:cutoff_limit])       # Saves the first half of RGB snap as .NPY
                except FileNotFoundError:
                    os.mkdir(PATH+mapa+cuts_RGB)
                    np.save(PATH+mapa+cuts_RGB+'cut_RGB_'+str(cut_count)+'.npy', rgb_snap[:,:,:cutoff_limit])
            if save_cut_as_faza:
                try:
                    np.save(PATH+mapa+cuts_faza+'cut_faza_'+str(cut_count)+'.npy', faza[:,:,:cutoff_limit])                 # Saves the first half of faza snap as .NPY
                except FileNotFoundError:
                    os.mkdir(PATH+mapa+cuts_faza)
                    np.save(PATH+mapa+cuts_faza+'cut_faza_'+str(cut_count)+'.npy', faza[:,:,:cutoff_limit])
                    
            yp_count+=1
            cut_count+=1
            
            faza = np.dstack((faza[:,:,cutoff_limit:],np.zeros((Z,X,cutoff_limit) )))
            for s in smeri:
                S[s] = np.dstack((S[s][:,:,cutoff_limit:],np.zeros((Z,X,cutoff_limit) )))
            FF=[]
            cas = np.dstack((cas[:,:,cutoff_limit:],np.zeros((Z,X,cutoff_limit) )))
                
            T = Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], i)[z_min:z_max]
            taljenje(T, Tmelt)
                    
        liquidus(T,Tmelt+dTliquidus)              # absolute liquidus line
        
        dTt =  T  -  Tmelt                                  # undercooling [K]
        
        interface = NP(dTt)                                              
        bulk = NV(dTt)
        live= nakljucje(taula)

        try:
            T_next = Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], i+time_shift)[z_min:z_max]
        except FileNotFoundError:
            print(); print('FileNotFound Exception !'); print()
            tm_count+=1
            T_next = Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], i+time_shift)[z_min:z_max]
        except IndexError:
            raise IndexError('Please, correct the time range!')

        ''' ======================================= N  U  C  L  E  A  T  I  O  N ============================================================== '''
        new_grains_ID = []
        
        for k in range(Z):
         for ii in range(X):
            for j in range(Y):
               if faza[k][ii][j]==0 and (ß<live[k][ii][j]<interface[k][ii][j] or (live[k][ii][j]<bulk[k][ii][j] and bulk[k][ii][j]>ß)) and T_next[k][ii][j] < T[k][ii][j]:              
               
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
                  oi = Rotate_the_Cube_Z(cub_xy, gama)
                  asc[grain_counter] ={'oi': oi, 'alfa':alfa, 'beta':beta, 'gama':gama, 'rgb':rgb, 'coords':(i,k,ii,j), 'temp': T[k,ii,j]-273, }    # ALL data about nuclei
        ''' ========================================================================================================================== '''
        try:
            Selection = Selection_Mechanism(grain_ID_, smeri_database, pick_selection_mechanism)   
            grain_ID = grain_ID_.copy()
        except NameError:
            grain_ID = np.array(list(asc.keys())); print('dolžina grain_ID: ', len(grain_ID))
            Selection = Selection_Mechanism(grain_ID, smeri_database, pick_selection_mechanism) 
        
        cas[np.isin(faza, new_grains_ID, invert=False)] = -1          # vrednost časovne matrike vseh novih nukleusov je -1
        vg = growth_speed(dTt)

    elif not avtomatska_nukleacija:
        Selection = Selection_Mechanism(grain_ID, smeri_database, pick_selection_mechanism)
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 

    # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   2.   T  I  M  E  -  MPI   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    if __name__ == '__main__':
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.submit(mef.NN, time_counter, Negatives)
            executor.submit(mef.MM, smeri, Negatives, faza, grain_ID, cas, S)
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 

    # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   3.   P  H  A  S   E   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

    F = faza.copy()
    FF.append(faza); r=[]
    
    if __name__ == '__main__':
        with concurrent.futures.ThreadPoolExecutor() as executor:
            r.append(executor.submit(mef.GG, GDSM, Selection, faza, asc, vg, S, taula, T_next, T, likvid, R, cell, W, cas, random_distances))

    SZP=[i.result() for i in r]
 
    faza = SZP[0][0]
    cas = SZP[0][1]

    '''
    if GDSM == 'new':
        #*********************************************************************** NEW ***************************************************************************************************
        grain_size = np.count_nonzero(faza == grain)
        asc[grain]['grain size'] = grain_size 
        # ----------------------------------------- this function selects inactive grains (i.e. grains, which no longer grow -
        # ----------------------------------------- either domain limits were reached and/or has reached another grain)
       
        if delete_inactive_grains_ID and len(FF)== FF_length:
            active=np.isin(faza-FF[0], grain)
            if np.all(active==False):
               inactive_grains.append(grain)
       
               #IG[grain].append(i)
            #else:
               #AG[grain].append(i)
        # -----------------------------------------
        #************************************************************************* NEW ***************************************************************************************************
    '''
    if GDSM == 'old':
        #------------------------------------------------------------------------------------  old  ---------------------------------------------------------------------------------------------------------------
        total=merging(seznam_premikov[0], seznam_premikov[1])
        if len(smeri)> 2:
            for move in range(2, len(smeri)):
                total=merging(total, seznam_premikov[move])
        faza = total.copy()
        grain_size = np.count_nonzero(faza == grain)
        asc[grain]['grain size'] = grain_size
        '''
        with open(PATH+mapa+track+'nuclei_data.json', 'w') as nuks:              # Writing data of ALL nuclei as .json file, but values must be list NOT np.array !!!
            asc_list =asc.copy()
            for nuk in asc:
                asc_list[nuk]['oi']=asc[nuk]['oi'].tolist()                                
                asc_list[nuk]['rgb']=asc[nuk]['rgb'].tolist()
            json.dump(asc_list, nuks)
        ''' 
        # ----------------------------------------- this function selects inactive grains (i.e. grains, which no longer grow - either domain limits were reached and/or has reached another grain)
        if delete_inactive_grains_ID and len(FF)==FF_length:
            active=np.isin(faza-FF[0], grain)
            if np.all(active==False):
               inactive_grains.append(grain)
               #IG[grain].append(i)
            #else:
               #AG[grain].append(i)
        # -----------------------------------------
        #------------------------------------------------------------------------------------  old  ---------------------------------------------------------------------------------------------------------------

    grain_ID = np.array(list(asc.keys()))
    #grain_ID_ = grain_ID
    
    if delete_inactive_grains_ID and len(FF)== FF_length:

        '''............................................ erasing inactive grains from grain_ID register ........................................'''
        #for ggg in grain_ID:
         #   active=np.isin(faza-FF[0], ggg)  
          #  if np.all(active==False):
           #     inactive_grains.append(ggg)
                
        grain_ID_ = np.array(list(set(grain_ID)-set(inactive_grains)))    #  array of active grains (grain_ID_)
        
        '''.......................................................................................................................................................................'''
        del FF[0]

        '''............................................ erasing Negatives items with values less than dt_thresh ........................................'''

        if avtomatska_nukleacija:
            Negatives ={key:val for key, val in Negatives.items() if val < dt_thresh}
            #Negatives ={key:val for key, val in Negatives.items() if val < 1000*dt}
            
        else:
            if negind<=critical_negind:
                Negatives ={key:val for key, val in Negatives.items() if key < (negind+negatives_thresh)}
    
        '''.......................................................................................................................................................................'''

    else:
        grain_ID_ = grain_ID
    
    if not np.all(faza==F):
        #negind -=1
        negind = min(list(Negatives.keys()))- 1             #  Here, negind should be defined as absolute counter of negative numbers to avoid erasing existing Negatives items
        for s in S:
            S[s][np.isin(faza-F, grain_ID_)]= negind
        Negatives[negind]=0

    if GDSM == 'old':
        #--------------------------------------------------------------------------------------------- TIME MATRIX part - OLD --------------------<<<
        times = np.array(list(S.values()))
        cas_total = merging(times[0], times[1])
        if len(smeri)> 2:
            for move in range(2, len(smeri)):
                cas_total = merging(cas_total, times[move])
        try:
            cas = cas_total.copy()
        except NameError:
            pass
        #--------------------------------------------------------------------------------------------- TIME MATRIX part - OLD --------------------<<<
    
    rgb_snap = np.zeros((faza.shape[0], faza.shape[1], faza.shape[2], 3)).astype('int')
    for zrno in asc:
        rgb_snap[faza==zrno]=asc[zrno]['rgb']
    
    if save_flashy_as_RGB:
        try:
            np.save(PATH+mapa+flashies_RGB+'flashy_RGB_'+str(i)+'.npy', rgb_snap)
        except FileNotFoundError:
            os.mkdir(PATH+mapa+flashies_RGB)
            np.save(PATH+mapa+flashies_RGB+'flashy_RGB_'+str(i)+'.npy', rgb_snap)
  
    if save_flashy_as_faza:
        try:
            np.save(PATH+mapa+flashies_faza+'flashy_faza_'+str(i)+'.npy', faza)
        except FileNotFoundError:
            os.mkdir(PATH+mapa+flashies_faza)
            np.save(PATH+mapa+flashies_faza+'flashy_faza_'+str(i)+'.npy', faza)


    step_cpu = round(time.time()- step_time_start, 3)
    try:
        rate = round((time.time()- start_time)/h, 4)
    except ZeroDivisionError:
        rate = 'inf.'

    print('TIME STEP # ',i,'    |    current step CPU: ',step_cpu,' s    |    RATE = ',rate,' s / step    |    active grains / ALL grains :  ',str(grain_ID_.shape[0]),' / ',grain_counter)
    print(34*' ','length = ',round(FEM_scanning_speed*dt*h, 4),u'\u00B5m              |      time = ',round(1000*dt*h, 4),' msec.','      |    inactive grains: ',len(list(set(inactive_grains))))
    print(30*' .  . ')
    
    if not bool(grain_ID.shape[0]) and grain_counter:
       break
    

process_time = time.time() - start_time
print();print(); print()
print(20*' ',100*'*')
print(25*' ','DONE!!  Total computing time =  ',round(process_time, 3),'  seconds.   (RATE = ', round(process_time/h, 4),' s / step)'); print(25*' ',100*'-')
print(25*' ','Printed length real dimension = ',round(FEM_scanning_speed*dt*h, 5),u'\u00B5m      |    Printing time = ',round(1000*dt*h, 5),' msec.')
print(25*' ',100*'-')
print(65*' ','Number of grains =  ',grain_counter,' !'); print(20*' ',100*'*'); print(); print(); print()

# ......................................... PRIKAZ in IZPIS REZULTATOV ..................................................#

#plt.imshow(rgb_snap[0]); plt.figure(); plt.imshow(faza[0])
np.save(PATH+mapa+cuts_RGB+'cut_RGB_'+str(cut_count)+'.npy', rgb_snap[:,:,cutoff_limit:,:])
np.save(PATH+mapa+cuts_faza+'cut_faza_'+str(cut_count)+'.npy', faza[:,:,cutoff_limit:])
if save_kickoff:
    Save_KickOff()
if save_last_cut_figure:
    np.save(PATH+mapa+cuts_RGB+'cut_RGB_'+str(cut_count)+'.npy', rgb_snap[:,:,:cutoff_limit,:])
    np.save(PATH+mapa+cuts_RGB+'cut_RGB_'+str(cut_count+1)+'.npy', rgb_snap[:,:,cutoff_limit:,:])
    np.save(PATH+mapa+cuts_faza+'cut_faza_'+str(cut_count)+'.npy', faza[:,:,:cutoff_limit])
    np.save(PATH+mapa+cuts_faza+'cut_faza_'+str(cut_count+1)+'.npy', faza[:,:,cutoff_limit:])

''' ---------------------------------------------------------------------------------- the end -------------------------------------------------------------------------------------------------------------- '''






