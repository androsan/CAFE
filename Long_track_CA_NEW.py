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
import psutil
plt.ion()
#np.random.seed(909175683)

""" =====================================  F  U  N  C  T  I  O  N  S ================================================================"""
def f(n):
    f=np.load('C:/sm-2018-w64-0-3/WORK/SLM_2D_Source_Files/post processing database/0002/'+
              'INTER  time=1, space=8  Z[0-9], X[15-27], Y[12-97], 1500°C, N=12/2D 1st order Moore, real field/flashies_RGB/flashy_RGB_{}.npy'.format(n))
    return f


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
    ns_max =                      5e10                                          #  Number of new nuclei at optimal dT; unit: per SQUARE METER  [m^-2]
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
    #x = np.random.uniform(0,1)
    #y = np.random.uniform(0,x)
    xran1 = np.random.randint(0, posibilities+1)
    xran2 = np.random.randint(0, posibilities+1)
    xran3 = np.random.randint(0, posibilities+1)
    xran = (xran1+xran2+xran3)/3  
    x = xran / posibilities
    y = np.random.randint(0, xran+1)/posibilities
    return x,y

def merging(prva,druga):
    print('MERGING!')
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
                #print(20*'-')
                #print('zrno: ',grain)
                #print('smer: ',s)
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

pick_selection_mechanism =         5            # pick 7 for fully random grains and directions;  5 for segmented directions

from_beginning =      True                           # True to start from beginning, False to continue from previously saved simulation (KickOff)
save_kickoff =          True                            # if True it saves kickoff parameters (runs Save_KickOff())
save_last_cut_figure = True

avtomatska_nukleacija   =              False
delete_inactive_grains_ID =         False     ;  FF_length = 30     #  Inactive grains are deleted every FF_length time step

save_flashy_as_RGB =               False
save_flashy_as_faza =                False

save_cut_as_RGB =                    True
save_cut_as_faza =                     True

run_1st = True
run_2nd = True
run_3rd = True

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
#ß =   0.95
ß =     0.95

'''****************************'''
START_step =          0                               # Starting time step number
END_step     =          50                             # Ending time step number

yp_count =         0
tm_count =         1                                    #int((START_step+1)/(time_factor*space_factor))

'''****************************'''
time_shift =   1

cut_count =   0
Cut_Off_Percent =   50                    #  Percent of Stack_2 domain lenght at which this domain should be cut off  ---> the left part is saved (.npy and/or .png)the right goes on to CA and so on and on and on..
cutoff_limit =  64

h=0                                                  #  Relative time step counter
'''======================================================================================================================================================================'''
'''                                                     Moore Neighbourhood - Crystallographic Orientations, Directions & GROUPS                                                                                                                                                    '''
'''======================================================================================================================================================================'''
# Crystallographic GROUPS Booleans

# 2D --- Moore 1st Order Neighbourhood:
'''***********************'''
Moore_I_2D = True
'''***********************'''
group_1 = True if Moore_I_2D else False                 # .............................. [001], [010]
group_4 = True if Moore_I_2D else False                 # .............................. [011]

# 2D --- Moore 2nd Order Neighbourhood:
'''***********************'''
Moore_II_2D = True
'''***********************'''
group_9 =   True if Moore_II_2D else False              # .............................. [012], [021]
group_14 = True if Moore_II_2D else False              # .............................. [002], [020]
group_17 = True if Moore_II_2D else False              # .............................. [022]

# 3D === Moore 1st Order Neighbourhood:
'''***********************'''
Moore_I_3D = True
'''***********************'''
group_2 = True if Moore_I_3D else False                # .............................. [100]
group_3 = True if Moore_I_3D else False                # .............................. [_100]
group_5 = True if Moore_I_3D else False                # .............................. [101], [110]
group_6 = True if Moore_I_3D else False                # .............................. [_101], [_110]
group_7 = True if Moore_I_3D else False                # .............................. [111]
group_8 = True if Moore_I_3D else False                # .............................. [_111]

# 3D === Moore 2nd Order Neighbourhood:
'''***********************'''
Moore_II_3D = True
'''***********************'''
group_10 = True if Moore_II_3D else False              # .............................. [102], [120]
group_11 = True if Moore_II_3D else False              # .............................. [_102], [_120]
group_12 = True if Moore_II_3D else False              # .............................. [201], [210]
group_13 = True if Moore_II_3D else False              # .............................. [_201], [_210]
group_15 = True if Moore_II_3D else False              # .............................. [200]
group_16 = True if Moore_II_3D else False              # .............................. [_200]
group_18 = True if Moore_II_3D else False              # .............................. [202], [220]
group_19 = True if Moore_II_3D else False              # .............................. [_202], [_220]
group_20 = True if Moore_II_3D else False              # .............................. [112], [121]
group_21 = True if Moore_II_3D else False              # .............................. [_112], [_121]
group_22 = True if Moore_II_3D else False              # .............................. [211]
group_23 = True if Moore_II_3D else False              # .............................. [_211]
group_24 = True if Moore_II_3D else False              # .............................. [122]
group_25 = True if Moore_II_3D else False              # .............................. [_122]
group_26 = True if Moore_II_3D else False              # .............................. [212], [221]
group_27 = True if Moore_II_3D else False              # .............................. [_212], [_221]
group_28 = True if Moore_II_3D else False              # .............................. [222]
group_29 = True if Moore_II_3D else False              # .............................. [_222]


'''........................................................................................................................................... I. order Moore neighbourhood '''
G1  =  np.array(['001', '00_1', '010', '0_10']) if group_1 else None
G2  =  np.array(['100']) if group_2 else None
G3  =  np.array(['_100']) if group_3 else None
G4  =  np.array(['011', '01_1', '0_11', '0_1_1'])if group_4 else None
G5  =  np.array(['101', '10_1', '110', '1_10'])if group_5 else None
G6  =  np.array(['_101', '_10_1', '_110', '_1_10'])if group_6 else None
G7  =  np.array(['111', '11_1', '1_11', '1_1_1'])if group_7 else None
G8  =  np.array(['_111', '_11_1', '_1_11', '_1_1_1'])if group_8 else None

'''.......................................................................................................................................... II. order Moore neighbourhood '''
G9  =  np.array(['012', '01_2', '0_12', '0_1_2', '021', '02_1', '0_21', '0_2_1'])if group_9 else None
G10 =  np.array(['102', '10_2', '120', '1_20'])if group_10 else None
G11 =  np.array(['_102', '_10_2', '_120', '_1_20'])if group_11 else None
G12 =  np.array(['201', '20_1', '210', '2_10'])if group_12 else None
G13 =  np.array(['_201', '_20_1', '_210', '_2_10'])if group_13 else None
G14 = np.array(['002', '00_2', '020', '0_20'])if group_14 else None
G15 = np.array(['200'])if group_15 else None
G16 = np.array(['_200'])if group_16 else None
G17 = np.array(['022', '02_2', '0_22', '0_2_2'])if group_17 else None
G18 = np.array(['202', '20_2', '220', '2_20'])if group_18 else None
G19 = np.array(['_202', '_20_2', '_220', '_2_20'])if group_19 else None
G20 = np.array(['112', '11_2', '1_12', '1_1_2',    '121', '12_1', '1_21', '1_2_1'])if group_20 else None
G21 = np.array(['_112', '_11_2', '_1_12', '_1_1_2',    '_121', '_12_1', '_1_21', '_1_2_1'])if group_21 else None
G22 = np.array(['211', '21_1', '2_11', '2_1_1'])if group_22 else None
G23 = np.array(['_211', '_21_1', '_2_11', '_2_1_1'])if group_23 else None
G24 = np.array(['122', '12_2', '1_22', '1_2_2'])if group_24 else None
G25 = np.array(['_122', '_12_2', '_1_22', '_1_2_2'])if group_25 else None

G26 = np.array(['212', '21_2', '2_12', '2_1_2',    '221', '22_1', '2_21', '2_2_1'])if group_26 else None
G27 = np.array(['_212', '_21_2', '_2_12', '_2_1_2',    '_221', '_22_1', '_2_21', '_2_2_1'])if group_27 else None
G28 = np.array(['222', '22_2', '2_22', '2_2_2'])if group_28 else None
G29 = np.array(['_222', '_22_2', '_2_22', '_2_2_2'])if group_29 else None

''' CONSTRAINTS of the Groups '''

constrains_of_group_9 =   False
constrains_of_group_14 = False
constrains_of_group_17 = False

''' . . . . . . . . . . . . . . . . . . . . . STRUCTURING (subarrays)groups (G1,..)within smeri_database to form segments . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .'''

#smeri_database = np.array([np.concatenate((G1,))])
#smeri_database = np.array([np.concatenate((G1,G4))])                                                                             # 1 segment --- 2D , I. order Moore
smeri_database = np.array([np.concatenate((G1,G4, G9, G14, G17))])                                                       # 1 segment --- 2D , II. order Moore
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
#smeri_database = np.array([['001', '00_1',]])
#smeri_database = np.array([['001', '00_1', '010', '0_10', ]])
smeri = smeri_database.flatten()
#smeri = np.array([np.concatenate((G1,G4, G9, G14, G17))]).flatten()

#smeri = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,9)]))]).flatten()
#smeri = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,30)]))]).flatten()


"""
smeri = np.array(['001', '00_1', '010', '0_10', '011', '01_1', '0_11', '0_1_1',
                           '012', '01_2', '0_12', '0_1_2', '021', '02_1', '0_21', '0_2_1',
                           '002', '00_2', '020', '0_20',
                           '022', '02_2', '0_22', '0_2_2',
                   ])
                   
smeri_database = np.array([

                           np.array(['001', '00_1', '010', '0_10',       # 1st (cell * 1)- - - z= 0 .... 4 sites (GROUP 1)
                           #'100',                                      # 1st (cell * 1)- - - z= 1 .... 1 site (GROUP 2) 
                           #'_100',                                   # 1st (cell * 1)- - - z= -1 .... 1 site (GROUP 3)
                           #.....................................................................................................................................................................===> TOTAL:  6 sites

                  
                           '011', '01_1', '0_11', '0_1_1']),                      # 2nd (cell * sqrt[2] )- - - z= 0 .... 4 sites (GROUP 4)
                           #'101', '10_1', '110', '1_10',                          # 2nd (cell * sqrt[2] )- - - z= 1 .... 4 sites (GROUP 5)
                           #'_101', '_10_1', '_110', '_1_10',              # 2nd (cell * sqrt[2] )- - - z= -1 .... 4 sites (GROUP 6)
                           #.....................................................................................................................................................................===> TOTAL:  12 sites


                           #'111', '11_1', '1_11', '1_1_1',                    # 3rd (cell * sqrt[3])- - - z= 1 .... 4 sites (GROUP 7)
                           #'_111', '_11_1', '_1_11', '_1_1_1',        # 3rd (cell * sqrt[3])- - - z= -1 .... 4 sites (GROUP 8)
                           #.....................................................................................................................................................................===> TOTAL:  8 sites

                           
                           np.array(['012',
                           '01_2',
                           '0_12',
                           '0_1_2',

                           '021',
                           '02_1',
                           '0_21',
                           '0_2_1',
                           


                           # 4th (cell * sqrt[5])- - - z= 0 .... 8 sites (GROUP 9)
                           #'102', '10_2', '120', '1_20',                                                                # 4th (cell * sqrt[5])- - - z= 1 .... 4 sites (GROUP 10)
                           #'_102', '_10_2', '_120', '_1_20',                                                   # 4th (cell * sqrt[5])- - - z= -1 .... 4 sites (GROUP 11)
                           #'201', '20_1', '210', '2_10',                                                               # 4th (cell * sqrt[5])- - - z= 2 .... 4 sites (GROUP 12)
                           #'_201', '_20_1', '_210', '_2_10',                                                   # 4th (cell * sqrt[5])- - - z= -2 .... 4 sites(GROUP 13)
                           #.....................................................................................................................................................................===> TOTAL:  24 sites
                                                
                           #np.array([
                           '002',
                           '00_2',
                           '020',
                           '0_20',
                                     #]),                                                                                               # 5th (cell * 2)- - - z= 0 .... 4 sites (GROUP 14)
                        
                           #'200',                                                                                                # 5th (cell * 2)- - - z= 2 .... 1 site (GROUP 15)
                           #'_200',                                                                                             # 5th (cell * 2)- - - z= -2 .... 1 site (GROUP 16)
                           #.....................................................................................................................................................................===> TOTAL:  6 sites

                        
                           #np.array([
                           '022',
                           '02_2',
                           '0_22',
                           '0_2_2',
                                        ]),                                                                                         # 6th (cell * sqrt[8])- - - z= 0 .... 4 sites (GROUP 17)
                           
                           #'202', '20_2', '220', '2_20',                                                               # 6th (cell * sqrt[8])- - - z= 2 .... 4 sites (GROUP 18)
                           #'_202', '_20_2', '_220', '_2_20',                                                   # 6th (cell * sqrt[8])- - - z= -2 .... 4 sites(GROUP 19)
                           #.....................................................................................................................................................................===> TOTAL:  12 sites

                           
                           #'112', '11_2', '1_12', '1_1_2',    '121', '12_1', '1_21', '1_2_1',                                    # 7th (cell * sqrt[6])- - - z= 1 .... 8 sites (GROUP 20)
                           #'_112', '_11_2', '_1_12', '_1_1_2',    '_121', '_12_1', '_1_21', '_1_2_1',           # 7th (cell * sqrt[6])- - - z= -1 .... 8 sites(GROUP 21)
                           #'211', '21_1', '2_11', '2_1_1',                                                                                      # 7th (cell * sqrt[6])- - - z= 2.... 4 sites (GROUP 22)
                           #'_211', '_21_1', '_2_11', '_2_1_1',                                                                          # 7th (cell * sqrt[6])- - - z= -2.... 4 sites (GROUP 23)
                           #.....................................................................................................................................................................===> TOTAL:  24 sites
                  

                           #'122', '12_2', '1_22', '1_2_2',                                                                                        # 8th (cell * 3)- - - z= 1 .... 4 sites (GROUP 24) 
                           #'_122', '_12_2', '_1_22', '_1_2_2',                                                                           # 8th (cell * 3)- - - z= -1 .... 4 sites(GROUP 25)
                           #'212', '21_2', '2_12', '2_1_2',    '221', '22_1', '2_21', '2_2_1',                                     # 8th (cell * 3)- - - z= 2 .... 8 sites(GROUP 26)
                           #'_212', '_21_2', '_2_12', '_2_1_2',    '_221', '_22_1', '_2_21', '_2_2_1',             # 8th (cell * 3)- - - z=-2 .... 8 sites(GROUP 27)
                           #.....................................................................................................................................................................===> TOTAL:  24 sites


                           #'222', '22_2', '2_22', '2_2_2',                                        # 9th (cell * sqrt[12])- - - z= 2 .... 4 sites (GROUP 28)
                           #'_222', '_22_2', '_2_22', '_2_2_2',                            # 9th (cell * sqrt[12])- - - z= -2 .... 4 sites (GROUP 29)
                           #.....................................................................................................................................................................===> TOTAL:  8 sites

                           # === TOTAL SITES === 124
                  ])
"""

# `````` Material properties ``````

''' SS 316 at melting point, i.e. 1650 Kelvin '''
Lambda = 35
Rho = 7284
Cp = 678
TherDiff = Lambda / (Rho*Cp)
DTime = cell**2 / TherDiff

''' Melting point '''
Tmelt_Celsius =                            1507             #     case: SLM_2D_Source

''' Absolute liquidus temperature '''
dTliquidus =    50          

Tmelt= Tmelt_Celsius + 273                                 #   Melting point; unit:  KELVIN [K]

''' Number of Possible Cubic Unit Cell Random Orientations '''
rp = 100                                                                #  Number of possible random alfa, beta, gama choices for grain orientation randomization

''' ..........................NEIGHBOURHOOD  DISTANCES   (fixed or randomized)......................... '''

'''********************'''
epsilon =   0     # (it's a float between zero and 0.49)
'''********************'''

if epsilon == 0:
   random_distances = False
elif epsilon > 0 and epsilon <= 0.49:
   random_distances = True


def Dsr_1st(r):
   if random_distances:
      cell_1st = cell * (0.5 + 1.158312*r)                              #"""------------------1st shell ::: [001], [010], [100] ::: 6 neighbours ------------------"""
   else:
      cell_1st=cell
   return cell_1st

def Dsr_2nd(r):
   if random_distances:
      #cell_2nd=cell*math.sqrt(2)*(0.5 + 1.081*r)           #''' ------------------2nd shell ::: [011], [101], [110] ::: 12 neighbours ------------------ '''
      cell_2nd=cell*(0.5*math.sqrt(2)+1.4723*r)
   else:
      cell_2nd=cell*math.sqrt(2)
   return cell_2nd

def Dsr_3rd(r):
   if random_distances:
      cell_3rd= None         # !!! to be constructed                #''' ------------------3rd shell ::: [111] ::: 8 neighbours --------------------------- '''
   else:
      cell_3rd=cell*math.sqrt(3)                                      
   return cell_3rd

def Dsr_4th(r):
   if random_distances:
      cell_4th=cell*(1.581 + 1.376*r)                                #''' ------------------4th shell ::: [012], [021], [102], [120], [201], [210] ::: 24 neighbours ------------------ '''
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
      cell_6th=cell*(2.12132 + 1.44939*r)                          #''' ------------------6th shell ::: [022], [202], [220] ::: 12 neighbours ------------------ '''
   else:
      cell_6th=cell*math.sqrt(8)
   return cell_6th

def Dsr_7th(r):
   if random_distances:
      cell_7th= None         # !!! to be constructed                  #''' ------------------7th shell ::: [112], [121], [211] ::: 24 neighbours ------------------ '''
   else:
      cell_7th=cell*math.sqrt(6)                                      
   return cell_7th

def Dsr_8th(r):
   if random_distances:
      cell_8th= None         # !!! to be constructed                  #''' ------------------8th shell ::: [122], [212], [221] ::: 24 neighbours ------------------ '''
   else:
      cell_8th=cell*3                                 
   return cell_8th

def Dsr_9th(r):
   if random_distances:
      cell_9th= None         # !!! to be constructed                  #''' ------------------9th shell ::: [222] ::: 8 neighbours ------------------------------ '''
   else:
      cell_9th=cell*math.sqrt(12)                                 
   return cell_9th


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

''' ======================================= N  U  C  L  E  A  T  I  O  N ============================================================== '''
#  ~~~~ Manual Nucleation (MN)~~~~
if not avtomatska_nukleacija:
    '''||||||||||||||||||||||||||| nucleation - manual |||||||||||||||||||||||||| '''
    Z,X,Y = 1, 101, 101                                              # Size of domain in terms of cells in Z,X,Y directions, respectively, for testing and code development
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
                   2: {'ß':(0, 40, 40), 'Ł': (0,0,45)},
                   3: {'ß':(0, 60, 60), 'Ł': (1,1,0)},                    #  data of manually created nuclei, 'ß' are the (Z,X,Y) coordinates, 'Ł' are tilting parameters (x,y,gama)
                   4: {'ß':(0, 50, 30), 'Ł': (0,0,45*0.5)},
                   5: {'ß':(0, 32, 42), 'Ł': (0,0,45*0.75)},
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



if random_distances:
       R=Dissipate_Distances(epsilon)
else:
       R=1

fn =     1                                                               # Additional weight (factor) of orientation weight (W) for first neighbours (fn)
sn =    1                                                               # Additional weight (factor) of orientation weight (W) for second neighbours (sn)
en =    1  #1/math.sqrt(2)                                   # Additional weight (factor) of orientation weight (W) for extra neighbours (en)
en2 =  1
en3 =  1


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

# Plotting Data (PD)
PD = {'i': [], 'current step CPU': [], 'ALL grains': []}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ... STARTING the CA ..
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Growth (time FOR-LOOP, grains FOR-LOOP, directions FOR-LOOP)~~~~~~~~~~~~~~~~~~~~~~~~~
start_time = time.time()

for i in range(START_step, END_step+1):
    step_time_start = time.time(); h+=1
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

        '''............. the following code can be written with numpy broadcasting to avoid outrageously slow for-loops.. .....................................'''
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

        try:
            Selection = Selection_Mechanism(grain_ID_, smeri_database, pick_selection_mechanism)   

            #print('dolžina grain_ID_: ', len(grain_ID_))
            grain_ID = grain_ID_.copy()
        except NameError:
            grain_ID = np.array(list(asc.keys())); print('dolžina grain_ID: ', len(grain_ID))
            Selection = Selection_Mechanism(grain_ID, smeri_database, pick_selection_mechanism) 
        
        cas[np.isin(faza, new_grains_ID, invert=False)] = -1          # vrednost časovne matrike vseh novih nukleusov je -1
        vg = growth_speed(dTt)

    elif not avtomatska_nukleacija:
        Selection = Selection_Mechanism(grain_ID, smeri_database, pick_selection_mechanism)
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 

    # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   2.   T  I  M  E   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    if run_1st:
        for negind in Negatives:
            time_counter(negind)

    if run_2nd:
        for s in smeri[:]:
            for negind in Negatives:
                ''' ----------------------------------------- GROUP 1 ::: [001], [010] >>> 4 sites -------------------------------------'''
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
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 2 ::: [100] >>> 1 site -------------------------------------------'''
                elif s == '100':
                    if negind == -1:
                        cas100 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==negind), Negatives[negind], cas)
                        S[s]=cas100; del cas100
                    else:
                        cas100 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==negind), Negatives[negind], S[s])
                        S[s]=cas100; del cas100
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 3 ::: [_100] >>> 1 site -------------------------------------------'''
                elif s == '_100':
                    if negind == -1:
                        cas_100 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==negind), Negatives[negind], cas)
                        S[s]=cas_100; del cas_100
                    else:
                        cas_100 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==negind), Negatives[negind], S[s])
                        S[s]=cas_100; del cas_100
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 4 ::: [011] >>> 4 sites -------------------------------------------'''
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
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 5 ::: [101], [110] >>> 4 sites -------------------------------------------'''
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
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 6 ::: [_101], [_110] >>> 4 sites -------------------------------------------'''
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
                        
                        ''' ----------------------------------------- GROUP 7 ::: [111] >>> 4 sites -------------------------------------------'''
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
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 8 ::: [_111] >>> 4 sites -------------------------------------------'''
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
                        
                        ''' ----------------------------------------- GROUP 9 ::: [012], [021] >>> 8 sites -------------------------------------------'''
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
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 10 ::: [102], [120] >>> 4 sites -------------------------------------------'''
                elif s == '102':
                   if negind == -1:
                        cas102 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas102; del cas102
                   else:
                        cas102 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas102; del cas102

                elif s == '10_2':
                   if negind == -1:
                        cas10_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==negind), Negatives[negind], cas)
                        S[s]=cas10_2; del cas10_2
                   else:
                        cas10_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas10_2; del cas10_2

                elif s == '120':
                   if negind == -1:
                        cas120 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==negind), Negatives[negind], cas)
                        S[s]=cas120; del cas120
                   else:
                        cas120 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==negind), Negatives[negind], S[s])
                        S[s]=cas120; del cas120

                elif s == '1_20':
                   if negind == -1:
                        cas1_20 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==negind), Negatives[negind], cas)
                        S[s]=cas1_20; del cas1_20
                   else:
                        cas1_20 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==negind), Negatives[negind], S[s])
                        S[s]=cas1_20; del cas1_20
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 11 ::: [_102], [_120] >>> 4 sites -------------------------------------------'''
                elif s == '_102':
                   if negind == -1:
                        cas_102 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,0), (2,0)), 'constant')[1:,:,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas_102; del cas_102
                   else:
                        cas_102 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,0), (2,0)), 'constant')[1:,:,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas_102; del cas_102

                elif s == '_10_2':
                   if negind == -1:
                        cas_10_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,0), (0,2)), 'constant')[1:,:,2:]==negind), Negatives[negind], cas)
                        S[s]=cas_10_2; del cas_10_2
                   else:
                        cas_10_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,0), (0,2)), 'constant')[1:,:,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas_10_2; del cas_10_2

                elif s == '_120':
                   if negind == -1:
                        cas_120 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(2,0), (0,0)), 'constant')[1:,:-2,:]==negind), Negatives[negind], cas)
                        S[s]=cas_120; del cas_120
                   else:
                        cas_120 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(2,0), (0,0)), 'constant')[1:,:-2,:]==negind), Negatives[negind], S[s])
                        S[s]=cas_120; del cas_120

                elif s == '_1_20':
                   if negind == -1:
                        cas_1_20 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,2), (0,0)), 'constant')[1:,2:,:]==negind), Negatives[negind], cas)
                        S[s]=cas_1_20; del cas_1_20
                   else:
                        cas_1_20 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,2), (0,0)), 'constant')[1:,2:,:]==negind), Negatives[negind], S[s])
                        S[s]=cas_1_20; del cas_1_20
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 12 ::: [201], [210] >>> 4 sites -------------------------------------------'''
                elif s == '201':
                   if negind == -1:
                        cas201 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,0), (1,0)), 'constant')[:-2,:,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas201; del cas201
                   else:
                        cas201 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,0), (1,0)), 'constant')[:-2,:,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas201; del cas201

                elif s == '20_1':
                   if negind == -1:
                        cas20_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,0), (0,1)), 'constant')[:-2,:,1:]==negind), Negatives[negind], cas)
                        S[s]=cas20_1; del cas20_1
                   else:
                        cas20_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,0), (0,1)), 'constant')[:-2,:,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas20_1; del cas20_1

                elif s == '210':
                   if negind == -1:
                        cas210 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(1,0), (0,0)), 'constant')[:-2,:-1,:]==negind), Negatives[negind], cas)
                        S[s]=cas210; del cas210
                   else:
                        cas210 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(1,0), (0,0)), 'constant')[:-2,:-1,:]==negind), Negatives[negind], S[s])
                        S[s]=cas210; del cas210

                elif s == '2_10':
                   if negind == -1:
                        cas2_10 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,1), (0,0)), 'constant')[:-2,1:,:]==negind), Negatives[negind], cas)
                        S[s]=cas2_10; del cas2_10
                   else:
                        cas2_10 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,1), (0,0)), 'constant')[:-2,1:,:]==negind), Negatives[negind], S[s])
                        S[s]=cas2_10; del cas2_10
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 13 ::: [_201], [_210] >>> 4 sites -------------------------------------------'''
                elif s == '_201':
                   if negind == -1:
                        cas_201 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,0), (1,0)), 'constant')[2:,:,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas_201; del cas_201
                   else:
                        cas_201 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,0), (1,0)), 'constant')[2:,:,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas_201; del cas_201

                elif s == '_20_1':
                   if negind == -1:
                        cas_20_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,0), (0,1)), 'constant')[2:,:,1:]==negind), Negatives[negind], cas)
                        S[s]=cas_20_1; del cas_20_1
                   else:
                        cas_20_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,0), (0,1)), 'constant')[2:,:,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas_20_1; del cas_20_1

                elif s == '_210':
                   if negind == -1:
                        cas_210 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(1,0), (0,0)), 'constant')[2:,:-1,:]==negind), Negatives[negind], cas)
                        S[s]=cas_210; del cas_210
                   else:
                        cas_210 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(1,0), (0,0)), 'constant')[2:,:-1,:]==negind), Negatives[negind], S[s])
                        S[s]=cas_210; del cas_210

                elif s == '_2_10':
                   if negind == -1:
                        cas_2_10 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,1), (0,0)), 'constant')[2:,1:,:]==negind), Negatives[negind], cas)
                        S[s]=cas_2_10; del cas_2_10
                   else:
                        cas_2_10 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,1), (0,0)), 'constant')[2:,1:,:]==negind), Negatives[negind], S[s])
                        S[s]=cas_2_10; del cas_2_10
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 14 ::: [002], [020] >>> 4 sites -------------------------------------------'''
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
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 15 ::: [200] >>> 1 site -------------------------------------------'''
                elif s == '200':
                    if negind == -1:
                        cas200 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,0), (0,0)), 'constant')[:-2,:,:]==negind), Negatives[negind], cas)
                        S[s]=cas200; del cas200
                    else:
                        cas200 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,0), (0,0)), 'constant')[:-2,:,:]==negind), Negatives[negind], S[s])
                        S[s]=cas200; del cas200
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 16 ::: [_200] >>> 1 site -------------------------------------------'''
                elif s == '_200':
                    if negind == -1:
                        cas_200 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,0), (0,0)), 'constant')[2:,:,:]==negind), Negatives[negind], cas)
                        S[s]=cas_200; del cas_200
                    else:
                        cas_200 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,0), (0,0)), 'constant')[2:,:,:]==negind), Negatives[negind], S[s])
                        S[s]=cas_200; del cas_200
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 17 ::: [022] >>> 4 sites -------------------------------------------'''
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
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 18 ::: [202], [220] >>> 4 sites -------------------------------------------'''
                elif s == '202':
                   if negind == -1:
                        cas202 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,0), (2,0)), 'constant')[:-2,:,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas202; del cas202
                   else:
                        cas202 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,0), (2,0)), 'constant')[:-2,:,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas202; del cas202

                elif s == '20_2':
                   if negind == -1:
                        cas20_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,0), (0,2)), 'constant')[:-2,:,2:]==negind), Negatives[negind], cas)
                        S[s]=cas20_2; del cas20_2
                   else:
                        cas20_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,0), (0,2)), 'constant')[:-2,:,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas20_2; del cas20_2

                elif s == '220':
                   if negind == -1:
                        cas220 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(2,0), (0,0)), 'constant')[:-2,:-2,:]==negind), Negatives[negind], cas)
                        S[s]=cas220; del cas220
                   else:
                        cas220 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(2,0), (0,0)), 'constant')[:-2,:-2,:]==negind), Negatives[negind], S[s])
                        S[s]=cas220; del cas220

                elif s == '2_20':
                   if negind == -1:
                        cas2_20 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==negind), Negatives[negind], cas)
                        S[s]=cas2_20; del cas2_20
                   else:
                        cas2_20 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==negind), Negatives[negind], S[s])
                        S[s]=cas2_20; del cas2_20
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 19 ::: [_202], [_220] >>> 4 sites -------------------------------------------'''
                elif s == '_202':
                   if negind == -1:
                        cas_202 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,0), (2,0)), 'constant')[2:,:,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas_202; del cas_202
                   else:
                        cas_202 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,0), (2,0)), 'constant')[2:,:,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas_202; del cas_202

                elif s == '_20_2':
                   if negind == -1:
                        cas_20_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,0), (0,2)), 'constant')[2:,:,2:]==negind), Negatives[negind], cas)
                        S[s]=cas_20_2; del cas_20_2
                   else:
                        cas_20_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,0), (0,2)), 'constant')[2:,:,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas_20_2; del cas_20_2

                elif s == '_220':
                   if negind == -1:
                        cas_220 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(2,0), (0,0)), 'constant')[2:,:-2,:]==negind), Negatives[negind], cas)
                        S[s]=cas_220; del cas_220
                   else:
                        cas_220 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(2,0), (0,0)), 'constant')[2:,:-2,:]==negind), Negatives[negind], S[s])
                        S[s]=cas_220; del cas_220

                elif s == '_2_20':
                   if negind == -1:
                        cas_2_20 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,2), (0,0)), 'constant')[2:,2:,:]==negind), Negatives[negind], cas)
                        S[s]=cas_2_20; del cas_2_20
                   else:
                        cas_2_20 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,2), (0,0)), 'constant')[2:,2:,:]==negind), Negatives[negind], S[s])
                        S[s]=cas_2_20; del cas_2_20
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 20 ::: [112], [121] >>> 8 sites -------------------------------------------'''
                elif s == '112':
                   if negind == -1:
                        cas112 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas112; del cas112
                   else:
                        cas112 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas112; del cas112

                elif s == '11_2':
                   if negind == -1:
                        cas11_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==negind), Negatives[negind], cas)
                        S[s]=cas11_2; del cas11_2
                   else:
                        cas11_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas11_2; del cas11_2

                elif s == '1_12':
                   if negind == -1:
                        cas1_12 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas1_12; del cas1_12
                   else:
                        cas1_12 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas1_12; del cas1_12

                elif s == '1_1_2':
                   if negind == -1:
                        cas1_1_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==negind), Negatives[negind], cas)
                        S[s]=cas1_1_2; del cas1_1_2
                   else:
                        cas1_1_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas1_1_2; del cas1_1_2

                elif s == '121':
                   if negind == -1:
                        cas121 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas121; del cas121
                   else:
                        cas121 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas121; del cas121

                elif s == '12_1':
                   if negind == -1:
                        cas12_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==negind), Negatives[negind], cas)
                        S[s]=cas12_1; del cas12_1
                   else:
                        cas12_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas12_1; del cas12_1

                elif s == '1_21':
                   if negind == -1:
                        cas1_21 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas1_21; del cas1_21
                   else:
                        cas1_21 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas1_21; del cas1_21

                elif s == '1_2_1':
                   if negind == -1:
                        cas1_2_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==negind), Negatives[negind], cas)
                        S[s]=cas1_2_1; del cas1_2_1
                   else:
                        cas1_2_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas1_2_1; del cas1_2_1
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 21 ::: [_112], [_121] >>> 8 sites -------------------------------------------'''
                elif s == '_112':
                   if negind == -1:
                        cas_112 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(1,0), (2,0)), 'constant')[1:,:-1,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas_112; del cas_112
                   else:
                        cas_112 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(1,0), (2,0)), 'constant')[1:,:-1,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas_112; del cas_112

                elif s == '_11_2':
                   if negind == -1:
                        cas_11_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(1,0), (0,2)), 'constant')[1:,:-1,2:]==negind), Negatives[negind], cas)
                        S[s]=cas_11_2; del cas_11_2
                   else:
                        cas_11_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(1,0), (0,2)), 'constant')[1:,:-1,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas_11_2; del cas_11_2

                elif s == '_1_12':
                   if negind == -1:
                        cas_1_12 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,1), (2,0)), 'constant')[1:,1:,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas_1_12; del cas_1_12
                   else:
                        cas_1_12 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,1), (2,0)), 'constant')[1:,1:,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas_1_12; del cas_1_12

                elif s == '_1_1_2':
                   if negind == -1:
                        cas_1_1_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,1), (0,2)), 'constant')[1:,1:,2:]==negind), Negatives[negind], cas)
                        S[s]=cas_1_1_2; del cas_1_1_2
                   else:
                        cas_1_1_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,1), (0,2)), 'constant')[1:,1:,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas_1_1_2; del cas_1_1_2

                elif s == '_121':
                   if negind == -1:
                        cas_121 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(2,0), (1,0)), 'constant')[1:,:-2,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas_121; del cas_121
                   else:
                        cas_121 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(2,0), (1,0)), 'constant')[1:,:-2,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas_121; del cas_121

                elif s == '_12_1':
                   if negind == -1:
                        cas_12_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(2,0), (0,1)), 'constant')[1:,:-2,1:]==negind), Negatives[negind], cas)
                        S[s]=cas_12_1; del cas_12_1
                   else:
                        cas_12_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(2,0), (0,1)), 'constant')[1:,:-2,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas_12_1; del cas_12_1

                elif s == '_1_21':
                   if negind == -1:
                        cas_1_21 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,2), (1,0)), 'constant')[1:,2:,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas_1_21; del cas_1_21
                   else:
                        cas_1_21 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,2), (1,0)), 'constant')[1:,2:,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas_1_21; del cas_1_21

                elif s == '_1_2_1':
                   if negind == -1:
                        cas_1_2_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,2), (0,1)), 'constant')[1:,2:,1:]==negind), Negatives[negind], cas)
                        S[s]=cas_1_2_1; del cas_1_2_1
                   else:
                        cas_1_2_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,2), (0,1)), 'constant')[1:,2:,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas_1_2_1; del cas_1_2_1
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 22 ::: [211] >>> 4 sites -------------------------------------------'''
                elif s == '211':
                   if negind == -1:
                        cas211 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(1,0), (1,0)), 'constant')[:-2,:-1,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas211; del cas211
                   else:
                        cas211 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(1,0), (1,0)), 'constant')[:-2,:-1,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas211; del cas211

                elif s == '21_1':
                   if negind == -1:
                        cas21_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(1,0), (0,1)), 'constant')[:-2,:-1,1:]==negind), Negatives[negind], cas)
                        S[s]=cas21_1; del cas21_1
                   else:
                        cas21_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(1,0), (0,1)), 'constant')[:-2,:-1,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas21_1; del cas21_1

                elif s == '2_11':
                   if negind == -1:
                        cas2_11 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,1), (1,0)), 'constant')[:-2,1:,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas2_11; del cas2_11
                   else:
                        cas2_11 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,1), (1,0)), 'constant')[:-2,1:,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas2_11; del cas2_11

                elif s == '2_1_1':
                   if negind == -1:
                        cas2_1_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,1), (0,1)), 'constant')[:-2,1:,1:]==negind), Negatives[negind], cas)
                        S[s]=cas2_1_1; del cas2_1_1
                   else:
                        cas2_1_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,1), (0,1)), 'constant')[:-2,1:,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas2_1_1; del cas2_1_1
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 23 ::: [_211] >>> 4 sites -------------------------------------------'''
                elif s == '_211':
                   if negind == -1:
                        cas_211 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(1,0), (1,0)), 'constant')[2:,:-1,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas_211; del cas_211
                   else:
                        cas_211 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(1,0), (1,0)), 'constant')[2:,:-1,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas_211; del cas_211

                elif s == '_21_1':
                   if negind == -1:
                        cas_21_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(1,0), (0,1)), 'constant')[2:,:-1,1:]==negind), Negatives[negind], cas)
                        S[s]=cas_21_1; del cas_21_1
                   else:
                        cas_21_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(1,0), (0,1)), 'constant')[2:,:-1,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas_21_1; del cas_21_1

                elif s == '_2_11':
                   if negind == -1:
                        cas_2_11 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,1), (1,0)), 'constant')[2:,1:,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas_2_11; del cas_2_11
                   else:
                        cas_2_11 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,1), (1,0)), 'constant')[2:,1:,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas_2_11; del cas_2_11

                elif s == '_2_1_1':
                   if negind == -1:
                        cas_2_1_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,1), (0,1)), 'constant')[2:,1:,1:]==negind), Negatives[negind], cas)
                        S[s]=cas_2_1_1; del cas_2_1_1
                   else:
                        cas_2_1_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,1), (0,1)), 'constant')[2:,1:,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas_2_1_1; del cas_2_1_1
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 24 ::: [122] >>> 4 sites -------------------------------------------'''
                elif s == '122':
                   if negind == -1:
                        cas122 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas122; del cas122
                   else:
                        cas122 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas122; del cas122

                elif s == '12_2':
                   if negind == -1:
                        cas12_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==negind), Negatives[negind], cas)
                        S[s]=cas12_2; del cas12_2
                   else:
                        cas12_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas12_2; del cas12_2

                elif s == '1_22':
                   if negind == -1:
                        cas1_22 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas1_22; del cas1_22
                   else:
                        cas1_22 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas1_22; del cas1_22

                elif s == '1_2_2':
                   if negind == -1:
                        cas1_2_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==negind), Negatives[negind], cas)
                        S[s]=cas1_2_2; del cas1_2_2
                   else:
                        cas1_2_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas1_2_2; del cas1_2_2
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 25 ::: [_122] >>> 4 sites -------------------------------------------'''
                elif s == '_122':
                   if negind == -1:
                        cas_122 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(2,0), (2,0)), 'constant')[1:,:-2,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas_122; del cas_122
                   else:
                        cas_122 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(2,0), (2,0)), 'constant')[1:,:-2,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas_122; del cas_122

                elif s == '_12_2':
                   if negind == -1:
                        cas_12_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(2,0), (0,2)), 'constant')[1:,:-2,2:]==negind), Negatives[negind], cas)
                        S[s]=cas_12_2; del cas_12_2
                   else:
                        cas_12_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(2,0), (0,2)), 'constant')[1:,:-2,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas_12_2; del cas_12_2

                elif s == '_1_22':
                   if negind == -1:
                        cas_1_22 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,2), (2,0)), 'constant')[1:,2:,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas_1_22; del cas_1_22
                   else:
                        cas_1_22 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,2), (2,0)), 'constant')[1:,2:,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas_1_22; del cas_1_22

                elif s == '_1_2_2':
                   if negind == -1:
                        cas_1_2_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,1),(0,2), (0,2)), 'constant')[1:,2:,2:]==negind), Negatives[negind], cas)
                        S[s]=cas_1_2_2; del cas_1_2_2
                   else:
                        cas_1_2_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,1),(0,2), (0,2)), 'constant')[1:,2:,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas_1_2_2; del cas_1_2_2
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 26 ::: [212], [221] >>> 8 sites -------------------------------------------'''
                elif s == '212':
                   if negind == -1:
                        cas212 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(1,0), (2,0)), 'constant')[:-2,:-1,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas212; del cas212
                   else:
                        cas212 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(1,0), (2,0)), 'constant')[:-2,:-1,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas212; del cas212

                elif s == '21_2':
                   if negind == -1:
                        cas21_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(1,0), (0,2)), 'constant')[:-2,:-1,2:]==negind), Negatives[negind], cas)
                        S[s]=cas21_2; del cas21_2
                   else:
                        cas21_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(1,0), (0,2)), 'constant')[:-2,:-1,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas21_2; del cas21_2

                elif s == '2_12':
                   if negind == -1:
                        cas2_12 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,1), (2,0)), 'constant')[:-2,1:,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas2_12; del cas2_12
                   else:
                        cas2_12 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,1), (2,0)), 'constant')[:-2,1:,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas2_12; del cas2_12

                elif s == '2_1_2':
                   if negind == -1:
                        cas2_1_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,1), (0,2)), 'constant')[:-2,1:,2:]==negind), Negatives[negind], cas)
                        S[s]=cas2_1_2; del cas2_1_2
                   else:
                        cas2_1_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,1), (0,2)), 'constant')[:-2,1:,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas2_1_2; del cas2_1_2

                elif s == '221':
                   if negind == -1:
                        cas221 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(2,0), (1,0)), 'constant')[:-2,:-2,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas221; del cas221
                   else:
                        cas221 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(2,0), (1,0)), 'constant')[:-2,:-2,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas221; del cas221

                elif s == '22_1':
                   if negind == -1:
                        cas22_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(2,0), (0,1)), 'constant')[:-2,:-2,1:]==negind), Negatives[negind], cas)
                        S[s]=cas22_1; del cas22_1
                   else:
                        cas22_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(2,0), (0,1)), 'constant')[:-2,:-2,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas22_1; del cas22_1

                elif s == '2_21':
                   if negind == -1:
                        cas2_21 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,2), (1,0)), 'constant')[:-2,2:,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas2_21; del cas2_21
                   else:
                        cas2_21 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,2), (1,0)), 'constant')[:-2,2:,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas2_21; del cas2_21

                elif s == '2_2_1':
                   if negind == -1:
                        cas2_2_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==negind), Negatives[negind], cas)
                        S[s]=cas2_2_1; del cas2_2_1
                   else:
                        cas2_2_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas2_2_1; del cas2_2_1
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 27 ::: [_212], [_221] >>> 8 sites -------------------------------------------'''
                elif s == '_212':
                   if negind == -1:
                        cas_212 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(1,0), (2,0)), 'constant')[2:,:-1,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas_212; del cas_212
                   else:
                        cas_212 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(1,0), (2,0)), 'constant')[2:,:-1,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas_212; del cas_212

                elif s == '_21_2':
                   if negind == -1:
                        cas_21_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(1,0), (0,2)), 'constant')[2:,:-1,2:]==negind), Negatives[negind], cas)
                        S[s]=cas_21_2; del cas_21_2
                   else:
                        cas_21_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(1,0), (0,2)), 'constant')[2:,:-1,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas_21_2; del cas_21_2

                elif s == '_2_12':
                   if negind == -1:
                        cas_2_12 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,1), (2,0)), 'constant')[2:,1:,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas_2_12; del cas_2_12
                   else:
                        cas_2_12 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,1), (2,0)), 'constant')[2:,1:,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas_2_12; del cas_2_12

                elif s == '_2_1_2':
                   if negind == -1:
                        cas_2_1_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,1), (0,2)), 'constant')[2:,1:,2:]==negind), Negatives[negind], cas)
                        S[s]=cas_2_1_2; del cas_2_1_2
                   else:
                        cas_2_1_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,1), (0,2)), 'constant')[2:,1:,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas_2_1_2; del cas_2_1_2

                elif s == '_221':
                   if negind == -1:
                        cas_221 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(2,0), (1,0)), 'constant')[2:,:-2,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas_221; del cas_221
                   else:
                        cas_221 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(2,0), (1,0)), 'constant')[2:,:-2,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas_221; del cas_221

                elif s == '_22_1':
                   if negind == -1:
                        cas_22_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(2,0), (0,1)), 'constant')[2:,:-2,1:]==negind), Negatives[negind], cas)
                        S[s]=cas_22_1; del cas_22_1
                   else:
                        cas_22_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(2,0), (0,1)), 'constant')[2:,:-2,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas_22_1; del cas_22_1

                elif s == '_2_21':
                   if negind == -1:
                        cas_2_21 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,2), (1,0)), 'constant')[2:,2:,:-1]==negind), Negatives[negind], cas)
                        S[s]=cas_2_21; del cas_2_21
                   else:
                        cas_2_21 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,2), (1,0)), 'constant')[2:,2:,:-1]==negind), Negatives[negind], S[s])
                        S[s]=cas_2_21; del cas_2_21

                elif s == '_2_2_1':
                   if negind == -1:
                        cas_2_2_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,2), (0,1)), 'constant')[2:,2:,1:]==negind), Negatives[negind], cas)
                        S[s]=cas_2_2_1; del cas_2_2_1
                   else:
                        cas_2_2_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,2), (0,1)), 'constant')[2:,2:,1:]==negind), Negatives[negind], S[s])
                        S[s]=cas_2_2_1; del cas_2_2_1
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
                        ''' ----------------------------------------- GROUP 28 ::: [222] >>> 4 sites -------------------------------------------'''
                elif s == '222':
                   if negind == -1:
                        cas222 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(2,0), (2,0)), 'constant')[:-2,:-2,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas222; del cas222
                   else:
                        cas222 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(2,0), (2,0)), 'constant')[:-2,:-2,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas222; del cas222

                elif s == '22_2':
                   if negind == -1:
                        cas22_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(2,0), (0,2)), 'constant')[:-2,:-2,2:]==negind), Negatives[negind], cas)
                        S[s]=cas22_2; del cas22_2
                   else:
                        cas22_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(2,0), (0,2)), 'constant')[:-2,:-2,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas22_2; del cas22_2

                elif s == '2_22':
                   if negind == -1:
                        cas2_22 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,2), (2,0)), 'constant')[:-2,2:,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas2_22; del cas2_22
                   else:
                        cas2_22 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,2), (2,0)), 'constant')[:-2,2:,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas2_22; del cas2_22

                elif s == '2_2_2':
                   if negind == -1:
                        cas2_2_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==negind), Negatives[negind], cas)
                        S[s]=cas2_2_2; del cas2_2_2
                   else:
                        cas2_2_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas2_2_2; del cas2_2_2
                        ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        # under construction
                        ''' ----------------------------------------- GROUP 29 ::: [_222] >>> 4 sites -------------------------------------------'''
                elif s == '_222':
                   if negind == -1:
                        cas_222 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(2,0), (2,0)), 'constant')[2:,:-2,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas_222; del cas_222
                   else:
                        cas_222 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(2,0), (2,0)), 'constant')[2:,:-2,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas_222; del cas_222

                elif s == '_22_2':
                   if negind == -1:
                        cas_22_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(2,0), (0,2)), 'constant')[2:,:-2,2:]==negind), Negatives[negind], cas)
                        S[s]=cas_22_2; del cas_22_2
                   else:
                        cas_22_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(2,0), (0,2)), 'constant')[2:,:-2,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas_22_2; del cas_22_2

                elif s == '_2_22':
                   if negind == -1:
                        cas_2_22 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,2), (2,0)), 'constant')[2:,2:,:-2]==negind), Negatives[negind], cas)
                        S[s]=cas_2_22; del cas_2_22
                   else:
                        cas_2_22 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,2), (2,0)), 'constant')[2:,2:,:-2]==negind), Negatives[negind], S[s])
                        S[s]=cas_2_22; del cas_2_22

                elif s == '_2_2_2':
                   if negind == -1:
                        cas_2_2_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,2),(0,2), (0,2)), 'constant')[2:,2:,2:]==negind), Negatives[negind], cas)
                        S[s]=cas_2_2_2; del cas_2_2_2
                   else:
                        cas_2_2_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,2),(0,2), (0,2)), 'constant')[2:,2:,2:]==negind), Negatives[negind], S[s])
                        S[s]=cas_2_2_2; del cas_2_2_2

    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 

    # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   3.   P  H  A  S   E   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    if run_3rd:
        F = faza.copy()
        FF.append(faza)
        if GDSM == 'old':
           seznam_premikov=[]               # old
        for choice in Selection:
           if GDSM == 'new':
               seznam_premikov=[]           # NEW
           grain=choice[0]
           s=choice[1]
        
           ''' ----------------------------------------- GROUP 1 ::: [001], [010] >>> 4 sites -------------------------------------------'''
           if s == '001':
               dij = np.array([0,0,1])
               wij = fn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza001 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain) , grain, faza)
               seznam_premikov.append(faza001); del faza001

           elif s == '00_1':
               dij = np.array([0,0,-1])
               wij = fn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza00_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)& (lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain), grain, faza)
               seznam_premikov.append(faza00_1); del faza00_1

           elif s == '010':
               dij = np.array([0,1,0])
               wij = fn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza010 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain), grain, faza)
               seznam_premikov.append(faza010); del faza010

           elif s == '0_10':
               dij = np.array([0,-1,0])
               wij = fn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza0_10 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain), grain, faza)
               seznam_premikov.append(faza0_10); del faza0_10
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
               ''' ----------------------------------------- GROUP 2 ::: [100] >>> 1 site -------------------------------------------'''
           elif s == '100':
               dij = np.array([1,0,0])
               wij = fn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza100 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)& (lij>=Dsr_1st(R) )&(np.pad(faza,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain), grain, faza)
               seznam_premikov.append(faza100); del faza100
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
               ''' ----------------------------------------- GROUP 3 ::: [_100] >>> 1 site -------------------------------------------'''
           elif s == '_100':
               dij = np.array([-1,0,0])
               wij = fn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_100 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_1st(R) )&(np.pad(faza,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain), grain, faza)
               seznam_premikov.append(faza_100); del faza_100
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 4 ::: [011] >>> 4 sites -------------------------------------------'''
           elif s == '011':
               dij = np.array([0,1,1])
               wij = W(asc[grain]['oi'], dij)
               lij = sn*wij*vg*S[s]
               faza011 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R) )&(np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)& (
                                                                      #(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)         |              # 001
                                                                      #(np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)         |              # 010
                                                                      #(np.pad(faza,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)   |              # 111
                                                                      #(np.pad(faza,((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==grain)                   # _111
                                                True), grain, faza)
               seznam_premikov.append(faza011); del faza011

           elif s == '01_1':
               dij = np.array([0,1,-1])
               wij = sn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza01_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R) )&(np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)& (
                                                                         #(np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)            |          # 00_1
                                                                         #(np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)           |          # 010
                                                                         #(np.pad(faza,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)      |          # 11_1
                                                                         #(np.pad(faza,((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==grain)                  # _11_1
                                                True), grain, faza)
               seznam_premikov.append(faza01_1); del faza01_1

           elif s == '0_11':
               dij = np.array([0,-1,1])
               wij = sn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza0_11 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R) )&(np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)& (
                                                                         #(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)           |          # 001
                                                                         #(np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)           |          # 0_10
                                                                         #(np.pad(faza,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)     |           # 1_11
                                                                         #(np.pad(faza,((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==grain)                  # _1_11
                                                True), grain, faza)
               seznam_premikov.append(faza0_11); del faza0_11

           elif s == '0_1_1':
               dij = np.array([0,-1,-1])
               wij = sn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza0_1_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R) )&(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)& (
                                                                            #(np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)         |         # 00_1
                                                                            #(np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)       |         # 0_10
                                                                            #(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)    |         # 1_1_1
                                                                            #(np.pad(faza,((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==grain)               # _1_1_1
                                                True), grain, faza)
               seznam_premikov.append(faza0_1_1); del faza0_1_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 5 ::: [101], [110] >>> 4 sites -------------------------------------------'''
           elif s == '101':
               dij = np.array([1,0,1])
               wij = sn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza101 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==grain)& (
                                                                       #(np.pad(faza,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain)             |        # 100
                                                                       #(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)             |        # 001          
                                                                       #(np.pad(faza,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)       |        # 111      
                                                                       #(np.pad(faza,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)                 # 1_11         
                                                True), grain, faza)
               seznam_premikov.append(faza101); del faza101

           elif s == '10_1':
               dij = np.array([1,0,-1])
               wij = sn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza10_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==grain)& (
                                                                         #(np.pad(faza,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain)            |       # 100
                                                                         #(np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)             |       # 00_1          
                                                                         #(np.pad(faza,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)       |       # 11_1      
                                                                         #(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)                # 1_1_1
                                               True), grain, faza)
               seznam_premikov.append(faza10_1); del faza10_1

           elif s == '110':
               dij = np.array([1,1,0])
               wij = sn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza110 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==grain)& (
                                                                      #(np.pad(faza,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain)             |         # 100
                                                                      #(np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)             |         # 010          
                                                                      #(np.pad(faza,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)       |           # 111      
                                                                      #(np.pad(faza,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)                    # 11_1
                                               True), grain, faza)
               seznam_premikov.append(faza110); del faza110

           elif s == '1_10':
               dij = np.array([1,-1,0])
               wij = sn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza1_10 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)& (
                                                                         #(np.pad(faza,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain)          |      # 100
                                                                         #(np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)           |      # 0_10          
                                                                         #(np.pad(faza,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)     |       # 1_11      
                                                                         #(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)              # 1_1_1
                                                True), grain, faza)
               seznam_premikov.append(faza1_10); del faza1_10
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 6 ::: [_101], [_110] >>> 4 sites -------------------------------------------'''
           elif s == '_101':
               dij = np.array([-1,0,1])
               wij = sn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_101 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((0,1),(0,0), (1,0)), 'constant')[1:,:,:-1]==grain)& (
                                                                         #(np.pad(faza,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain)          |       # _100
                                                                         #(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)         |       # 001          
                                                                         #(np.pad(faza,((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==grain)    |        # _111      
                                                                         #(np.pad(faza,((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==grain)             # _1_11
                                                True), grain, faza)
               seznam_premikov.append(faza_101); del faza_101

           elif s == '_10_1':
               dij = np.array([-1,0,-1])
               wij = sn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_10_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((0,1),(0,0), (0,1)), 'constant')[1:,:,1:]==grain)& (
                                                                            #(np.pad(faza,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain)           |       # _100
                                                                            #(np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)           |       # 00_1          
                                                                            #(np.pad(faza,((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==grain)       |        # _11_1      
                                                                            #(np.pad(faza,((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==grain)                # _1_1_1
                                                True), grain, faza)
               seznam_premikov.append(faza_10_1); del faza_10_1

           elif s == '_110':
               dij = np.array([-1,1,0])
               wij = sn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_110 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((0,1),(1,0), (0,0)), 'constant')[1:,:-1,:]==grain)& (
                                                                         #(np.pad(faza,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain)         |         # _100
                                                                         #(np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)        |         # 010          
                                                                         #(np.pad(faza,((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==grain)    |        # _111      
                                                                         #(np.pad(faza,((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==grain)              # _11_1
                                                True), grain, faza)
               seznam_premikov.append(faza_110); del faza_110

           elif s == '_1_10':
               dij = np.array([-1,-1,0])
               wij = sn*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_1_10 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_2nd(R))&(np.pad(faza,((0,1),(0,1), (0,0)), 'constant')[1:,1:,:]==grain)& (
                                                                            #(np.pad(faza,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain)       |      # _100
                                                                            #(np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)         |      # 0_10             
                                                                            #(np.pad(faza,((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==grain)     |      # _1_11      
                                                                            #(np.pad(faza,((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==grain)             # _1_1_1
                                                True), grain, faza)
               seznam_premikov.append(faza_1_10); del faza_1_10
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 7 ::: [111] >>> 4 sites -------------------------------------------'''
           elif s == '111':
               dij = np.array([1,1,1])
               wij = W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza111 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)& (
                                                                      #(np.pad(faza,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==grain)    |        # 101
                                                                      #(np.pad(faza,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==grain)    |        # 110          
                                                                      #(np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)             # 011      
                                                True), grain, faza)
               seznam_premikov.append(faza111); del faza111

           elif s == '11_1':
               dij = np.array([1,1,-1])
               wij = W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza11_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)& (
                                                                         #(np.pad(faza,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==grain)     |       # 10_1
                                                                         #(np.pad(faza,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==grain)    |       # 110          
                                                                         #(np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)             # 01_1
                                                True), grain, faza)
               seznam_premikov.append(faza11_1); del faza11_1

           elif s == '1_11':
               dij = np.array([1,-1,1])
               wij = W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza1_11 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)& (
                                                                         #(np.pad(faza,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==grain)    |       # 101
                                                                         #(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)     |       # 1_10          
                                                                         #(np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)             # 0_11
                                                True), grain, faza)
               seznam_premikov.append(faza1_11); del faza1_11

           elif s == '1_1_1':
               dij = np.array([1,-1,-1])
               wij = W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza1_1_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)& (
                                                                            #(np.pad(faza,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==grain)     |        # 10_1
                                                                            #(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)     |        # 1_10          
                                                                            #(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)               # 0_1_1
                                                True), grain, faza)
               seznam_premikov.append(faza1_1_1); del faza1_1_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 8 ::: [_111] >>> 4 sites -------------------------------------------'''
           elif s == '_111':
               dij = np.array([-1,1,1])
               wij = W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_111 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==grain)& (
                                                                         #(np.pad(faza,((0,1),(0,0), (1,0)), 'constant')[1:,:,:-1]==grain)     |       # _101
                                                                         #(np.pad(faza,((0,1),(1,0), (0,0)), 'constant')[1:,:-1,:]==grain)     |       # _110          
                                                                         #(np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)            # 011
                                                True), grain, faza)
               seznam_premikov.append(faza_111); del faza_111

           elif s == '_11_1':
               dij = np.array([-1,1,-1])
               wij = W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_11_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==grain)& (
                                                                            #(np.pad(faza,((0,1),(0,0), (0,1)), 'constant')[1:,:,1:]==grain)        |      # _10_1
                                                                            #(np.pad(faza,((0,1),(1,0), (0,0)), 'constant')[1:,:-1,:]==grain)       |      # _110          
                                                                            #(np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)              # 01_1
                                                True), grain, faza)
               seznam_premikov.append(faza_11_1); del faza_11_1

           elif s == '_1_11':
               dij = np.array([-1,-1,1])
               wij = W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_1_11 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==grain)& (
                                                                            #(np.pad(faza,((0,1),(0,0), (1,0)), 'constant')[1:,:,:-1]==grain)       |      # _101
                                                                            #(np.pad(faza,((0,1),(0,1), (0,0)), 'constant')[1:,1:,:]==grain)        |      # _1_10          
                                                                            #(np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)              # 0_11                
                                                True), grain, faza)
               seznam_premikov.append(faza_1_11); del faza_1_11

           elif s == '_1_1_1':
               dij = np.array([-1,-1,-1])
               wij = W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_1_1_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_3rd(R))&(np.pad(faza,((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==grain)& (
                                                                               #(np.pad(faza,((0,1),(0,0), (0,1)), 'constant')[1:,:,1:]==grain)    |         # _10_1
                                                                               #(np.pad(faza,((0,1),(0,1), (0,0)), 'constant')[1:,1:,:]==grain)    |        # _1_10          
                                                                               #(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)             # 0_1_1
                                                 True), grain, faza)
               seznam_premikov.append(faza_1_1_1); del faza_1_1_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 9 ::: [012], [021] >>> 8 sites -------------------------------------------'''
           elif s == '012':
               dij = np.array([0,1,2])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_9:
                   CG9=(#(np.pad(faza,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)       |      # 002
                              #(np.pad(faza,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)    |      # 021
                              #(np.pad(faza,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==grain)    |      # 022
                              (np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)       |      # 001
                              (np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)     )   # 011
               else:
                   CG9 = True
               faza012 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)& CG9, grain, faza) 
               seznam_premikov.append(faza012); del faza012

           elif s == '01_2':
               dij = np.array([0,1,-2])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(faza,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)         |         # 00_2
                                #(np.pad(faza,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)       |         # 02_1
                                #(np.pad(faza,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==grain)       |         # 02_2
                                 (np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)          |        # 00_1
                                 (np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)       )      # 01_1      
               else:
                   CG9 = True
               faza01_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)& CG9, grain, faza) 
               seznam_premikov.append(faza01_2); del faza01_2

           elif s == '0_12':
               dij = np.array([0,-1,2])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(faza,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)      |         # 002
                                #(np.pad(faza,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)     |         # 0_21
                                #(np.pad(faza,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==grain)     |       # 0_22
                                (np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)       |         # 001
                                (np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)     )       # 0_11         
               else:
                   CG9 = True
               faza0_12 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)& CG9, grain, faza)
               seznam_premikov.append(faza0_12); del faza0_12

           elif s == '0_1_2':
               dij = np.array([0,-1,-2])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(faza,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)       |        # 00_2
                                #(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)      |        # 0_2_1
                                #(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)      |        # 0_2_2
                                (np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)        |        # 00_1
                                (np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)       )     # 0_1_1         
               else:
                   CG9 = True
               faza0_1_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)& CG9, grain, faza) 
               seznam_premikov.append(faza0_1_2); del faza0_1_2

           elif s == '021':
               dij = np.array([0,2,1])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(faza,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)        |         # 020
                                #(np.pad(faza,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==grain)      |         # 022
                                #(np.pad(faza,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)      |        # 012
                                (np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)         |         # 010
                                (np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)      )       # 011         
               else:
                   CG9 = True
               faza021 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)& CG9, grain, faza)
               seznam_premikov.append(faza021); del faza021

           elif s == '02_1':
               dij = np.array([0,2,-1])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(faza,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)       |           # 020
                                #(np.pad(faza,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==grain)      |           # 02_2
                                #(np.pad(faza,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)      |           # 01_2
                                (np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)         |          # 010
                                (np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)        )       # 01_1
               else:
                   CG9 = True
               faza02_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)& CG9, grain, faza) 
               seznam_premikov.append(faza02_1); del faza02_1

           elif s == '0_21':
               dij = np.array([0,-2,1])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)        |          # 0_20
                                #(np.pad(faza,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==grain)      |         # 0_22
                                #(np.pad(faza,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)    |          # 0_12
                                (np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)         |         # 0_10
                                (np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)      )       # 0_11       
               else:
                   CG9 = True
               faza0_21 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)& CG9, grain, faza)
               seznam_premikov.append(faza0_21); del faza0_21

           elif s == '0_2_1':
               dij = np.array([0,-2,-1])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)     |       # 0_20
                                #(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)    |      # 0_2_2
                                #(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |       # 0_1_2
                                (np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)      |       # 0_10
                                (np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    )     # 0_1_1             
               else:
                   CG9 = True
               faza0_2_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)& CG9, grain, faza) 
               seznam_premikov.append(faza0_2_1); del faza0_2_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 10 ::: [102], [120] >>> 4 sites -------------------------------------------'''
           elif s == '102':
               dij = np.array([1,0,2])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza102 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==grain)& (
                                                                       #(np.pad(faza,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==grain)       |       # 101
                                                                       #(np.pad(faza,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)          |       # 002
                                                                       #(np.pad(faza,((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==grain)    |       # 112
                                                                       #(np.pad(faza,((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==grain)     |       # 1_12
                                                                       #(np.pad(faza,((2,0),(0,0), (2,0)), 'constant')[:-2,:,:-2]==grain)               # 202
                                                True), grain, faza) 
               seznam_premikov.append(faza102); del faza102

           elif s == '10_2':
               dij = np.array([1,0,-2])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza10_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==grain)& (
                                                                       #(np.pad(faza,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==grain)       |       # 10_1
                                                                       #(np.pad(faza,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)          |       # 00_2
                                                                       #(np.pad(faza,((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==grain)    |       # 11_2
                                                                       #(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)     |       # 1_1_2
                                                                       #(np.pad(faza,((2,0),(0,0), (0,2)), 'constant')[:-2,:,2:]==grain)               # 20_2
                                                True), grain, faza) 
               seznam_premikov.append(faza10_2); del faza10_2

           elif s == '120':
               dij = np.array([1,2,0])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza120 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==grain)& (
                                                                       #(np.pad(faza,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==grain)       |       # 110
                                                                       #(np.pad(faza,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)          |       # 020
                                                                       #(np.pad(faza,((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==grain)    |       # 121
                                                                       #(np.pad(faza,((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==grain)     |       # 12_1
                                                                       #(np.pad(faza,((2,0),(2,0), (0,0)), 'constant')[:-2,:-2,:]==grain)               # 220
                                                True), grain, faza) 
               seznam_premikov.append(faza120); del faza120

           elif s == '1_20':
               dij = np.array([1,-2,0])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza1_20 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)& (
                                                                       #(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       #(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       #(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       #(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       #(np.pad(faza,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, faza) 
               seznam_premikov.append(faza1_20); del faza1_20
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
        
               ''' ----------------------------------------- GROUP 11 ::: [_102], [_120] >>> 4 sites -------------------------------------------'''
           elif s == '_102':
               dij = np.array([-1,0,2])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_102 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,1),(0,0), (2,0)), 'constant')[1:,:,:-2]==grain)& (
                                                                       #wrong(np.pad(faza,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==grain)       |       # 101
                                                                       #wrong(np.pad(faza,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)          |       # 002
                                                                       #wrong(np.pad(faza,((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==grain)    |       # 112
                                                                       #wrong(np.pad(faza,((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==grain)     |       # 1_12
                                                                       #wrong(np.pad(faza,((2,0),(0,0), (2,0)), 'constant')[:-2,:,:-2]==grain)               # 202
                                                True), grain, faza) 
               seznam_premikov.append(faza_102); del faza_102

           elif s == '_10_2':
               dij = np.array([-1,0,-2])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_10_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,1),(0,0), (0,2)), 'constant')[1:,:,2:]==grain)& (
                                                                       #wrong(np.pad(faza,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==grain)       |       # 10_1
                                                                       #wrong(np.pad(faza,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)          |       # 00_2
                                                                       #wrong(np.pad(faza,((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==grain)    |       # 11_2
                                                                       #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)     |       # 1_1_2
                                                                       #wrong(np.pad(faza,((2,0),(0,0), (0,2)), 'constant')[:-2,:,2:]==grain)               # 20_2
                                               True), grain, faza) 
               seznam_premikov.append(faza_10_2); del faza_10_2

           elif s == '_120':
               dij = np.array([-1,2,0])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_120 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,1),(2,0), (0,0)), 'constant')[1:,:-2,:]==grain)& (
                                                                       #(np.pad(faza,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==grain)       |       # 110
                                                                       #(np.pad(faza,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)          |       # 020
                                                                       #(np.pad(faza,((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==grain)    |       # 121
                                                                       #(np.pad(faza,((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==grain)     |       # 12_1
                                                                       #(np.pad(faza,((2,0),(2,0), (0,0)), 'constant')[:-2,:-2,:]==grain)               # 220
                                                True), grain, faza) 
               seznam_premikov.append(faza_120); del faza_120

           elif s == '_1_20':
               dij = np.array([-1,-2,0])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_1_20 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,1),(0,2), (0,0)), 'constant')[1:,2:,:]==grain)& (
                                                                       #(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       #(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       #(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       #(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       #(np.pad(faza,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, faza) 
               seznam_premikov.append(faza_1_20); del faza_1_20
               

               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 12 ::: [201], [210] >>> 4 sites -------------------------------------------'''
           elif s == '201':
               dij = np.array([2,0,1])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza201 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((2,0),(0,0), (1,0)), 'constant')[:-2,:,:-1]==grain)& (
                                                                       # wrong(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(faza,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, faza) 
               seznam_premikov.append(faza201); del faza201

           elif s == '20_1':
               dij = np.array([2,0,-1])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza20_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((2,0),(0,0), (0,1)), 'constant')[:-2,:,1:]==grain)& (
                                                                       # wrong(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(faza,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, faza) 
               seznam_premikov.append(faza20_1); del faza20_1

           elif s == '210':
               dij = np.array([2,1,0])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza210 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((2,0),(1,0), (0,0)), 'constant')[:-2,:-1,:]==grain)& (
                                                                       # wrong(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(faza,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, faza) 
               seznam_premikov.append(faza210); del faza210

           elif s == '2_10':
               dij = np.array([2,-1,0])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza2_10 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((2,0),(0,1), (0,0)), 'constant')[:-2,1:,:]==grain)& (
                                                                       # wrong(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(faza,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, faza) 
               seznam_premikov.append(faza2_10); del faza2_10

               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 13 ::: [_201], [_210] >>> 4 sites -------------------------------------------'''
           elif s == '_201':
               dij = np.array([-2,0,1])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_201 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,2),(0,0), (1,0)), 'constant')[2:,:,:-1]==grain)& (
                                                                       # wrong(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(faza,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, faza) 
               seznam_premikov.append(faza_201); del faza_201

           elif s == '_20_1':
               dij = np.array([-2,0,-1])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_20_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,2),(0,0), (0,1)), 'constant')[2:,:,1:]==grain)& (
                                                                       # wrong(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(faza,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, faza) 
               seznam_premikov.append(faza_20_1); del faza_20_1

           elif s == '_210':
               dij = np.array([-2,1,0])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_210 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,2),(1,0), (0,0)), 'constant')[2:,:-1,:]==grain)& (
                                                                       # wrong(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(faza,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, faza) 
               seznam_premikov.append(faza_210); del faza_210

           elif s == '_2_10':
               dij = np.array([-2,-1,0])
               wij = en*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_2_10 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )&(np.pad(faza,((0,2),(0,1), (0,0)), 'constant')[2:,1:,:]==grain)& (
                                                                       # wrong(np.pad(faza,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(faza,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, faza) 
               seznam_premikov.append(faza_2_10); del faza_2_10
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 14 ::: [002], [020] >>> 4 sites -------------------------------------------'''
           elif s == '002':
               dij = np.array([0,0,2])
               wij = en2*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_14:
                   CG14 = (#(np.pad(faza,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)      |        # 012           
                                  #(np.pad(faza,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)       |         # 0_12
                                  (np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)         |         # 001
                                  (np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)       |        # 011
                                  (np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)        )      # 0_11      
               else:
                   CG14 = True
               faza002 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )&(np.pad(faza,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)& CG14, grain, faza) 
               seznam_premikov.append(faza002); del faza002

           elif s == '00_2':
               dij = np.array([0,0,-2])
               wij = en2*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_14:
                   CG14 = (#(np.pad(faza,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)    |       # 01_2
                                  #(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)     |       # 0_1_2
                                  (np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)       |       # 00_1
                                  (np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)     |      # 01_1
                                  (np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)      )    # 0_1_1    
               else:
                   CG14 = True
               faza00_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )&(np.pad(faza,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)& CG14, grain, faza)
               seznam_premikov.append(faza00_2); del faza00_2

           elif s == '020':
               dij = np.array([0,2,0])
               wij = en2*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_14:
                   CG14 = (#(np.pad(faza,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)     |         # 021
                                  #(np.pad(faza,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)      |       # 02_1
                                  (np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)        |         # 010 
                                  (np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)     |         # 011
                                  (np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)      )       # 01_1
               else:
                   CG14 = True
               faza020 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )&(np.pad(faza,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)& CG14, grain, faza)
               seznam_premikov.append(faza020); del faza020

           elif s == '0_20':
               dij = np.array([0,-2,0])
               wij = en2*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_14:
                   CG14 = (#(np.pad(faza,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)    |        # 0_21
                                  #(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)     |        # 0_2_1
                                  (np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)        |        # 0_10
                                  (np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)     |        # 0_11 
                                  (np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)      )      # 0_1_1
               else:
                   CG14 = True
               faza0_20 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )&(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)& CG14, grain, faza)     
               seznam_premikov.append(faza0_20); del faza0_20
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 15 ::: [200] >>> 1 site -------------------------------------------'''
           elif s == '200':
               dij = np.array([2,0,0])
               wij = en2*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza200 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )&(np.pad(faza,((2,0),(0,0), (0,0)), 'constant')[:-2,:,:]==grain)& (
                                                True), grain, faza)     
               seznam_premikov.append(faza200); del faza200
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 16 ::: [_200] >>> 1 site -------------------------------------------'''
           elif s == '_200':
               dij = np.array([-2,0,0])
               wij = en2*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_200 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )&(np.pad(faza,((0,2),(0,0), (0,0)), 'constant')[2:,:,:]==grain)& (
                                                True), grain, faza)     
               seznam_premikov.append(faza_200); del faza_200
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 17 ::: [022] >>> 4 sites -------------------------------------------'''
           elif s == '022':
               dij = np.array([0,2,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_17:
                   CG17 = (#(np.pad(faza,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)      |       # 012
                                  #(np.pad(faza,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)      |       # 021
                                  #(np.pad(faza,((0,0),(1,0), (3,0)), 'constant')[:,:-1,:-3]==grain)      |       # 013
                                  #(np.pad(faza,((0,0),(3,0), (1,0)), 'constant')[:,:-3,:-1]==grain)      |       # 031
                                  (np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)       )    # 011
               else:
                   CG17 = True
               faza022 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )&(np.pad(faza,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==grain)& CG17, grain, faza) 
               seznam_premikov.append(faza022); del faza022

           elif s == '02_2':
               dij = np.array([0,2,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_17:
                   CG17 = (#(np.pad(faza,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)    |        # 01_2
                                  #(np.pad(faza,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)    |        # 02_1
                                  #(np.pad(faza,((0,0),(1,0), (0,3)), 'constant')[:,:-1,3:]==grain)    |        # 01_3
                                  #(np.pad(faza,((0,0),(3,0), (0,1)), 'constant')[:,:-3,1:]==grain)   |        # 03_1
                                  (np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)     )      # 01_1
               else:
                   CG17 = True
               faza02_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )&(np.pad(faza,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==grain)& CG17, grain, faza) 
               seznam_premikov.append(faza02_2); del faza02_2

           elif s == '0_22':
               dij = np.array([0,-2,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_17:
                   CG17 = (#(np.pad(faza,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)     |        # 0_12
                                  #(np.pad(faza,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)     |        # 0_21
                                  #(np.pad(faza,((0,0),(0,1), (3,0)), 'constant')[:,1:,:-3]==grain)     |        # 0_13
                                  #(np.pad(faza,((0,0),(0,3), (1,0)), 'constant')[:,3:,:-1]==grain)     |        # 0_31
                                  (np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)      )     # 0_11
               else:
                   CG17 = True
               faza0_22 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )&(np.pad(faza,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==grain)& CG17, grain, faza) 
               seznam_premikov.append(faza0_22); del faza0_22

           elif s == '0_2_2':
               dij = np.array([0,-2,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               if constrains_of_group_17:
                   CG17 = (#(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)     |      # 0_1_2                
                                  #(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)     |      # 0_2_1
                                  #(np.pad(faza,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)     |      # 0_1_3
                                  #(np.pad(faza,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)     |     # 0_3_1
                                  (np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)      )    # 0_1_1
               else:
                   CG17 = True
               faza0_2_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )&(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)& CG17, grain, faza)
               seznam_premikov.append(faza0_2_2); del faza0_2_2
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 18 ::: [202], [220] >>> 4 sites -------------------------------------------'''
           elif s == '202':
               dij = np.array([2,0,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza202 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )&(np.pad(faza,((2,0),(0,0), (2,0)), 'constant')[:-2,:,:-2]==grain)& (                                                 
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(faza,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, faza)
               seznam_premikov.append(faza202); del faza202

           elif s == '20_2':
               dij = np.array([2,0,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza20_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )&(np.pad(faza,((2,0),(0,0), (0,2)), 'constant')[:-2,:,2:]==grain)& (                                                 
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(faza,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, faza)
               seznam_premikov.append(faza20_2); del faza20_2

           elif s == '220':
               dij = np.array([2,2,0])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza220 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )&(np.pad(faza,((2,0),(2,0), (0,0)), 'constant')[:-2,:-2,:]==grain)& (                                                 
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(faza,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, faza)
               seznam_premikov.append(faza220); del faza220

           elif s == '2_20':
               dij = np.array([2,-2,0])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza2_20 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )&(np.pad(faza,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)& (                                                 
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(faza,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, faza)
               seznam_premikov.append(faza2_20); del faza2_20
           

               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 19 ::: [_202], [_220] >>> 4 sites -------------------------------------------'''
           elif s == '_202':
               dij = np.array([-2,0,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_202 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )&(np.pad(faza,((0,2),(0,0), (2,0)), 'constant')[2:,:,:-2]==grain)& (                                                 
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(faza,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, faza)
               seznam_premikov.append(faza_202); del faza_202

           elif s == '_20_2':
               dij = np.array([-2,0,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_20_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )&(np.pad(faza,((0,2),(0,0), (0,2)), 'constant')[2:,:,2:]==grain)& (                                                 
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(faza,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, faza)
               seznam_premikov.append(faza_20_2); del faza_20_2

           elif s == '_220':
               dij = np.array([-2,2,0])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_220 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )&(np.pad(faza,((0,2),(2,0), (0,0)), 'constant')[2:,:-2,:]==grain)& (                                                 
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(faza,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, faza)
               seznam_premikov.append(faza_220); del faza_220

           elif s == '_2_20':
               dij = np.array([-2,-2,0])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_2_20 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )&(np.pad(faza,((0,2),(0,2), (0,0)), 'constant')[2:,2:,:]==grain)& (                                                 
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(faza,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(faza,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, faza)
               seznam_premikov.append(faza_2_20); del faza_2_20
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 20 ::: [112], [121] >>> 8 sites -------------------------------------------'''
           elif s == '112':
               dij = np.array([1,1,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza112 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==grain)& (                                                 
                                                                           #(np.pad(faza,((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==grain)       |      # 102                
                                                                           #(np.pad(faza,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)    |      # 111
                                                                           #(np.pad(faza,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)       |      # 012
                                                                           #(np.pad(faza,((2,0),(1,0), (2,0)), 'constant')[:-2,:-1,:-2]==grain)    |      # 212
                                                                           #(np.pad(faza,((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==grain)           # 122
                                               True ), grain, faza)
               seznam_premikov.append(faza112); del faza112

           elif s == '11_2':
               dij = np.array([1,1,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza11_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==grain)& (                                                 
                                                                           #(np.pad(faza,((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==grain)       |      # 10_2                
                                                                           #(np.pad(faza,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)    |      # 11_1
                                                                           #(np.pad(faza,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)       |      # 01_2
                                                                           #(np.pad(faza,((2,0),(1,0), (0,2)), 'constant')[:-2,:-1,2:]==grain)    |      # 21_2
                                                                           #(np.pad(faza,((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==grain)           # 12_2
                                                True), grain, faza)
               seznam_premikov.append(faza11_2); del faza11_2

           elif s == '1_12':
               dij = np.array([1,-1,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza1_12 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==grain)& (                                                 
                                                                           #(np.pad(faza,((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==grain)      |      # 102                
                                                                           #(np.pad(faza,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)    |      # 1_11
                                                                           #(np.pad(faza,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)       |      # 0_12
                                                                           #(np.pad(faza,((2,0),(0,1), (2,0)), 'constant')[:-2,1:,:-2]==grain)    |      # 2_12
                                                                           #(np.pad(faza,((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==grain)           # 1_22
                                                True), grain, faza)
               seznam_premikov.append(faza1_12); del faza1_12

           elif s == '1_1_2':
               dij = np.array([1,-1,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza1_1_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)& (                                                 
                                                                           #(np.pad(faza,((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==grain)      |      # 10_2                
                                                                           #(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)    |      # 1_1_1
                                                                           #(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)       |      # 0_1_2
                                                                           #(np.pad(faza,((2,0),(0,1), (0,2)), 'constant')[:-2,1:,2:]==grain)    |      # 2_1_2
                                                                           #(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)           # 1_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza1_1_2); del faza1_1_2

           elif s == '121':
               dij = np.array([1,2,1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza121 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==grain)& (                                                 
                                                                           #(np.pad(faza,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)      |      # 111                
                                                                           #(np.pad(faza,((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==grain)         |      # 120
                                                                           #(np.pad(faza,((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==grain)      |      # 122
                                                                           #(np.pad(faza,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)         |      # 021
                                                                           #(np.pad(faza,((2,0),(2,0), (1,0)), 'constant')[:-2,:-2,:-1]==grain)             # 221
                                                True), grain, faza)
               seznam_premikov.append(faza121); del faza121

           elif s == '12_1':
               dij = np.array([1,2,-1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza12_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==grain)& (                                                 
                                                                           #(np.pad(faza,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)      |      # 11_1                
                                                                           #(np.pad(faza,((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==grain)        |      # 120
                                                                           #(np.pad(faza,((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==grain)      |      # 12_2
                                                                           #(np.pad(faza,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)         |      # 02_1
                                                                           #(np.pad(faza,((2,0),(2,0), (0,1)), 'constant')[:-2,:-2,1:]==grain)             # 22_1
                                                True), grain, faza)
               seznam_premikov.append(faza12_1); del faza12_1

           elif s == '1_21':
               dij = np.array([1,-2,1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza1_21 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)& (                                                 
                                                                           #(np.pad(faza,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)      |      # 1_11                
                                                                           #(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)         |      # 1_20
                                                                           #(np.pad(faza,((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==grain)      |      # 1_22
                                                                           #(np.pad(faza,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)         |      # 0_21
                                                                           #(np.pad(faza,((2,0),(0,2), (1,0)), 'constant')[:-2,2:,:-1]==grain)             # 2_21
                                                True), grain, faza)
               seznam_premikov.append(faza1_21); del faza1_21

           elif s == '1_2_1':
               dij = np.array([1,-2,-1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza1_2_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)& (                                                 
                                                                           #(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #(np.pad(faza,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, faza)
               seznam_premikov.append(faza1_2_1); del faza1_2_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 21 ::: [_112], [_121] >>> 8 sites -------------------------------------------'''
           elif s == '_112':
               dij = np.array([-1,1,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_112 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((0,1),(1,0), (2,0)), 'constant')[1:,:-1,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==grain)       |      # 102                
                                                                           #wrong(np.pad(faza,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)    |      # 111
                                                                           #wrong(np.pad(faza,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)       |      # 012
                                                                           #wrong(np.pad(faza,((2,0),(1,0), (2,0)), 'constant')[:-2,:-1,:-2]==grain)    |      # 212
                                                                           #wrong(np.pad(faza,((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==grain)           # 122
                                               True ), grain, faza)
               seznam_premikov.append(faza_112); del faza_112

           elif s == '_11_2':
               dij = np.array([-1,1,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_11_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((0,1),(1,0), (0,2)), 'constant')[1:,:-1,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==grain)       |      # 10_2                
                                                                           #wrong(np.pad(faza,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)    |      # 11_1
                                                                           #wrong(np.pad(faza,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)       |      # 01_2
                                                                           #wrong(np.pad(faza,((2,0),(1,0), (0,2)), 'constant')[:-2,:-1,2:]==grain)    |      # 21_2
                                                                           #wrong(np.pad(faza,((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==grain)           # 12_2
                                                True), grain, faza)
               seznam_premikov.append(faza_11_2); del faza_11_2

           elif s == '_1_12':
               dij = np.array([-1,-1,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_1_12 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((0,1),(0,1), (2,0)), 'constant')[1:,1:,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==grain)      |      # 102                
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)    |      # 1_11
                                                                           #wrong(np.pad(faza,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)       |      # 0_12
                                                                           #wrong(np.pad(faza,((2,0),(0,1), (2,0)), 'constant')[:-2,1:,:-2]==grain)    |      # 2_12
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==grain)           # 1_22
                                                True), grain, faza)
               seznam_premikov.append(faza_1_12); del faza_1_12

           elif s == '_1_1_2':
               dij = np.array([-1,-1,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_1_1_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((0,1),(0,1), (0,2)), 'constant')[1:,1:,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==grain)      |      # 10_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)    |      # 1_1_1
                                                                           #wrong(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)       |      # 0_1_2
                                                                           #wrong(np.pad(faza,((2,0),(0,1), (0,2)), 'constant')[:-2,1:,2:]==grain)    |      # 2_1_2
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)           # 1_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_1_1_2); del faza_1_1_2

           elif s == '_121':
               dij = np.array([-1,2,1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_121 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((0,1),(2,0), (1,0)), 'constant')[1:,:-2,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)      |      # 111                
                                                                           #wrong(np.pad(faza,((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==grain)         |      # 120
                                                                           #wrong(np.pad(faza,((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==grain)      |      # 122
                                                                           #wrong(np.pad(faza,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)         |      # 021
                                                                           #wrong(np.pad(faza,((2,0),(2,0), (1,0)), 'constant')[:-2,:-2,:-1]==grain)             # 221
                                                True), grain, faza)
               seznam_premikov.append(faza_121); del faza_121

           elif s == '_12_1':
               dij = np.array([-1,2,-1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_12_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((0,1),(2,0), (0,1)), 'constant')[1:,:-2,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)      |      # 11_1                
                                                                           #wrong(np.pad(faza,((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==grain)        |      # 120
                                                                           #wrong(np.pad(faza,((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==grain)      |      # 12_2
                                                                           #wrong(np.pad(faza,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)         |      # 02_1
                                                                           #wrong(np.pad(faza,((2,0),(2,0), (0,1)), 'constant')[:-2,:-2,1:]==grain)             # 22_1
                                                True), grain, faza)
               seznam_premikov.append(faza_12_1); del faza_12_1

           elif s == '_1_21':
               dij = np.array([-1,-2,1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_1_21 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((0,1),(0,2), (1,0)), 'constant')[1:,2:,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)      |      # 1_11                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)         |      # 1_20
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==grain)      |      # 1_22
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)         |      # 0_21
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (1,0)), 'constant')[:-2,2:,:-1]==grain)             # 2_21
                                                True), grain, faza)
               seznam_premikov.append(faza_1_21); del faza_1_21

           elif s == '_1_2_1':
               dij = np.array([-1,-2,-1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_1_2_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((0,1),(0,2), (0,1)), 'constant')[1:,2:,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, faza)
               seznam_premikov.append(faza_1_2_1); del faza_1_2_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 22 ::: [211] >>> 4 sites -------------------------------------------'''
           elif s == '211':
               dij = np.array([2,1,1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza211 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((2,0),(1,0), (1,0)), 'constant')[:-2,:-1,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, faza)
               seznam_premikov.append(faza211); del faza211

           elif s == '21_1':
               dij = np.array([2,1,-1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza21_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((2,0),(1,0), (0,1)), 'constant')[:-2,:-1,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, faza)
               seznam_premikov.append(faza21_1); del faza21_1

           elif s == '2_11':
               dij = np.array([2,-1,1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza2_11 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((2,0),(0,1), (1,0)), 'constant')[:-2,1:,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, faza)
               seznam_premikov.append(faza2_11); del faza2_11

           elif s == '2_1_1':
               dij = np.array([2,-1,-1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza2_1_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((2,0),(0,1), (0,1)), 'constant')[:-2,1:,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, faza)
               seznam_premikov.append(faza2_1_1); del faza2_1_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 23 ::: [_211] >>> 4 sites -------------------------------------------'''
           elif s == '_211':
               dij = np.array([-2,1,1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_211 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((0,2),(1,0), (1,0)), 'constant')[2:,:-1,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, faza)
               seznam_premikov.append(faza_211); del faza_211

           elif s == '_21_1':
               dij = np.array([-2,1,-1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_21_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((0,2),(1,0), (0,1)), 'constant')[2:,:-1,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, faza)
               seznam_premikov.append(faza_21_1); del faza_21_1

           elif s == '_2_11':
               dij = np.array([-2,-1,1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_2_11 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((0,2),(0,1), (1,0)), 'constant')[2:,1:,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, faza)
               seznam_premikov.append(faza_2_11); del faza_2_11

           elif s == '_2_1_1':
               dij = np.array([-2,-1,-1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_2_1_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_7th(R) )&(np.pad(faza,((0,2),(0,1), (0,1)), 'constant')[2:,1:,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, faza)
               seznam_premikov.append(faza_2_1_1); del faza_2_1_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 24 ::: [122] >>> 4 sites -------------------------------------------'''
           elif s == '122':
               dij = np.array([1,2,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza122 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==grain)& (                                                 
                                                                           #(np.pad(faza,((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==grain)      |      # 112                
                                                                           #(np.pad(faza,((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==grain)      |      # 121
                                                                           #(np.pad(faza,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==grain)         |      # 022
                                                                           #(np.pad(faza,((2,0),(2,0), (2,0)), 'constant')[:-2,:-2,:-2]==grain)             # 222
                                                True), grain, faza)
               seznam_premikov.append(faza122); del faza122

           elif s == '12_2':
               dij = np.array([1,2,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza12_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==grain)& (                                                 
                                                                           #(np.pad(faza,((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==grain)      |      # 11_2                
                                                                           #(np.pad(faza,((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==grain)      |      # 12_1
                                                                           #(np.pad(faza,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==grain)         |      # 02_2
                                                                           #(np.pad(faza,((2,0),(2,0), (0,2)), 'constant')[:-2,:-2,2:]==grain)             # 22_2
                                                True), grain, faza)
               seznam_premikov.append(faza12_2); del faza12_2

           elif s == '1_22':
               dij = np.array([1,-2,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza1_22 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==grain)& (                                                 
                                                                           #(np.pad(faza,((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==grain)      |      # 1_12                
                                                                           #(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)      |      # 1_21
                                                                           #(np.pad(faza,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==grain)         |      # 0_22
                                                                           #(np.pad(faza,((2,0),(0,2), (2,0)), 'constant')[:-2,2:,:-2]==grain)             # 2_22
                                                True), grain, faza)
               seznam_premikov.append(faza1_22); del faza1_22

           elif s == '1_2_2':
               dij = np.array([1,-2,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza1_2_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)& (                                                 
                                                                           #(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza1_2_2); del faza1_2_2
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 25 ::: [_122] >>> 4 sites -------------------------------------------'''
           elif s == '_122':
               dij = np.array([-1,2,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_122 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((0,1),(2,0), (2,0)), 'constant')[1:,:-2,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==grain)      |      # 112                
                                                                           #wrong(np.pad(faza,((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==grain)      |      # 121
                                                                           #wrong(np.pad(faza,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==grain)         |      # 022
                                                                           #wrong(np.pad(faza,((2,0),(2,0), (2,0)), 'constant')[:-2,:-2,:-2]==grain)             # 222
                                                True), grain, faza)
               seznam_premikov.append(faza_122); del faza_122

           elif s == '_12_2':
               dij = np.array([-1,2,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_12_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((0,1),(2,0), (0,2)), 'constant')[1:,:-2,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==grain)      |      # 11_2                
                                                                           #wrong(np.pad(faza,((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==grain)      |      # 12_1
                                                                           #wrong(np.pad(faza,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==grain)         |      # 02_2
                                                                           #wrong(np.pad(faza,((2,0),(2,0), (0,2)), 'constant')[:-2,:-2,2:]==grain)             # 22_2
                                                True), grain, faza)
               seznam_premikov.append(faza_12_2); del faza_12_2

           elif s == '_1_22':
               dij = np.array([-1,-2,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_1_22 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((0,1),(0,2), (2,0)), 'constant')[1:,2:,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==grain)      |      # 1_12                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)      |      # 1_21
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==grain)         |      # 0_22
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (2,0)), 'constant')[:-2,2:,:-2]==grain)             # 2_22
                                                True), grain, faza)
               seznam_premikov.append(faza_1_22); del faza_1_22

           elif s == '_1_2_2':
               dij = np.array([-1,-2,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_1_2_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((0,1),(0,2), (0,2)), 'constant')[1:,2:,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_1_2_2); del faza_1_2_2

               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 26 ::: [212], [221] >>> 8 sites -------------------------------------------'''
           elif s == '212':
               dij = np.array([2,1,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza212 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((2,0),(1,0), (2,0)), 'constant')[:-2,:-1,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza212); del faza212

           elif s == '21_2':
               dij = np.array([2,1,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza21_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((2,0),(1,0), (0,2)), 'constant')[:-2,:-1,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza21_2); del faza21_2

           elif s == '2_12':
               dij = np.array([2,-1,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza2_12 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((2,0),(0,1), (2,0)), 'constant')[:-2,1:,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza2_12); del faza2_12

           elif s == '2_1_2':
               dij = np.array([2,-1,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza2_1_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((2,0),(0,1), (0,2)), 'constant')[:-2,1:,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza2_1_2); del faza2_1_2

           elif s == '221':
               dij = np.array([2,2,1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza221 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((2,0),(2,0), (1,0)), 'constant')[:-2,:-2,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza221); del faza221

           elif s == '22_1':
               dij = np.array([2,2,-1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza22_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((2,0),(2,0), (0,1)), 'constant')[:-2,:-2,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza22_1); del faza22_1

           elif s == '2_21':
               dij = np.array([2,-2,1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza2_21 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((2,0),(0,2), (1,0)), 'constant')[:-2,2:,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza2_21); del faza2_21

           elif s == '2_2_1':
               dij = np.array([2,-2,-1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza2_2_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza2_2_1); del faza2_2_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 27 ::: [_212], [_221] >>> 8 sites -------------------------------------------'''
           elif s == '_212':
               dij = np.array([-2,1,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_212 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((0,2),(1,0), (2,0)), 'constant')[2:,:-1,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_212); del faza_212

           elif s == '_21_2':
               dij = np.array([-2,1,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_21_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((0,2),(1,0), (0,2)), 'constant')[2:,:-1,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_21_2); del faza_21_2

           elif s == '_2_12':
               dij = np.array([-2,-1,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_2_12 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((0,2),(0,1), (2,0)), 'constant')[2:,1:,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_2_12); del faza_2_12

           elif s == '_2_1_2':
               dij = np.array([-2,-1,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_2_1_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((0,2),(0,1), (0,2)), 'constant')[2:,1:,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_2_1_2); del faza_2_1_2

           elif s == '_221':
               dij = np.array([-2,2,1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_221 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((0,2),(2,0), (1,0)), 'constant')[2:,:-2,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_221); del faza_221

           elif s == '_22_1':
               dij = np.array([-2,2,-1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_22_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((0,2),(2,0), (0,1)), 'constant')[2:,:-2,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_22_1); del faza_22_1

           elif s == '_2_21':
               dij = np.array([-2,-2,1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_2_21 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((0,2),(0,2), (1,0)), 'constant')[2:,2:,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_2_21); del faza_2_21

           elif s == '_2_2_1':
               dij = np.array([-2,-2,-1])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_2_2_1 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_8th(R) )&(np.pad(faza,((0,2),(0,2), (0,1)), 'constant')[2:,2:,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_2_2_1); del faza_2_2_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 28 ::: [222] >>> 4 sites -------------------------------------------'''
           elif s == '222':
               dij = np.array([2,2,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza222 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_9th(R) )&(np.pad(faza,((2,0),(2,0), (2,0)), 'constant')[:-2,:-2,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza222); del faza222

           elif s == '22_2':
               dij = np.array([2,2,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza22_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_9th(R) )&(np.pad(faza,((2,0),(2,0), (0,2)), 'constant')[:-2,:-2,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza22_2); del faza22_2

           elif s == '2_22':
               dij = np.array([2,-2,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza2_22 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_9th(R) )&(np.pad(faza,((2,0),(0,2), (2,0)), 'constant')[:-2,2:,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza2_22); del faza2_22

           elif s == '2_2_2':
               dij = np.array([2,-2,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza2_2_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_9th(R) )&(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza2_2_2); del faza2_2_2
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               # *** UNDER CONSTRUCTION ! ***
               ''' ----------------------------------------- GROUP 29 ::: [_222] >>> 4 sites -------------------------------------------'''
           elif s == '_222':
               dij = np.array([-2,2,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_222 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_9th(R) )&(np.pad(faza,((0,2),(2,0), (2,0)), 'constant')[2:,:-2,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_222); del faza_222

           elif s == '_22_2':
               dij = np.array([-2,2,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_22_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_9th(R) )&(np.pad(faza,((0,2),(2,0), (0,2)), 'constant')[2:,:-2,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_22_2); del faza_22_2

           elif s == '_2_22':
               dij = np.array([-2,-2,2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_2_22 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_9th(R) )&(np.pad(faza,((0,2),(0,2), (2,0)), 'constant')[2:,2:,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_2_22); del faza_2_22

           elif s == '_2_2_2':
               dij = np.array([-2,-2,-2])
               wij = en3*W(asc[grain]['oi'], dij)
               lij = wij*vg*S[s]
               faza_2_2_2 = np.where((T_next < T)&(faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_9th(R) )&(np.pad(faza,((0,2),(0,2), (0,2)), 'constant')[2:,2:,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(faza,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(faza,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(faza,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, faza)
               seznam_premikov.append(faza_2_2_2); del faza_2_2_2



           if GDSM == 'new':
               #*********************************************************************** NEW ***************************************************************************************************
               faza = seznam_premikov[0]
               cas = S[s]
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
            negind -=1
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
        print('CPU (%): ', psutil.cpu_percent(),' %       RAM (%): ',psutil.virtual_memory().percent,' %')
        print(30*' .  . ')

        PD['i'].append(i)
        PD['current step CPU'].append(step_cpu)
        PD['ALL grains'].append(grain_counter)

    #if grain_counter and not bool(grain_ID.shape[0])  :
     #  break

    if i>=END_step:
        break
    

process_time = time.time() - start_time
print();print(); print()
print(20*' ',100*'*')
print(25*' ','DONE!!  Total computing time =  ',round(process_time, 3),'  seconds.   (RATE = ', round(process_time/h, 4),' s / step)'); print(25*' ',100*'-')
print(25*' ','Printed length real dimension = ',round(FEM_scanning_speed*dt*h, 5),u'\u00B5m      |    Printing time = ',round(1000*dt*h, 5),' msec.')
print(25*' ',100*'-')
print(65*' ','Number of grains =  ',grain_counter,' !'); print(20*' ',100*'*'); print(); print(); print()

# ......................................... PRIKAZ in IZPIS REZULTATOV ..................................................#
"""
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
"""
''' ---------------------------------------------------------------------------------- the end -------------------------------------------------------------------------------------------------------------- '''






