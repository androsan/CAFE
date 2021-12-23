import json
import numpy as np
from numpy.linalg import multi_dot
import random, math, time
import matplotlib.pyplot as plt
plt.ion()

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
    xran = np.random.randint(0, posibilities+1)
    x = xran / posibilities
    y = np.random.randint(0, xran+1)/posibilities
    return x,y

rp = 90

Z=1
X=96
Y=128

case =      'SLM_2D_Source'   ;  subcase = '0002'
PATH =     'C:/sm-2018-w64-0-3/WORK/'+case+'_Files/post processing database/'+subcase+'/'

mapa    =   'INTER  time=40, space=8  Z[7-8], X[15-27], Y[12-97], 1500°C, N=12/'

cut_BOTTOM = 'flashies_faza'   ;  cut_UP =   ''

X_cutoff_percent = 75       ;  X_cutoff = int(X_cutoff_percent * X / 100)      # Analogy with Hatch Spacing

yp_count =  0
cutoff_limit = 64

FS_left   =  np.load(PATH+mapa+cut_BOTTOM+'/results/cut_'+str(yp_count)   +'.npy')[:,:X_cutoff,:]            # FS stands for "flashy snap"
FS_right =  np.load(PATH+mapa+cut_BOTTOM+'/results/cut_'+str(yp_count+1)+'.npy')[:,:X_cutoff,:cutoff_limit]
FS = np.dstack((FS_left, FS_right))
grain_ID = list(np.unique(FS).astype('int')); grain_ID.remove(0)

HA = np.zeros((Z,X-X_cutoff, Y))

faza = np.hstack((HA,FS))



''' Construction of arbitrary asc dictionary (asc_arb) --- in real process, for each item in grain_ID data from nuclei_data.json will be imported to form new asc '''

asc_arb = {}
for j in grain_ID:
    """ generation of random grain orientation """
    x,y = random_xy_tilt(rp)
    rgb = get_color(x,y)
    alfa  = math.degrees(math.atan(x))- 9.7356103173*x*y
    beta = math.degrees(math.atan(y))- 9.7356103173*x*y
    cub_xy = Rotate_the_Cube_XY(alfa, beta)
    gama = math.degrees(math.atan(np.random.randint(0, rp+1)/rp))
    oi = Rotate_the_Cube_Z(cub_xy, gama)
    asc_arb[j] ={'oi': oi, 'alfa':alfa, 'beta':beta, 'gama':gama, 'rgb':rgb, }    # ALL data about nuclei


with open(PATH+mapa+'nuclei_data.json', 'w') as nuks:              # Writing data of ALL nuclei as .json file, but values must be list NOT np.array !!!
    asc_list =asc_arb.copy()
    for nuk in asc_list:
        asc_list[nuk]['oi']=asc_arb[nuk]['oi'].tolist()                                
        asc_list[nuk]['rgb']=asc_arb[nuk]['rgb'].tolist()
    asc_list = {str(keyy): vall for keyy, vall in asc_list.items()}
    json.dump(asc_list, nuks)


''' Loading asc from nuclei_data.json'''
with open(PATH+mapa+'nuclei_data.json', 'r') as nuks:
    asc = json.load(nuks)
    asc = {int(keyy): vall for keyy, vall in asc.items()}
    for nuk in asc:
        asc[nuk]['oi']=np.array(asc[nuk]['oi'])
        asc[nuk]['rgb']=np.array(asc[nuk]['rgb'])


plt.imshow(faza[0])
