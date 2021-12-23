import numpy as np
import math
import matplotlib.pyplot as plt
from numpy.linalg import multi_dot
plt.ion()


Oi001 = np.array([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
    ])

Oi101 = np.array([
    [math.cos(math.radians(45)), 0, math.cos(math.radians(45))],
    [0, 1, 0],
    [-math.cos(math.radians(45)), 0, math.cos(math.radians(45))],
    ])

k=35.2643896827

Oi111 = np.array([
    [math.cos(math.radians(k))* math.cos(math.radians(30)), - math.cos(math.radians(k))* math.cos(math.radians(60)), math.cos(math.radians(90-k))],
    [0, math.cos(math.radians(k)), math.cos(math.radians(90-k))],
    [- math.cos(math.radians(k))* math.cos(math.radians(30)), -math.cos(math.radians(k))* math.cos(math.radians(60)), math.cos(math.radians(90-k))],
    ])

'''--------------------------------------------------------------------------------------------------------------------------------- funkcije --------------------------------------------------------------'''
def Get_001_101_111_vectors(m, axis):
    global red, green, blue
    ax = {'x':0, 'y':1, 'z':2}
    Qi = np.array([
        m[2],
        np.array([m[0,0]+m[2,0], m[0,1]+m[2,1], m[0,2]+m[2,2]]),
        np.array([m[0,0]+m[1,0]+m[2,0], m[0,1]+m[1,1]+m[2,1], m[0,2]+m[1,2]+m[2,2]]),
        ])
    Qi = np.around(Qi, decimals=4)
    Qi_lengths = np.linalg.norm(Qi, axis=1)
    Qi_axed = Qi[:,ax[axis]]

    R_vector = np.sum(Qi, axis=0); print('R_vector is:  ', R_vector)

    return Qi, Qi_lengths, R_vector


def Get_Vector_Sum(m):
    sum_m = np.array([m[0,0]+m[1,0]+m[2,0], m[0,1]+m[1,1]+m[2,1], m[0,2]+m[1,2]+m[2,2]])
    return sum_m

def Get_Projection_Height(v1, v2):
    ph = np.dot(v1,v2)
    return ph


def Get_Vector_Angle(v1, v2, acute):
    # v1 is first vector
    # v2 is second vector
    angle_rad = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))    # angle between vectors v1 and v2 in radians!
    angle = np.degrees(angle_rad)
    if (acute == True):
        return np.around(angle, decimals=2)
    else:
        return np.around(2 * np.pi - angle, decimals=2)


def Rotate_the_Cube(alpha, beta):
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


def Custom_Grain(alpha,beta):
    global Oicus, Qicus
    Oicus = Rotate_the_Cube(alpha,beta)
    Qicus = Get_001_101_111_vectors(Oicus, 'z')
    koticus=[]; phcus = []; rvcus = Qicus[2]
    for i in range(3):
        koticus.append(Get_Vector_Angle(os_opazovanja, Qicus[0][i], True)  )
        phcus.append(Get_Projection_Height(os_opazovanja, Qicus[0][i]))
    
    return np.array(phcus)


    
'''----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'''

os_opazovanja = np.array([0,0,1])

Qi001 = Get_001_101_111_vectors(Oi001,'z')
Qi101 = Get_001_101_111_vectors(Oi101, 'z')
Qi111 = Get_001_101_111_vectors(Oi111, 'z')

koti001=[]; koti101=[]; koti111=[]
ph001 = []; ph101=[]; ph111=[]
rv001 = Qi001[2]; rv101 = Qi101[2]; rv111 = Qi111[2]; 

for i in range(3):
    koti001.append(Get_Vector_Angle(os_opazovanja, Qi001[0][i], True)  )
    koti101.append(Get_Vector_Angle(os_opazovanja, Qi101[0][i], True)  )
    koti111.append(Get_Vector_Angle(os_opazovanja, Qi111[0][i], True)  )
    
    ph001.append(Get_Projection_Height(os_opazovanja, Qi001[0][i]))
    ph101.append(Get_Projection_Height(os_opazovanja, Qi101[0][i]))
    ph111.append(Get_Projection_Height(os_opazovanja, Qi111[0][i]))

    

P = np.transpose(np.array([np.array(ph001), np.array(ph101), np.array(ph111)]))
#P = np.transpose(np.array([np.array(rv001), np.array(rv101), np.array(rv111)]))

''' ~~~~~~~~~~~~~~~~~~~ custom crystal orientation ~~~~~~~~~~~~~~~~~~~~~ '''

alfas =  [22.5]              #range(0,50,5)
betas = [0]             #range(0,36)

asc = {}
grain_counter=0

R,G,B = [], [], []
all_counter=0; poz_counter=0
for i in alfas:
    for j in betas:
        grain_counter+=1
        asc[grain_counter] ={     'oi': Rotate_the_Cube(i,j), 'alpha':i, 'beta':j,}
        C = Custom_Grain(i,j)
        RGB = np.linalg.solve(np.array([P[0], P[1], P[2]]), C)
        if np.all(RGB>=0):
            #print(RGB)
            #print(i,j); print()
            poz_counter+=1
        R.append(RGB[0]); G.append(RGB[1]); B.append(RGB[2])
        all_counter+=1

#plt.plot(R, color='#ff0000', marker='o')
#plt.plot(G, color='#00ff00', marker='o')
#plt.plot(B, color='#0000ff', marker='o')


slika=np.zeros((3,3,3)).astype('int')
rgb = np.around(255*RGB, decimals=0)
slika[1,1]=rgb
plt.figure()
plt.imshow(slika)

















""" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ """
'''
x_orient = ['001', '101', '111']

plt.plot(x_orient, ph001, label='RED --- Oi001',marker='s', color='#ff0000'); #print(koti001)
plt.plot(x_orient, ph101, label='GREEN --- Oi101',marker='^', color='#00ff00'); #print(koti101)
plt.plot(x_orient, ph111, label='BLUE --- Oi111',marker='o', color='#0000ff'); #print(koti111)

plt.plot(x_orient, phcus, label='custom',marker='.', color='magenta'); #print(koticus)

plt.legend()


K = np.array([np.array(koti001), np.array(koti101), np.array(koti111), np.array(koticus)])
K_delta = np.array([
                np.absolute(K[:3,0]-K[3,0]),
                np.absolute(K[:3,1]-K[3,1]),
                np.absolute(K[:3,2]-K[3,2]),
                ])
'''

"""-------------------------------------- ListedColormap for Coloring grains according to orientation --------------------------------------------------"""
"""
from matplotlib import cm
from matplotlib.colors import ListedColormap

viridis = cm.get_cmap('viridis', 5)
newcolors = viridis(np.linspace(0, 1, 5))

slika = np.zeros((20,20))

slika[15,15]=1
slika[5,5]=2
slika[15,7]=3
slika[10,10]=4


asc = {1: {'rgb': (51, 204, 255), 'name': 'pastel_blue', 'hex': '#33ccff'},
            2: {'rgb': (204, 153, 255), 'name': 'pastel_purple', 'hex': '#cc99ff'},
            3:  {'rgb': (255, 204, 153), 'name': 'pastel_orange', 'hex': '#ffcc99'},
            4:  {'rgb': (153, 255, 153), 'name': 'pastel_green', 'hex': '#99ff99'},
       }

stevc = 1
for grain in asc:
    c = asc[grain]['rgb']
    color = np.array([c[0]/255, c[1]/255, c[2]/255, 1])
    newcolors[stevc, :] = color
    stevc+=1

newcmap = ListedColormap(newcolors)
plt.imshow(slika, cmap=newcmap)
"""





