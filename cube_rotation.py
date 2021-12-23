import math
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import multi_dot
from mpl_toolkits.mplot3d import Axes3D
plt.ion()

cube = np.array([
                              [1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1]
                            ])

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


"""********************"""
alfa =      45
beta =     45
"""**************"""

c = Rotate_the_Cube(alfa,beta)
ca = np.array([[0,0,0], c[0], c[1], c[2], c[0]+c[1], c[1]+c[2], c[0]+c[2], c[0]+c[1]+c[2]])
ogli = [(0,1), (0,2), (0,3), (1,4), (2,4), (2,5), (4,7),  (1,6), (3,6), (3,5), (5,7), (6,7)]

def edges(p,q):
    edge = [[ca[p][0], ca[q][0]], [ca[p][1], ca[q][1]], [ca[p][2], ca[q][2]]]
    return edge

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

'''------------- fake bounding box ----------------'''
ax.set_aspect('equal')
max_range = 2
Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() 
Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten()
Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() 
for xb, yb, zb in zip(Xb, Yb, Zb):
   ax.plot([xb], [yb], [zb], 'w')
''' -----------------------------------------------------'''

ax.scatter(ca[:,0], ca[:,1], ca[:,2], marker='o', s=80, color='blue')
for i in ogli:
    rob = edges(i[0],i[1])
    ax.plot(rob[0], rob[1], rob[2], color='red')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
