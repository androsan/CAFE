from numpy.linalg import multi_dot
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.ion()

k=35.2643896827

def get_color(x,y,z):
    u = z-x
    v = x-y
    w = y
    rgb = np.array([u,v,w])
    RGB = 255*rgb/np.max(rgb)
    return RGB.astype('int')

def random_xy_tilt():
    x = np.random.uniform(0,1)
    y = np.random.uniform(0,x)
    return x,y
   
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

    return oiz

#++++++++++++++++++++++ exe line ++++++++++++++++++++++++++++++++++++++++

''' tilting the cube '''
#x,y = random_xy_tilt()
x,y = 1,0
RGB = get_color(x,y,1)
alfa  = math.degrees(math.atan(x))
beta = math.degrees(math.atan(y))- 9.7356103173*y
cub_xy = Rotate_the_Cube_XY(alfa, beta)
gama = math.degrees(math.atan(np.random.uniform(0,1)))
c = Rotate_the_Cube_Z(cub_xy, gama)


''' showing the tilted cube '''

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


''' get the grain oriented color '''
slika=np.zeros((3,3,3)).astype('int')
slika[1,1]=RGB
plt.figure()
plt.imshow(slika)


#--------------------------------------------------------------------------------------------------------------------
'''
R=get_color(0,0,1)
G=get_color(1,0,1)
B=get_color(1,1,1)

Y=get_color(0.5,0,1)
P=get_color(0.5,0.5,1)
T=get_color(1,0.5,1)
'''

    
    
