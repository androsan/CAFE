import numpy as np
import math
import matplotlib.pyplot as plt
plt.ion()


rp = 100

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


#x,y = random_xy_tilt(rp)

x = 0.5
y = 0

rgb = get_color(x,y)
alfa  = math.degrees(math.atan(x))- 9.7356103173*x*y
beta = math.degrees(math.atan(y))- 9.7356103173*x*y

A = np.zeros((1,1,3)).astype('int')
A[0,0]=rgb
plt. imshow(A)
plt.title(u"\u03B1 = "+str(round(alfa, 3))+u"\u00B0"+'  ,   '+u"\u03B2 = "+str(round(beta, 3))+u"\u00B0", fontsize=18)

