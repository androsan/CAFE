import numpy as np
'''
a= np.zeros((5,5,5))

a[0,0,1]=1
a[1,0,1]=1
a[1,1,1]=1
a[1,1,0]=1
a[1,1,2]=1

n = []
for i in range(a.shape[0]):
    t=np.where(np.all(a[i]==0, axis=0))
    n.append([t[0], len(t[0])])


m=np.array(n)
v = m[np.min(m[:,1])==m[:,1]][0][0]

#v=np.array([0,1,2,3,8,9,10,11])
#v=np.array([8,9,10,11])

if v[0]!=0:
    v=np.insert(v, 0, 0)
    dx = np.diff(v)
    bol=v[:-1][dx!=1]
    limits=(None, v[bol[0]+1])

else:
    dx = np.diff(v)
    bol=v[:-1][dx!=1]
    limits=(v[bol[0]]+1, v[bol[0]+1])



b = a[:,:, limits[0]: limits[1]]


for i in range(b.shape[0]):
    t=np.where(np.all(b[i]==0, axis=1))
    print(t[0])
'''

d=np.zeros((3,3,3))

d[0,1,1]=1
d[0,0,1]=1

d[1,1,1]=1
d[1,1,2]=1

d[2,1,2]=1
d[2,1,1]=1
d[2,1,0]=1



ax2=np.where(np.any(d==1, axis=2))
ax1=np.where(np.any(d==1, axis=1))
ax0=np.where(np.any(d==1, axis=0))

Z_lim = (np.min(ax2[0]), np.max(ax2[0])+1)
X_lim = (np.min(ax2[1]), np.max(ax2[1])+1)
Y_lim = (np.min(ax1[1]), np.max(ax1[1])+1)

f= d[Z_lim[0]:Z_lim[1], X_lim[0]:X_lim[1], Y_lim[0]:Y_lim[1]]











