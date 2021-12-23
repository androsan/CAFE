import numpy as np

Z=2
X=3
Y=3

z_ = np.arange(Z)
x_ = np.arange(X)
y_ = np.arange(Y)

z, x, y = np.meshgrid(z_, x_, y_, indexing='ij')
assert np.all(z[:,0,0] == z_)
assert np.all(x[0,:,0] == x_)
assert np.all(y[0,0,:] == y_)

ZXY = np.column_stack([z.flat, x.flat, y.flat])
ZXY = ZXY.reshape((Z, X, Y, 3))

CA = ZXY+0.5                                                        #  Creating matrix of coordinates (z,x,y) of nodes in the center of CA cells

epsilon = 0.49                                                         #  Scaling value, should be in between 0 and 0.49
rand = np.random.uniform(0,1)                              # Random float between 0 and 1
R = np.random.uniform(0, 1, (Z,X,Y,3))              # Matrix of random float values between 0 and 1

PA = np.zeros((Z,X,Y,3))                                 #  Creating matrix of coordinates (z,x,y) of nodes displaced from the CA cells center for random z,x,y values
PA[:,:,:,0] = CA[:,:,:,0] + epsilon*(2*R[:,:,:,0] - 1)
PA[:,:,:,1] = CA[:,:,:,1] + epsilon*(2*R[:,:,:,1] - 1)
PA[:,:,:,2] = CA[:,:,:,2] + epsilon*(2*R[:,:,:,2] - 1)



