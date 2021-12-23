from multiprocessing import Pool, Manager
from functools import partial
import numpy as np

def G(a, b):
    Q=b[0]
    R=b[1]
    
    if R == '001':
        a= np.where((a==0)&(np.pad(a,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==Q) , Q, a)
    elif R == '00_1':
        a= np.where((a==0)&(np.pad(a,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==Q), Q, a)
    elif R == '010':
        a= np.where((a==0)&(np.pad(a,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==Q), Q, a)
    elif R == '0_10':
        a= np.where((a==0)&(np.pad(a,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==Q), Q, a)
    

    return a


Selection = [(8, '001'), (8, '00_1'), (8, '010'), (8, '0_10')]

"""
# 1st option ----------------------------------------------------
A = np.zeros((1,5,5)); A[0,2,2]=8
for s in Selection:
    A = G(A, s)
print('Result is OK, but 1st method is NOT:'); print('A=',A); print()
#------------------------------------------------------------------

# 2nd option ----------------------------------------------------
A = np.zeros((1,5,5)); A[0,2,2]=8
GP = partial(G, A)
for s in Selection:
    A = GP(s)
print('Result is wrong and 2nd method, too:'); print('A=',A); print()
#------------------------------------------------------------------

# 3rd option ----------------------------------------------------
A = np.zeros((1,5,5)); A[0,2,2]=8
for s in Selection:
    GP = partial(G, A)
    A = GP(s)
print('Result is OK, but 3rd option is NOT:'); print('A=',A); print()
#------------------------------------------------------------------

# 4th option ----------------------------------------------------
A = np.zeros((1,5,5)); A[0,2,2]=8
GP = partial(G, A)
A = map(GP, Selection)
print('Result is NOT, but the 4th method is OK:'); print('A=',list(A)); print()
#------------------------------------------------------------------
"""

# 5th option is MPI (calculated using multi-cores by multiprocessing module)
A = np.zeros((1,5,5)); A[0,2,2]=8
if __name__ =='__main__':
    pool=Pool()
    manager=Manager()
    A=manager.list(A)
    
    for i in range(3):
        GP = partial(G, A)
        GP_list = pool.map(GP, Selection)
        for j in range(len(GP_list)-1):
            if j==0:
                A=np.where(GP_list[j+1]==0, GP_list[j], GP_list[j+1])
            else:
                A=np.where(GP_list[j+1]==0, A, GP_list[j+1])
    print('Result and method are OK:'); print('A=',A); print()

    pool.close()
    pool.join()
    










    
