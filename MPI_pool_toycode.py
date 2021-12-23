from multiprocessing import Pool
from functools import partial
import itertools
from inspect import signature
import numpy as np
import time

# to obtain arguments of callable:
# signature(p).parameters    # where p is callable


H=np.full((3,100,100,),1)
G=np.full((3,100,100,),2)
counter=0

def g_pool(y):
    global H,G, counter
    H=(H+y[0])
    G=(y[1]+H)
    counter+=1
    return H,G

def g_exe(Y, H, G):
    global counter
    for y in Y:
        H=(H+y[0])
        G=(y[1]+H)
        counter+=1
    return H,G


a = range(1,1000)
b = range(1,1000)

c = list(itertools.product(a, b))                  # it replaces double (nested) for loops, i.e. for i in a: for j in b


SINGLE_executor = False            # CPU: 39.4 s    |     H[0,0,0] =  499000501     |     G[0,0,0] =  499001500     |     counter =  998001
SINGLE_pool =        False           # CPU: 40.7 s    |     H[0,0,0] =  499000501     |     G[0,0,0] =  499001500     |     counter =  998001

TPE =       False                            # CPU: 39.4 s    |     H[0,0,0] =  499000501     |     G[0,0,0] =  499001500     |     counter =  998001
POOL =     True



if SINGLE_executor:
    start = time.time()
    #----------------------------------------------------------------------------------------------- process
    H,G = g_exe(c, H, G)
    #----------------------------------------------------------------------------------------------------------
    end = time.time()
    print('Single executor process FINISHED in   ', round(end-start, 1), ' seconds.') 
    print('H[0,0,0] = ', H[0,0,0])
    print('G[0,0,0] = ', G[0,0,0])
    print('counter = ',counter)

    
if SINGLE_pool:
    start = time.time()
    #----------------------------------------------------------------------------------------------- process
    #L=[]
    for i in c:
        H,G = g_pool(i)
        #L.append((H,G))
    #----------------------------------------------------------------------------------------------------------
    end = time.time()
    print('Single pool process FINISHED in   ', round(end-start, 1), ' seconds.') 
    print('H[0,0,0] = ', H[0,0,0])
    print('G[0,0,0] = ', G[0,0,0])
    print('counter = ',counter)
    

if True:#POOL:
    start = time.time()
    #----------------------------------------------------------------------------------------------- process
    if __name__ == '__main__':
        pool = Pool()
        #basic_partial = partial( basic, H,G)
        #w=pool.map(basic_partial, c)

        #g_part = partial(g, HH, GG)
        #w=pool.map(g_part, c)
        print('bua bana')
        w=pool.map(g_pool, c)
        pool.close()
        pool.join()
    #----------------------------------------------------------------------------------------------------------
    end = time.time()
    print('Pool process FINISHED in   ', round(end-start, 1), ' seconds.') 
    #print('H[0,0,0] = ', w[0][0][0,0,0])
    #print('G[0,0,0] = ', w[0][1][0,0,0])
    print('counter = ',counter)



if TPE:
    import multiprocessing
    import concurrent.futures
    start = time.time()
    #----------------------------------------------------------------------------------------------- process
    if __name__ == '__main__':
        with concurrent.futures.ThreadPoolExecutor() as executor:
            r=[]
            r.append(executor.submit(g_exe, c, H, G))
                    
    q=[i.result() for i in r]
    #----------------------------------------------------------------------------------------------------------
    end = time.time()
    print('ThreadPoolExecutor process FINISHED in   ', round(end-start, 1), ' seconds.') 
    print('H[0,0,0] = ', q[0][0][0,0,0])
    print('G[0,0,0] = ', q[0][1][0,0,0])
    print('counter = ',counter)





