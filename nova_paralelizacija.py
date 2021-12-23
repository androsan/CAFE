from multiprocessing import Pool, Process, freeze_support, set_start_method
from functools import partial
import itertools
import time
import numpy as np
import matplotlib.pyplot as plt
import os
plt.ion()
def merging(prva,druga):
    tot = np.where(druga==0, prva, druga)
    return tot
def GG_(attr3, attr2):      #  attr3= faza, attr2= Selection
       grain=attr2[0];s=attr2[1]
       ''' ----------------------------------------- GROUP 1 ::: [001], [010] >>> 4 sites -------------------------------------------'''
       if s ==    '001':     attr3= np.where((attr3==0)&(np.pad(attr3,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain) , grain, attr3)
       elif s == '00_1':  attr3= np.where((attr3==0)&(np.pad(attr3,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain), grain, attr3)
       elif s == '010':     attr3= np.where((attr3==0)&(np.pad(attr3,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain), grain, attr3)
       elif s == '0_10':  attr3= np.where((attr3==0)&(np.pad(attr3,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain), grain, attr3)
       return attr3


from multiprocessing import sharedctypes

Z,X,Y = 1, 1001, 1001
faza = np.zeros((Z,X,Y))

#faza0 = np.ctypeslib.as_ctypes(np.zeros((Z,X,Y)))
#faza1 = sharedctypes.Array(faza0._type_, faza0, lock=True)
#faza1 = sharedctypes.copy(faza0)
#faza = np.ctypeslib.as_array(faza1)


M = { 1: {'ß':(0, 500, 500)}, 
           2: {'ß':(0, 800, 500)}, 
           3: {'ß':(0, 350, 300)},                  
           4: {'ß':(0, 200, 100)},
           5: {'ß':(0, 750, 750)}, 
        }
for i in M:
    faza[M[i]['ß'][0],M[i]['ß'][1],M[i]['ß'][2]]=i  

grain_ID = np.array([1,2,3,4,5])
smeri=np.array(['001', '00_1', '010', '0_10'])
Selection = list((itertools.product(grain_ID, smeri)))


def single_core():
    global faza
    single_start = time.time()
    for i in range(20):
        for j in Selection:
            faza=GG_(faza,j)
    single_stop = time.time()
    print('Single Core computing time = ', int(single_stop-single_start))


def MPI_Pool():
    mpi_start = time.time()
    os.system("taskset -p 0xff %d" % os.getpid())
    #freeze_support()
    #set_start_method('spawn')
    pool = Pool()
    for i in range(20):
        od=0; do=16
        for j in range(int(len(Selection)/16)):
            
            cpu_list = pool.map(partial(GG_, faza), Selection[od:do])

            '''
            total=merging(cpu_list[0], cpu_list[1])
            if len(cpu_list)> 2:
                for move in range(2, len(cpu_list)):
                    total=merging(total, cpu_list[move])
            faza = total.copy()
            od+=16; do+=16
            '''
    pool.close()
    pool.join()

    mpi_stop = time.time()
    print('MPI computing time = ', int(mpi_stop-mpi_start))


def MPI_Process():
    mpi_start = time.time()
    freeze_support()
    for i in range(2):
        procesi=[]
        for j in Selection:
            p = Process(target=GG_, args=[faza,j])
            p.start()
            procesi.append(p)

        for q in procesi:
            q.join()


import concurrent.futures

def MPI_Concurrent_Futures():
    global faza
    mpi_start = time.time()
    '''
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for i in range(20):
            for j in Selection[:]:
                results=[]
                results.append(executor.submit(partial(GG_, faza),j))
                faza=results[0].result()
    '''
    with concurrent.futures.ProcessPoolExecutor() as executor:
        Q=partial(GG_, faza)
        for i in range(20):
            for j, faza in zip(Selection[:], executor.map(Q, Selection[:])):
                Q=partial(GG_, faza)
    

    '''
        #for i in range(1):
            #for j, faza in zip(range(int(len(Selection))), executor.map(partial(GG_, faza), Selection) ):
                #print(j)
                #executor.map(partial(GG_, faza), Selection)

                
                total=merging(cpu_list[0], cpu_list[1])
                if len(cpu_list)> 2:
                    for move in range(2, len(cpu_list)):
                        total=merging(total, cpu_list[move])
                faza = total.copy()
                od+=16; do+=16
                '''
    mpi_stop = time.time()
    print('MPI computing time = ', int(mpi_stop-mpi_start))
    


if __name__ == '__main__':
    #single_core()
    MPI_Concurrent_Futures()
    









'''
# The tutorial from StackOverFlow, how Python for loop is paralelized
# setup output lists
output1 = list()
output2 = list()
output3 = list()

for j in range(0, 10):
    # calc individual parameter value
    parameter = j * offset
    # call the calculation
    out1, out2, out3 = calc_stuff(parameter = parameter)

    # put results into correct output list
    output1.append(out1)
    output2.append(out2)
    output3.append(out3)


pool = multiprocessing.Pool(4)
out1, out2, out3 = zip(*pool.map(calc_stuff, range(0, 10 * offset, offset)))
'''


    





        

