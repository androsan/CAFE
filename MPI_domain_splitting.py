from multiprocessing import Pool, Process, Queue, Manager, Array, sharedctypes, set_start_method
from functools import partial
import itertools
import time
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

def merging(prva,druga):
    tot = np.where(druga==0, prva, druga)
    return tot

def GG_pool(attr2, attr3):      #  attr3= faza, attr2= Selection
       grain=attr2[0]; s=attr2[1]
       ''' ----------------------------------------- GROUP 1 ::: [001], [010] >>> 4 sites -------------------------------------------'''
       if s == '001': attr3= np.where((attr3==0)&(np.pad(attr3,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain) , grain, attr3)
       elif s == '00_1': attr3= np.where((attr3==0)&(np.pad(attr3,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain), grain, attr3)
       elif s == '010': attr3= np.where((attr3==0)&(np.pad(attr3,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain), grain, attr3)
       elif s == '0_10': attr3= np.where((attr3==0)&(np.pad(attr3,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain), grain, attr3)
       elif s == '100': attr3= np.where((attr3==0)&(np.pad(attr3,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain), grain, attr3)
       elif s == '_100': attr3= np.where((attr3==0)&(np.pad(attr3,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain), grain, attr3)
       return attr3


def GG_process(Attr2, attr3, q):      #  attr3= faza, attr2= Selection,  q= Queue()
       for attr2 in Attr2:
           grain=attr2[0]; s=attr2[1]
           ''' ----------------------------------------- GROUP 1 ::: [001], [010] >>> 4 sites -------------------------------------------'''
           if s == '001': attr3= np.where((attr3==0)&(np.pad(attr3,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain) , grain, attr3)
           elif s == '00_1': attr3= np.where((attr3==0)&(np.pad(attr3,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain), grain, attr3)
           elif s == '010': attr3= np.where((attr3==0)&(np.pad(attr3,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain), grain, attr3)
           elif s == '0_10': attr3= np.where((attr3==0)&(np.pad(attr3,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain), grain, attr3)
           elif s == '100': attr3= np.where((attr3==0)&(np.pad(attr3,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain), grain, attr3)
           elif s == '_100': attr3= np.where((attr3==0)&(np.pad(attr3,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain), grain, attr3)
       q.put(attr3)
       




Z,X,Y = 1, 501, 1001

A = np.zeros((Z, X, Y))# definition of initial, i.e. basic domain

A[0,250,450] = 1                        # "nucleation", state of cell 1 denotes solid state (grain_ID == 1), whereas 0 is for liquid
A[0,330,660] = 2


grain_ID = np.array([1,2])
smeri=np.array(['001', '00_1', '010', '0_10', '100', '_100',])[:4]
Selection = list((itertools.product(grain_ID, smeri)))

#from Long_track_CA_NEW_GUI import Selection_Mechanism, Random_Grains_n_Directions
#Selection = Selection_Mechanism(grain_ID, smeri[:4], 7)



domain_splitting_sequential = False
domain_splitting_PARALLEL = True

if domain_splitting_sequential:
    # .................. domain splitting .....................................................................
    start_ds_sequential = time.time()
    Ycut_indices = [0, int(Y/2), Y]            # splitting basic matrix into 2 partitions of equal size in Y-direction

    for i in range(20):
        # domain partitioning
        A1=A[:,:,Ycut_indices[0]:Ycut_indices[1]+1]                        # 1 st partition (in Y-direction)
        A2=A[:,:,Ycut_indices[1]-1:Ycut_indices[2]]                         # 2nd partition (in Y-direction)
        # transition function on 1st partition ( CPU1)
        for j in Selection:
            A1=GG_pool(j, A1)
        # transition function on 2nd partition ( CPU2)
        for j in Selection:
            A2=GG_pool(j, A2)
        A=np.dstack((A1[:,:,Ycut_indices[0]:Ycut_indices[1]], A2[:,:,1:]))   #  Stacking of arrays horizontally (in Y-direction)
    end_ds_sequential = time.time()
    print('Domain splitting computing time - sequential = ', round(end_ds_sequential-start_ds_sequential, 3),' sec.')




if domain_splitting_PARALLEL:
    # .................. transformation of numpy arrays into sharedctypes .....................................................................
    start_ds_parallel = time.time()
    Ycut_indices = [0, int(Y/2), Y]            # splitting basic matrix into 2 partitions of equal size in Y-direction

    #A_ctype = np.ctypeslib.as_ctypes(A)
    #A = sharedctypes.Array(A_ctype._type_, A_ctype, lock=False)

    #A1=A[:,:,Ycut_indices[0]:Ycut_indices[1]+1]                        # 1st partition (in Y-direction)
    #A2=A[:,:,Ycut_indices[1]-1:Ycut_indices[2]]                         # 2nd partition (in Y-direction)
    
    #A1=A1.copy() 
    #A1_ctype = np.ctypeslib.as_ctypes(A1)
    #A1 = sharedctypes.Array(A1_ctype._type_, A1_ctype, lock=True)
    
    #A2=A2.copy()
    #A2_ctype = np.ctypeslib.as_ctypes(A2)
    #A2 = sharedctypes.Array(A2_ctype._type_, A2_ctype, lock=True)
    
    if __name__ == '__main__':
        #pool = Pool(processes=2, initializer=initProcess, initargs=(A1, ))
        #pool = Pool(processes=2)
        set_start_method('spawn')
        manager=Manager()
        Q1=manager.Queue()
        Q2=manager.Queue()
        
        #work=[]
        for i in range(20):
            # --- domain partitioning ---
            #A1=np.frombuffer(A1_ctype).reshape((Z,X,int(Y/2)+1))
            #A2=np.frombuffer(A2_ctype).reshape((Z,X,int(Y/2)+2))
            
            A1=A[:,:,Ycut_indices[0]:Ycut_indices[1]+1]                        # 1st partition (in Y-direction)
            A2=A[:,:,Ycut_indices[1]-1:Ycut_indices[2]]                         # 2nd partition (in Y-direction)
            '''
            A1=A1.copy(); A2=A2.copy()
            A1_ctype = np.ctypeslib.as_ctypes(A1)
            A1 = sharedctypes.Array(A1_ctype._type_, A1_ctype, lock=False)
            A2_ctype = np.ctypeslib.as_ctypes(A2)
            A2 = sharedctypes.Array(A2_ctype._type_, A2_ctype, lock=False)
            '''
            '''
            # transition function on 1st partition ( CPU1)
            for j in Selection:
                A1=GG_pool(j, A1)
            # transition function on 2nd partition ( CPU2)
            for j in Selection:
                A2=GG_pool(j, A2)
            '''

            #for j in Selection:
             #   A_list = pool.map(partial(GG_pool,j), [A1,A2] )
              #  A1,A2 = A_list[0], A_list[1]

            #work.append(pool.apply_async(GG_process,(Selection, A1)))
            
            p1 = Process(target=GG_process, args=(Selection, A1,Q1))
            p2 = Process(target=GG_process, args=(Selection, A2,Q2))

            p1.start()
            p2.start()
            
            A1=Q1.get()
            A2=Q2.get()
            
            p1.join()
            p2.join()
            
            A=np.dstack((A1[:,:,Ycut_indices[0]:Ycut_indices[1]], A2[:,:,1:]))   #  Stacking of arrays horizontally (in Y-direction)


        end_ds_parallel = time.time()
        print('Domain splitting computing time - parallel = ', round(end_ds_parallel-start_ds_parallel, 3),' sec.')
    
    

    #A=np.vstack((A1[:2],A2[1:]))   #  Stacking of arrays vertically (in Z-direction)
    '''
    ds_start = time.time()
    if __name__ == '__main__':
        pool = Pool()
        for i in range(2):
            b=d[:3]; c=d[1:]
            for j in Selection:
                faza_list = pool.map(partial(GG_pool,j), [b,c] )
                b,c=faza_list[0], faza_list[1]
            d=np.vstack((b[:2],c[1:]))   
        ds_stop=time.time()
        print('MPI domain splitting computing time = ', int(ds_stop-ds_start),' sec.')
        pool.close()
        pool.join()
    '''















single_core = False
MPI_Pool = False
MPI_Process = False

if single_core:
    single_start = time.time()
    for i in range(100):
        for j in Selection:
            A=GG_pool(j, A)
    single_stop = time.time()
    print('Single Core computing time = ', round(single_stop-single_start, 3),' sec.')





if MPI_Pool:
    mpi_start = time.time()
    if __name__ == '__main__':
        pool = Pool()    
        for i in range(1):
            od=0; do=4
            for j in range(int(len(Selection)/4)):
                cpu_list = pool.map(partial(GG_pool, faza), Selection[od:do])
                total=merging(cpu_list[0], cpu_list[1])
                if len(cpu_list)> 2:
                    for move in range(2, len(cpu_list)):
                        total=merging(total, cpu_list[move])
                faza = total.copy()
                od+=4; do+=4
        pool.close()
        pool.join()

    mpi_stop = time.time()
    print('MPI computing time = ', int(mpi_stop-mpi_start))





if MPI_Process:
    mpi_start = time.time()
    if __name__ == '__main__':
        
        for i in range(2):
            procesi=[]
            for j in Selection:
                p = Process(target=GG_pool, args=[faza,j])
                p.start()
                procesi.append(p)

            for q in procesi:
                q.join()
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


    





        

