from multiprocessing import Pool, Manager, freeze_support
from functools import partial
import numpy as np
import itertools
import matplotlib.pyplot as plt
plt.ion()


#--------------------------------------- GUI Tkinter part ----------------------------------------------------------------------------------
from tkinter import *

root=Tk();
var0 = StringVar();
lab0 = Label(root, textvariable=var0, font=('Arial', '22')); lab0.grid(row=100, column=100)

def stop_process():
    root.after_cancel(pol)
#--------------------------------------- GUI Tkinter part ----------------------------------------------------------------------------------
run_1st = True
run_2nd = True
run_3rd = True

MPI = True

#  ~~~~  Manual Nucleation (MN)~~~~~~
Z,X,Y = 1, 11, 11
faza = np.zeros((Z,X,Y)) 
cas = np.zeros((Z,X,Y))
dt = 1

Negatives = {-1:0}; S = {}

M = { 1: {'ß':(0, 3, 2),},
          3: {'ß':(0, 5, 6),},

      }

for i in M:
    faza[M[i]['ß'][0],M[i]['ß'][1],M[i]['ß'][2]]=i                               # define nucleus ID in faza matrix
    
grain_ID = np.array([1, 3])
cas[np.isin(faza, grain_ID, invert=False)] = -1



''' ......................................... time counter PART :: function. . . '''
def time_counter(attr0, neg_index):   #  attr0= Negatives
    attr0[neg_index] += dt
    return attr0

def time_counter_original(neg_index):
    global Negatives
    Negatives[neg_index] += dt
    return Negatives




#Selection=[(1, '0_10'), (1, '01_1'), (1, '0_1_1'), (1, '001'), (1, '0_11'), (1, '010'), (1, '011'), (1, '00_1')]
Selection=[(1, '001'), (1, '00_1'), (1, '010'), (1, '0_10'),
                (3, '001'), (3, '00_1'), (3, '010'), (3, '0_10'),
           ]
#Selection=[(1, '001'), (1, '00_1')]

#smeri=np.array(['001', '00_1', '010', '0_10', '011', '01_1', '0_11', '0_1_1'], dtype='<U5')
smeri=np.array(['001', '00_1', '010', '0_10'], dtype='<U4')
#smeri=np.array(['001', '00_1'], dtype='<U4')

i=-1
END_step=100

import MPI_pool_functions as mpf       # MPI pool functions
def Run_CA(p):
    global i, Negatives, TC_list, S, MM_list, pol, faza, F, cas, GG_list

    i+=1
    
    ''' ......................................... 1st RUN --- time counter PART :: execution . . . '''
    if run_1st:
        if not MPI:
            for negind in Negatives:
                time_counter_original(negind)
                
        elif MPI:
            nega = sorted(list(Negatives.keys()))
            TC_list = p.map(partial(time_counter, Negatives), nega)   # --------> returns :: Negatives
            Negatives = TC_list[-1]
            print('Negatives after time_counter: ', Negatives)



    ''' ......................................... 2nd RUN --- MM_ multiprocessing pool function PART :: execution . . . '''
    if run_2nd:
        if not MPI:
            q=list(itertools.product(smeri, Negatives))
            MM_partial = partial( mpf.MM, Negatives, faza, grain_ID, cas, S)
            for qq in q:
                MM_list = [MM_partial(qq)]
                S = MM_list[-1]
                
        elif MPI:
            nega = sorted(list(Negatives.keys()))
            q=list(itertools.product(smeri, nega))
            MM_list = p.map(partial( mpf.MM, Negatives, faza, grain_ID, cas, S), q)   # --------> returns :: S
            #print('S[001]: ',MM_list[-1]['001'][0,2,3])
            #print('S[00_1]: ',MM_list[-1]['00_1'][0,2,1])
            #S = MM_list[-1]
            for j in MM_list:
                S.update(j)

    ''' ......................................... 3rd RUN --- GG_ multiprocessing pool function PART :: execution . . . '''
    if run_3rd:
        if not MPI:
            F = faza.copy()
            #FF.append(faza)

            # Single-Core using Partial function
            GG_partial=partial( mpf.GG_, faza)
            for selsel in Selection:
                #GG_partial=partial( mpf.GG_, faza, cas, S)
                #GG_partial=partial( mpf.GG_, faza)
                faza = GG_partial(selsel)
                
                #cas = GG_list[1]
                #S = GG_list[2]
            if not np.all(faza==F):
                print('NOVA')
                negind = min(list(Negatives.keys()))- 1             #  Here, negind should be defined as absolute counter of negative numbers to avoid erasing existing Negatives items
                for s in S:
                    S[s][np.isin(faza-F, grain_ID)]= negind
                Negatives[negind]=0
            print(faza)
                
            '''
            # Single-Core without Partial function
            GG_=mpf.GG_
            for selsel in Selection:
                GG_list = GG_(faza,cas,S, selsel)
                faza=GG_list[0]
                cas=GG_list[1]
                S=GG_list[2]
                if not np.all(faza==F):
                    print('NOVA')
                    negind = min(list(Negatives.keys()))- 1             #  Here, negind should be defined as absolute counter of negative numbers to avoid erasing existing Negatives items
                    for s in S:
                        S[s][np.isin(faza-F, grain_ID)]= negind
                    Negatives[negind]=0
                print(faza)
            '''
            
        elif MPI:
            F=faza
            GG_list=p.map(partial( mpf.GG_, faza), Selection)   # --------> returns :: (faza)
            for gp in range(len(GG_list)-1):
                if gp==0:
                    faza=np.where(GG_list[gp+1]==0, GG_list[gp], GG_list[gp+1])
                else:
                    faza=np.where(GG_list[gp+1]==0, faza, GG_list[gp+1])


            #cas=np.where(S['00_1']==0, S['001'], S['00_1'])
            cas=S['001']
            
            #l.acquire()
            #A[np.array([[[True,False,True]]])]=7
            #A=np.where(A==7,np.array([[[True,False,True]]]),7)
            #A=np.where(A==0,np.array([[[0,9,0]]]),7)
            #l.release()
            #cas=GG_list[0][1]
            #S=GG_list[-1][2]
            
            if not np.all(faza==F):
                print('faza and F are not equal')
                #negind -=1
                negind = min(list(Negatives.keys()))- 1             #  Here, negind should be defined as absolute counter of negative numbers to avoid erasing existing Negatives items
                #for s in S._getvalue():
                for s in S:
                    S[s][np.isin(faza-F, grain_ID)]= negind
                Negatives[negind]=0
            

    if i>=END_step:
        root.after_cancel(pol)
    elif i<END_step:
        pol=root.after(200, lambda:Run_CA(p))
    
    

if __name__ =='__main__':
    if MPI:
        fs=freeze_support()
        pool=Pool()
        manager=Manager()
        Negatives=manager.dict(Negatives)
        #S=manager.dict(S)
        faza = manager.list(faza)
        cas = manager.list(cas)
  
        
    elif not MPI:
        pool=None

    Run_CA(pool)
        
    #pool.close()
    #pool.join()


#--------------------------------------- GUI Tkinter part -----------------------------------------------------------------------------------------------------------------------------------
but0 = Button(root, text='Run Process', command=lambda:Run_CA(pool), font=('Arial', '22'), fg='green'); but0.grid(row=101, column=100)
but1 = Button(root, text='Stop Process', command=stop_process, font=('Arial', '22'), fg='red'); but1.grid(row=102, column=100)
#--------------------------------------- GUI Tkinter part -----------------------------------------------------------------------------------------------------------------------------------

def show_S():
    plt.figure()
    plt.title('MM_list[0][0][0] --- S 1st item')
    plt.imshow(MM_list[0][0][0])

    plt.figure()
    plt.title('GG_list[1][0][0] --- faza 2nd item')
    plt.imshow(GG_list[1][0][0])

    plt.figure()
    plt.title('GG_list[2][0][0] --- faza 3rd item')
    plt.imshow(GG_list[2][0][0])

    plt.figure()
    plt.title('GG_list[3][0][0] --- faza 4th item')
    plt.imshow(GG_list[3][0][0])


def show_faza():
    plt.figure()
    plt.title('GG_list[0][0][0] --- faza 1st item')
    plt.imshow(GG_list[0][0][0])

    plt.figure()
    plt.title('GG_list[1][0][0] --- faza 2nd item')
    plt.imshow(GG_list[1][0][0])

    plt.figure()
    plt.title('GG_list[2][0][0] --- faza 3rd item')
    plt.imshow(GG_list[2][0][0])

    plt.figure()
    plt.title('GG_list[3][0][0] --- faza 4th item')
    plt.imshow(GG_list[3][0][0])

def show_cas():    
    plt.figure()
    plt.title('GG_list[0][1][0] --- cas 1st item')
    plt.imshow(GG_list[0][1][0])

    plt.figure()
    plt.title('GG_list[1][1][0] --- cas 2nd item')
    plt.imshow(GG_list[1][1][0])

    plt.figure()
    plt.title('GG_list[2][1][0] --- cas 3rd item')
    plt.imshow(GG_list[2][1][0])

    plt.figure()
    plt.title('GG_list[3][1][0] --- cas 4th item')
    plt.imshow(GG_list[3][1][0])
    




    
