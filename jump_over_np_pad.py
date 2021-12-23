import random
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

def time_counter(neg_index):
   Negatives[neg_index] += dt

def random_selection(x):
    global item
    if len(x)> 0:
        item = random.choice(x)
        x=np.delete(x,np.where(x==item))
    else:
        item = x[0]
        x=np.delete(x,np.where(x==item))
    return x

def merging(prva,druga):
    tot = np.where(druga==0, prva, druga)
    return tot

def W(O, L):
    O = np.array(O)
    cos_Z = np.dot(O[0], L) / (np.linalg.norm(O[0]* np.linalg.norm(L)))
    cos_X = np.dot(O[1], L) / (np.linalg.norm(O[1] * np.linalg.norm(L)))
    cos_Y = np.dot(O[2], L) / (np.linalg.norm(O[2] * np.linalg.norm(L)))
    return np.max(np.absolute(np.array([cos_Z, cos_X, cos_Y])))


'''[001]'''
def pad001(zrnec):pad001=(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==zrnec);return pad001
def pad00_1(zrnec):pad00_1=(np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==zrnec); return pad00_1

'''[010]'''
def pad010(zrnec):pad010=(np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==zrnec);return pad010
def pad0_10(zrnec):pad0_10=(np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==zrnec);return pad0_10

'''[011]'''
def pad011(zrnec):pad011=(np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==zrnec); return pad011
def pad01_1(zrnec):pad01_1=(np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==zrnec);return pad01_1
def pad0_11(zrnec):pad0_11=(np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==zrnec); return pad0_11
def pad0_1_1(zrnec):pad0_1_1=(np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==zrnec); return pad0_1_1

'''[002]'''
def pad002(zrnec):pad002=(np.pad(faza,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==zrnec); return pad002
def pad00_2(zrnec):pad00_2=(np.pad(faza,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==zrnec); return pad00_2

'''[020]'''
def pad020(zrnec):pad020=(np.pad(faza,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==zrnec); return pad020
def pad0_20(zrnec):pad0_20=(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==zrnec); return pad0_20

'''[012]'''
def pad012(zrnec):pad012=(np.pad(faza,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==zrnec); return pad012
def pad01_2(zrnec):pad01_2=(np.pad(faza,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==zrnec);return pad01_2
def pad0_12(zrnec):pad0_12=(np.pad(faza,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==zrnec);return pad0_12
def pad0_1_2(zrnec):pad0_1_2=(np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==zrnec);return pad0_1_2

'''[021]'''
def pad021(zrnec):pad021=(np.pad(faza,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==zrnec);return pad021
def pad02_1(zrnec):pad02_1=(np.pad(faza,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==zrnec);return pad02_1
def pad0_21(zrnec):pad0_21=(np.pad(faza,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==zrnec);return pad0_21
def pad0_2_1(zrnec):pad0_2_1=(np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==zrnec);return pad0_2_1

'''[022]'''
def pad022(zrnec):pad022=(np.pad(faza,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==zrnec);return pad022
def pad02_2(zrnec):pad02_2=(np.pad(faza,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==zrnec);return pad02_2
def pad0_22(zrnec):pad0_22=(np.pad(faza,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==zrnec);return pad0_22
def pad0_2_2(zrnec):pad0_2_2=(np.pad(faza,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==zrnec);return pad0_2_2

'''[030]'''
def pad030(zrnec):pad030=(np.pad(faza,((0,0),(3,0), (0,0)), 'constant')[:,:-3,:]==zrnec);return pad030
def pad0_30(zrnec):pad0_30=(np.pad(faza,((0,0),(0,3), (0,0)), 'constant')[:,3:,:]==zrnec);return pad0_30

'''[013]'''
def pad013(zrnec):pad013=(np.pad(faza,((0,0),(1,0), (3,0)), 'constant')[:,:-1,:-3]==zrnec);return pad013
def pad01_3(zrnec):pad01_3=(np.pad(faza,((0,0),(1,0), (0,3)), 'constant')[:,:-1,3:]==zrnec);return pad01_3
def pad0_13(zrnec):pad0_13=(np.pad(faza,((0,0),(0,1), (3,0)), 'constant')[:,1:,:-3]==zrnec);return pad0_13
def pad0_1_3(zrnec):pad0_1_3=(np.pad(faza,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==zrnec);return pad0_1_3

'''[031]'''
def pad031(zrnec):pad031=(np.pad(faza,((0,0),(3,0), (1,0)), 'constant')[:,:-3,:-1]==zrnec);return pad031
def pad03_1(zrnec):pad03_1=(np.pad(faza,((0,0),(3,0), (0,1)), 'constant')[:,:-3,1:]==zrnec);return pad03_1
def pad0_31(zrnec):pad0_31=(np.pad(faza,((0,0),(0,3), (1,0)), 'constant')[:,3:,:-1]==zrnec);return pad0_31
def pad0_3_1(zrnec):pad0_3_1=(np.pad(faza,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==zrnec);return pad0_3_1


def Dsr_1st(r):
   if random_distances:
      cell_1st = cell * (0.5 + 1.158*r)                              #"""------------------1st shell ::: [001], [010], [100] ::: 6 neighbours ------------------"""
   else:
      cell_1st=cell
   return cell_1st

def Dsr_5th(r):
   if random_distances:
      cell_5th=cell*(1.5 + 1.098*r)                                    #''' ------------------5th shell ::: [002], [020], [200] ::: 6 neighbours ------------------ '''
   else:
      cell_5th=cell*2
   return cell_5th


Z=1
X=11
Y=11

faza = np.zeros((Z,X,Y))
taula = np.zeros((Z,X,Y))
likvid = np.zeros((Z,X,Y))
cas = np.zeros((Z,X,Y))                                # časovna matrika
vg = np.zeros((Z,X,Y))                                  # matrika hitrosti rasti
vg =   1                                                               # Value of homogeneous 'growing velocity' field, a.u., for testing and code development


random_distances = False; R=None

cell=1                       # for mesh dependency development (MDD)
dt=cell/6

Negatives = {-1:0}; asc = {}; grain_counter = 0; S = {}; np.random.seed(75689803)
M = { 1: {'ß':(0, 5, 5)},
           2: {'ß':(0, 6, 6)},
             }
for i in M:
    faza[M[i]['ß'][0],M[i]['ß'][1],M[i]['ß'][2]]=i                               # define nucleus ID in faza matrix
    asc[i] = str(i)

grain_ID = np.array(list(asc.keys()))
cas[np.isin(faza, grain_ID, invert=False)] = -1
taula=0;likvid=0; grain_counter=len(grain_ID)



smeri = np.array(['001', '00_1', '010', '0_10', '002', '00_2', '020', '0_20', ])
#smeri = np.array(['001', '002'])

'''**************************'''
START_step =         0                       # Starting time step number
END_step     =         14                     # Ending time step number
'''**************************'''

for i in range(START_step, END_step):
    for negind in Negatives:
        time_counter(negind)

    for s in smeri[:]:
        for negind in Negatives:
            ''' ----------------------------------------- 1st shell -------------------> 6 nearest neighbours '''
            if s == '001':
                if negind == -1 and i==START_step:
                    cas001 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==negind), Negatives[negind], cas)
                    S[s]=cas001; del cas001
                else:
                    cas001 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==negind), Negatives[negind], S[s])
                    S[s]=cas001; del cas001

            elif s == '00_1':
                if negind == -1 and i==START_step:
                    cas00_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==negind), Negatives[negind], cas)
                    S[s]=cas00_1; del cas00_1
                else:
                    cas00_1 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==negind), Negatives[negind], S[s])
                    S[s]=cas00_1; del cas00_1

            elif s == '010':
                if negind == -1 and i==START_step:
                    cas010 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==negind), Negatives[negind], cas)
                    S[s]=cas010; del cas010
                else:
                    cas010 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==negind), Negatives[negind], S[s])
                    S[s]=cas010; del cas010

            elif s == '0_10':
                if negind == -1 and i==START_step:
                    cas0_10 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==negind), Negatives[negind], cas)
                    S[s]=cas0_10; del cas0_10
                else:
                    cas0_10 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==negind), Negatives[negind], S[s])
                    S[s]=cas0_10; del cas0_10

            elif s == '002':
                if negind == -1 and i==START_step:
                    cas002 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==negind), Negatives[negind], cas)
                    S[s]=cas002; del cas002
                else:
                    cas002 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==negind), Negatives[negind], S[s])
                    S[s]=cas002; del cas002

            elif s == '00_2':
                if negind == -1 and i==START_step:
                    cas00_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==negind), Negatives[negind], cas)
                    S[s]=cas00_2; del cas00_2
                else:
                    cas00_2 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==negind), Negatives[negind], S[s])
                    S[s]=cas00_2; del cas00_2

            elif s == '020':
                if negind == -1 and i==START_step:
                    cas020 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==negind), Negatives[negind], cas)
                    S[s]=cas020; del cas020
                else:
                    cas020 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==negind), Negatives[negind], S[s])
                    S[s]=cas020; del cas020

            elif s == '0_20':
                if negind == -1 and i==START_step:
                    cas0_20 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(cas,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==negind), Negatives[negind], cas)
                    S[s]=cas0_20; del cas0_20
                else:
                    cas0_20 = np.where((np.isin(faza, grain_ID, invert=True)) & (np.pad(S[s],((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==negind), Negatives[negind], S[s])
                    S[s]=cas0_20; del cas0_20            


    F = faza.copy()
    for g in grain_ID:
           grain_ID = random_selection(grain_ID)
           grain = item
           directions = smeri.copy()
           seznam_premikov = []
           for d in directions:
               directions = random_selection(directions)
               s = item
               
               ''' ----------------------------------------- 1st shell -------------------> 6 nearest neighbours '''
               if s == '001':
                   dij = np.array([0,0,1])
                   lij = vg*S[s]
                   faza001 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain), grain, faza)
                   seznam_premikov.append(faza001); del faza001

               elif s == '00_1':
                   dij = np.array([0,0,-1])
                   lij = vg*S[s]
                   faza00_1 = np.where((faza==0)& (taula==0)& (likvid==0)& (lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain), grain, faza)
                   seznam_premikov.append(faza00_1); del faza00_1

               elif s == '010':
                   dij = np.array([0,1,0])
                   lij = vg*S[s]
                   faza010 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain), grain, faza)
                   seznam_premikov.append(faza010); del faza010

               elif s == '0_10':
                   dij = np.array([0,-1,0])
                   lij = vg*S[s]
                   faza0_10 = np.where((faza==0)& (taula==0)& (likvid==0)&(lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain), grain, faza)
                   seznam_premikov.append(faza0_10); del faza0_10

               elif s == '002':
                   dij = np.array([0,0,2])
                   lij = vg*S[s]
                   faza002 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )
                                                   &(np.pad(faza,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)
                                                   & ((np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)                  # 001
                                                          | (np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)             # 011
                                                          |  (np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)             # 0_11
                                                          | (np.pad(faza,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)             # 012           
                                                          |  (np.pad(faza,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain))          # 0_12
                                                , grain, faza) 
                   seznam_premikov.append(faza002); del faza002

               elif s == '00_2':
                   dij = np.array([0,0,-2])
                   lij = vg*S[s]
                   faza00_2 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )
                                                      &(np.pad(faza,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)
                                                      & ((np.pad(faza,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)                # 00_1
                                                             | (np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)          # 01_1
                                                             |  (np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)          # 0_1_1
                                                             | (np.pad(faza,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)          # 01_2
                                                             |  (np.pad(faza,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain))       # 0_1_2
                                                     , grain, faza)
                   seznam_premikov.append(faza00_2); del faza00_2

               elif s == '020':
                   dij = np.array([0,2,0])
                   lij = vg*S[s]
                   faza020 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )
                                                   &(np.pad(faza,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)
                                                   & ((np.pad(faza,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)                  # 010 
                                                          | (np.pad(faza,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)             # 011
                                                          |  (np.pad(faza,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)             # 01_1
                                                          | (np.pad(faza,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)             # 021
                                                          |  (np.pad(faza,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain))          # 02_1  
                                                , grain, faza)
                   seznam_premikov.append(faza020); del faza020

               elif s == '0_20':
                   dij = np.array([0,-2,0])
                   lij = vg*S[s]
                   faza0_20 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )
                                                      &(np.pad(faza,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)
                                                      & ((np.pad(faza,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)               # 0_10
                                                             | (np.pad(faza,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)         # 0_11
                                                             |  (np.pad(faza,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)         # 0_1_1
                                                             | (np.pad(faza,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)         # 0_21
                                                             |  (np.pad(faza,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain))      # 0_2_1
                                                      , grain, faza)     
                   seznam_premikov.append(faza0_20); del faza0_20


           total=merging(seznam_premikov[0], seznam_premikov[1])

           if len(smeri)> 2:
               for move in range(2, len(smeri)):
                   total=merging(total, seznam_premikov[move])
           
           faza = total.copy()

    grain_ID = np.array(list(asc.keys()))
    
    if not np.all(faza==F):
        negind -=1
        for s in S:
            S[s][np.isin(faza-F, grain_ID)]= negind
        Negatives[negind]=0
    
    times = np.array(list(S.values()))
    cas_total = merging(times[0], times[1])
    if len(smeri)> 2:
        for move in range(2, len(smeri)):
            cas_total = merging(cas_total, times[move])
    
    try:
        cas  = cas_total.copy()
    except NameError:
        pass
    



plt.imshow(faza[0])
plt.figure(); plt.imshow(cas[0])



"""
'''........................................ [012], [021] families ::: Dsr_4th(R) ::: .........................................  8  neighbours '''                 

faza012 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )& pad012(grain)
                                &(pad001 (grain)| pad011 (grain)| pad002 (grain)| pad021 (grain)| pad022(grain)), grain, faza)

faza01_2 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )& pad01_2(grain)
                                   & (pad00_1(grain) | pad01_1(grain) | pad00_2(grain) | pad02_1(grain) | pad02_2(grain)), grain, faza)

faza0_12 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )& pad0_12(grain)
                                   & (pad001(grain) | pad0_11(grain) | pad002(grain) | pad0_21(grain) | pad0_22(grain)), grain, faza)

faza0_1_2 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )& pad0_1_2(grain)
                                      & (pad00_1(grain) | pad0_1_1(grain) | pad00_2(grain) | pad0_2_1(grain) | pad0_2_2(grain)), grain, faza)


faza021 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )& pad021(grain)
                                & (pad010(grain) | pad011(grain) | pad020(grain) | pad030(grain) | pad012(grain)), grain, faza)

faza02_1 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )& pad02_1(grain)
                                   & (pad010(grain) | pad01_1(grain) | pad020(grain) | pad030(grain) | pad01_2(grain)), grain, faza)

faza0_21 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )& pad0_21(grain)
                                   & (pad0_10(grain) | pad0_11(grain) | pad0_20(grain) | pad0_30(grain) | pad0_12(grain)), grain, faza)

faza0_2_1 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_4th(R) )& pad0_2_1(grain)
                                      & (pad0_10(grain) | pad0_1_1(grain) | pad0_20(grain) | pad0_30(grain) | pad0_1_2(grain)), grain, faza)

''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''


'''........................................ [002], [020] families ::: Dsr_5th(R) ::: .........................................  4  neighbours '''                 

faza002 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )& pad002(grain)
                                & (pad001(grain) | pad011(grain) | pad0_11(grain) | pad012(grain) | pad0_12(grain)), grain, faza)
              
faza00_2 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )& pad00_2(grain)
                                   & (pad00_1(grain) | pad01_1(grain) | pad0_1_1(grain) | pad01_2(grain) | pad0_1_2(grain)), grain, faza)

faza020 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )& pad020(grain)
                                & (pad010(grain) | pad011(grain) | pad01_1(grain) | pad021(grain) | pad02_1(grain)), grain, faza)

faza0_20 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_5th(R) )& pad0_20(grain)
                                   & (pad0_10(grain) | pad0_11(grain) | pad0_1_1(grain) | pad0_21(grain) | pad0_2_1(grain)), grain, faza) )
                    
'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''


'''................................................ [022] family ::: Dsr_6th(R) ::: ....................................................  4  neighbours '''    
                
faza022 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )& pad022(grain)
                                & (pad012(grain) | pad011(grain) | pad021(grain) | pad013(grain) | pad031(grain)), grain, faza)

faza02_2 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )& pad02_2(grain)
                                   & (pad01_2(grain) | pad01_1(grain) | pad02_1(grain) | pad01_3(grain) | pad03_1(grain)) , grain, faza)

faza0_22 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )& pad0_22(grain)
                                   & (pad0_12(grain) | pad0_11(grain) | pad0_21(grain) | pad0_13(grain) | pad0_31(grain)), grain, faza)

faza0_2_2 = np.where((faza==0)& (taula==0)& (likvid==0) &(lij>=Dsr_6th(R) )& pad0_2_2(grain)
                                      & (pad0_1_2(grain) | pad0_1_1(grain) | pad0_2_1(grain) | pad0_1_3(grain) | pad0_3_1(grain)), grain, faza)
                     
'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

"""







