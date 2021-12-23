import numpy as np
import math

def Dsr_1st(r, cell, random_distances):
   if random_distances:
      cell_1st = cell * (0.5 + 1.158312*r)                              #"""------------------1st shell ::: [001], [010], [100] ::: 6 neighbours ------------------"""
   else:
      cell_1st=cell
   return cell_1st

def Dsr_2nd(r, cell, random_distances):
   if random_distances:
      #cell_2nd=cell*math.sqrt(2)*(0.5 + 1.081*r)           #''' ------------------2nd shell ::: [011], [101], [110] ::: 12 neighbours ------------------ '''
      cell_2nd=cell*(0.5*math.sqrt(2)+1.4723*r)
   else:
      cell_2nd=cell*math.sqrt(2)
   return cell_2nd

def Dsr_3rd(r, cell, random_distances):
   if random_distances:
      cell_3rd= None         # !!! to be constructed                #''' ------------------3rd shell ::: [111] ::: 8 neighbours --------------------------- '''
   else:
      cell_3rd=cell*math.sqrt(3)                                      
   return cell_3rd

def Dsr_4th(r, cell, random_distances):
   if random_distances:
      cell_4th=cell*(1.581 + 1.376*r)                                #''' ------------------4th shell ::: [012], [021], [102], [120], [201], [210] ::: 24 neighbours ------------------ '''
   else:
      cell_4th=cell*math.sqrt(5)
   return cell_4th

def Dsr_5th(r, cell, random_distances):
   if random_distances:
      cell_5th=cell*(1.5 + 1.098*r)                                    #''' ------------------5th shell ::: [002], [020], [200] ::: 6 neighbours ------------------ '''
   else:
      cell_5th=cell*2
   return cell_5th

def Dsr_6th(r, cell, random_distances):
   if random_distances:
      cell_6th=cell*(2.12132 + 1.44939*r)                          #''' ------------------6th shell ::: [022], [202], [220] ::: 12 neighbours ------------------ '''
   else:
      cell_6th=cell*math.sqrt(8)
   return cell_6th

def Dsr_7th(r, cell, random_distances):
   if random_distances:
      cell_7th= None         # !!! to be constructed                  #''' ------------------7th shell ::: [112], [121], [211] ::: 24 neighbours ------------------ '''
   else:
      cell_7th=cell*math.sqrt(6)                                      
   return cell_7th

def Dsr_8th(r, cell, random_distances):
   if random_distances:
      cell_8th= None         # !!! to be constructed                  #''' ------------------8th shell ::: [122], [212], [221] ::: 24 neighbours ------------------ '''
   else:
      cell_8th=cell*3                                 
   return cell_8th

def Dsr_9th(r, cell, random_distances):
   if random_distances:
      cell_9th= None         # !!! to be constructed                  #''' ------------------9th shell ::: [222] ::: 8 neighbours ------------------------------ '''
   else:
      cell_9th=cell*math.sqrt(12)                                 
   return cell_9th


''' CONSTRAINTS of the Groups '''
constrains_of_group_9 =   False
constrains_of_group_14 = False
constrains_of_group_17 = False

        
def GG_process(attr1, attr3, attr4, attr5, attr6, attr7, attr8, attr9, attr10, attr11, attr12, attr13, attr14, attr15, ATTR2):
       # q= Queue(), attr1= GDSM, attr2= Selection, attr3= faza, attr4= asc, attr5= vg, attr6= S, attr7= taula, attr8= T_next, attr9= T, attr10= likvid, attr11= R, attr12= cell, attr13= W, attr14= cas, attr15= random_distances 
       if attr1 == 'new':
           seznam_premikov= []           # NEW

       for attr2 in ATTR2:
           grain=attr2[0]
           s=attr2[1]
           ''' ----------------------------------------- GROUP 1 ::: [001], [010] >>> 4 sites -------------------------------------------'''
           if s == '001':
               dij = np.array([0,0,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza001 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_1st(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain) , grain, attr3)
               seznam_premikov.append(faza001); del faza001

           elif s == '00_1':
               dij = np.array([0,0,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza00_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)& (lij>=Dsr_1st(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain), grain, attr3)
               seznam_premikov.append(faza00_1); del faza00_1

           elif s == '010':
               dij = np.array([0,1,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza010 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_1st(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain), grain, attr3)
               seznam_premikov.append(faza010); del faza010

           elif s == '0_10':
               dij = np.array([0,-1,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza0_10 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_1st(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain), grain, attr3)
               seznam_premikov.append(faza0_10); del faza0_10
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
               ''' ----------------------------------------- GROUP 2 ::: [100] >>> 1 site -------------------------------------------'''
           elif s == '100':
               dij = np.array([1,0,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza100 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)& (lij>=Dsr_1st(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain), grain, attr3)
               seznam_premikov.append(faza100); del faza100
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                        
               ''' ----------------------------------------- GROUP 3 ::: [_100] >>> 1 site -------------------------------------------'''
           elif s == '_100':
               dij = np.array([-1,0,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_100 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_1st(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain), grain, attr3)
               seznam_premikov.append(faza_100); del faza_100
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 4 ::: [011] >>> 4 sites -------------------------------------------'''
           elif s == '011':
               dij = np.array([0,1,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza011 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_2nd(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)& (
                                                                      #(np.pad(attr3,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)         |              # 001
                                                                      #(np.pad(attr3,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)         |              # 010
                                                                      #(np.pad(attr3,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)   |              # 111
                                                                      #(np.pad(attr3,((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==grain)                   # _111
                                                True), grain, attr3)
               seznam_premikov.append(faza011); del faza011

           elif s == '01_1':
               dij = np.array([0,1,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza01_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_2nd(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)& (
                                                                         #(np.pad(attr3,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)            |          # 00_1
                                                                         #(np.pad(attr3,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)           |          # 010
                                                                         #(np.pad(attr3,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)      |          # 11_1
                                                                         #(np.pad(attr3,((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==grain)                  # _11_1
                                                True), grain, attr3)
               seznam_premikov.append(faza01_1); del faza01_1

           elif s == '0_11':
               dij = np.array([0,-1,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza0_11 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_2nd(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)& (
                                                                         #(np.pad(attr3,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)           |          # 001
                                                                         #(np.pad(attr3,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)           |          # 0_10
                                                                         #(np.pad(attr3,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)     |           # 1_11
                                                                         #(np.pad(attr3,((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==grain)                  # _1_11
                                                True), grain, attr3)
               seznam_premikov.append(faza0_11); del faza0_11

           elif s == '0_1_1':
               dij = np.array([0,-1,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza0_1_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_2nd(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)& (
                                                                            #(np.pad(attr3,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)         |         # 00_1
                                                                            #(np.pad(attr3,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)       |         # 0_10
                                                                            #(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)    |         # 1_1_1
                                                                            #(np.pad(attr3,((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==grain)               # _1_1_1
                                                True), grain, attr3)
               seznam_premikov.append(faza0_1_1); del faza0_1_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 5 ::: [101], [110] >>> 4 sites -------------------------------------------'''
           elif s == '101':
               dij = np.array([1,0,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza101 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_2nd(attr11, attr12, attr15))&(np.pad(attr3,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==grain)& (
                                                                       #(np.pad(attr3,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain)             |        # 100
                                                                       #(np.pad(attr3,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)             |        # 001          
                                                                       #(np.pad(attr3,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)       |        # 111      
                                                                       #(np.pad(attr3,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)                 # 1_11         
                                                True), grain, attr3)
               seznam_premikov.append(faza101); del faza101

           elif s == '10_1':
               dij = np.array([1,0,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza10_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_2nd(attr11, attr12, attr15))&(np.pad(attr3,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==grain)& (
                                                                         #(np.pad(attr3,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain)            |       # 100
                                                                         #(np.pad(attr3,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)             |       # 00_1          
                                                                         #(np.pad(attr3,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)       |       # 11_1      
                                                                         #(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)                # 1_1_1
                                               True), grain, attr3)
               seznam_premikov.append(faza10_1); del faza10_1

           elif s == '110':
               dij = np.array([1,1,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza110 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_2nd(attr11, attr12, attr15))&(np.pad(attr3,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==grain)& (
                                                                      #(np.pad(attr3,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain)             |         # 100
                                                                      #(np.pad(attr3,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)             |         # 010          
                                                                      #(np.pad(attr3,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)       |           # 111      
                                                                      #(np.pad(attr3,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)                    # 11_1
                                               True), grain, attr3)
               seznam_premikov.append(faza110); del faza110

           elif s == '1_10':
               dij = np.array([1,-1,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza1_10 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_2nd(attr11, attr12, attr15))&(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)& (
                                                                         #(np.pad(attr3,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==grain)          |      # 100
                                                                         #(np.pad(attr3,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)           |      # 0_10          
                                                                         #(np.pad(attr3,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)     |       # 1_11      
                                                                         #(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)              # 1_1_1
                                                True), grain, attr3)
               seznam_premikov.append(faza1_10); del faza1_10
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 6 ::: [_101], [_110] >>> 4 sites -------------------------------------------'''
           elif s == '_101':
               dij = np.array([-1,0,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_101 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_2nd(attr11, attr12, attr15))&(np.pad(attr3,((0,1),(0,0), (1,0)), 'constant')[1:,:,:-1]==grain)& (
                                                                         #(np.pad(attr3,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain)          |       # _100
                                                                         #(np.pad(attr3,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)         |       # 001          
                                                                         #(np.pad(attr3,((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==grain)    |        # _111      
                                                                         #(np.pad(attr3,((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==grain)             # _1_11
                                                True), grain, attr3)
               seznam_premikov.append(faza_101); del faza_101

           elif s == '_10_1':
               dij = np.array([-1,0,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_10_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_2nd(attr11, attr12, attr15))&(np.pad(attr3,((0,1),(0,0), (0,1)), 'constant')[1:,:,1:]==grain)& (
                                                                            #(np.pad(attr3,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain)           |       # _100
                                                                            #(np.pad(attr3,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)           |       # 00_1          
                                                                            #(np.pad(attr3,((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==grain)       |        # _11_1      
                                                                            #(np.pad(attr3,((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==grain)                # _1_1_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_10_1); del faza_10_1

           elif s == '_110':
               dij = np.array([-1,1,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_110 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_2nd(attr11, attr12, attr15))&(np.pad(attr3,((0,1),(1,0), (0,0)), 'constant')[1:,:-1,:]==grain)& (
                                                                         #(np.pad(attr3,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain)         |         # _100
                                                                         #(np.pad(attr3,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)        |         # 010          
                                                                         #(np.pad(attr3,((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==grain)    |        # _111      
                                                                         #(np.pad(attr3,((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==grain)              # _11_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_110); del faza_110

           elif s == '_1_10':
               dij = np.array([-1,-1,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_1_10 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_2nd(attr11, attr12, attr15))&(np.pad(attr3,((0,1),(0,1), (0,0)), 'constant')[1:,1:,:]==grain)& (
                                                                            #(np.pad(attr3,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==grain)       |      # _100
                                                                            #(np.pad(attr3,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)         |      # 0_10             
                                                                            #(np.pad(attr3,((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==grain)     |      # _1_11      
                                                                            #(np.pad(attr3,((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==grain)             # _1_1_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_1_10); del faza_1_10
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 7 ::: [111] >>> 4 sites -------------------------------------------'''
           elif s == '111':
               dij = np.array([1,1,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza111 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_3rd(attr11, attr12, attr15))&(np.pad(attr3,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)& (
                                                                      #(np.pad(attr3,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==grain)    |        # 101
                                                                      #(np.pad(attr3,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==grain)    |        # 110          
                                                                      #(np.pad(attr3,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)             # 011      
                                                True), grain, attr3)
               seznam_premikov.append(faza111); del faza111

           elif s == '11_1':
               dij = np.array([1,1,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza11_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_3rd(attr11, attr12, attr15))&(np.pad(attr3,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)& (
                                                                         #(np.pad(attr3,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==grain)     |       # 10_1
                                                                         #(np.pad(attr3,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==grain)    |       # 110          
                                                                         #(np.pad(attr3,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)             # 01_1
                                                True), grain, attr3)
               seznam_premikov.append(faza11_1); del faza11_1

           elif s == '1_11':
               dij = np.array([1,-1,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza1_11 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_3rd(attr11, attr12, attr15))&(np.pad(attr3,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)& (
                                                                         #(np.pad(attr3,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==grain)    |       # 101
                                                                         #(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)     |       # 1_10          
                                                                         #(np.pad(attr3,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)             # 0_11
                                                True), grain, attr3)
               seznam_premikov.append(faza1_11); del faza1_11

           elif s == '1_1_1':
               dij = np.array([1,-1,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza1_1_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_3rd(attr11, attr12, attr15))&(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)& (
                                                                            #(np.pad(attr3,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==grain)     |        # 10_1
                                                                            #(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)     |        # 1_10          
                                                                            #(np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)               # 0_1_1
                                                True), grain, attr3)
               seznam_premikov.append(faza1_1_1); del faza1_1_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 8 ::: [_111] >>> 4 sites -------------------------------------------'''
           elif s == '_111':
               dij = np.array([-1,1,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_111 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_3rd(attr11, attr12, attr15))&(np.pad(attr3,((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==grain)& (
                                                                         #(np.pad(attr3,((0,1),(0,0), (1,0)), 'constant')[1:,:,:-1]==grain)     |       # _101
                                                                         #(np.pad(attr3,((0,1),(1,0), (0,0)), 'constant')[1:,:-1,:]==grain)     |       # _110          
                                                                         #(np.pad(attr3,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)            # 011
                                                True), grain, attr3)
               seznam_premikov.append(faza_111); del faza_111

           elif s == '_11_1':
               dij = np.array([-1,1,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_11_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_3rd(attr11, attr12, attr15))&(np.pad(attr3,((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==grain)& (
                                                                            #(np.pad(attr3,((0,1),(0,0), (0,1)), 'constant')[1:,:,1:]==grain)        |      # _10_1
                                                                            #(np.pad(attr3,((0,1),(1,0), (0,0)), 'constant')[1:,:-1,:]==grain)       |      # _110          
                                                                            #(np.pad(attr3,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)              # 01_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_11_1); del faza_11_1

           elif s == '_1_11':
               dij = np.array([-1,-1,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_1_11 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_3rd(attr11, attr12, attr15))&(np.pad(attr3,((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==grain)& (
                                                                            #(np.pad(attr3,((0,1),(0,0), (1,0)), 'constant')[1:,:,:-1]==grain)       |      # _101
                                                                            #(np.pad(attr3,((0,1),(0,1), (0,0)), 'constant')[1:,1:,:]==grain)        |      # _1_10          
                                                                            #(np.pad(attr3,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)              # 0_11                
                                                True), grain, attr3)
               seznam_premikov.append(faza_1_11); del faza_1_11

           elif s == '_1_1_1':
               dij = np.array([-1,-1,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_1_1_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0)&(lij>=Dsr_3rd(attr11, attr12, attr15))&(np.pad(attr3,((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==grain)& (
                                                                               #(np.pad(attr3,((0,1),(0,0), (0,1)), 'constant')[1:,:,1:]==grain)    |         # _10_1
                                                                               #(np.pad(attr3,((0,1),(0,1), (0,0)), 'constant')[1:,1:,:]==grain)    |        # _1_10          
                                                                               #(np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)             # 0_1_1
                                                 True), grain, attr3)
               seznam_premikov.append(faza_1_1_1); del faza_1_1_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 9 ::: [012], [021] >>> 8 sites -------------------------------------------'''
           elif s == '012':
               dij = np.array([0,1,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_9:
                   CG9=(#(np.pad(attr3,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)       |      # 002
                              #(np.pad(attr3,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)    |      # 021
                              #(np.pad(attr3,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==grain)    |      # 022
                              (np.pad(attr3,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)       |      # 001
                              (np.pad(attr3,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)     )   # 011
               else:
                   CG9 = True
               faza012 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)& CG9, grain, attr3) 
               seznam_premikov.append(faza012); del faza012

           elif s == '01_2':
               dij = np.array([0,1,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(attr3,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)         |         # 00_2
                                #(np.pad(attr3,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)       |         # 02_1
                                #(np.pad(attr3,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==grain)       |         # 02_2
                                 (np.pad(attr3,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)          |        # 00_1
                                 (np.pad(attr3,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)       )      # 01_1      
               else:
                   CG9 = True
               faza01_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)& CG9, grain, attr3) 
               seznam_premikov.append(faza01_2); del faza01_2

           elif s == '0_12':
               dij = np.array([0,-1,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(attr3,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)      |         # 002
                                #(np.pad(attr3,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)     |         # 0_21
                                #(np.pad(attr3,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==grain)     |       # 0_22
                                (np.pad(attr3,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)       |         # 001
                                (np.pad(attr3,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)     )       # 0_11         
               else:
                   CG9 = True
               faza0_12 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)& CG9, grain, attr3)
               seznam_premikov.append(faza0_12); del faza0_12

           elif s == '0_1_2':
               dij = np.array([0,-1,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(attr3,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)       |        # 00_2
                                #(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)      |        # 0_2_1
                                #(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)      |        # 0_2_2
                                (np.pad(attr3,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)        |        # 00_1
                                (np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)       )     # 0_1_1         
               else:
                   CG9 = True
               faza0_1_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)& CG9, grain, attr3) 
               seznam_premikov.append(faza0_1_2); del faza0_1_2

           elif s == '021':
               dij = np.array([0,2,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(attr3,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)        |         # 020
                                #(np.pad(attr3,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==grain)      |         # 022
                                #(np.pad(attr3,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)      |        # 012
                                (np.pad(attr3,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)         |         # 010
                                (np.pad(attr3,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)      )       # 011         
               else:
                   CG9 = True
               faza021 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)& CG9, grain, attr3)
               seznam_premikov.append(faza021); del faza021

           elif s == '02_1':
               dij = np.array([0,2,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(attr3,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)       |           # 020
                                #(np.pad(attr3,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==grain)      |           # 02_2
                                #(np.pad(attr3,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)      |           # 01_2
                                (np.pad(attr3,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)         |          # 010
                                (np.pad(attr3,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)        )       # 01_1
               else:
                   CG9 = True
               faza02_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)& CG9, grain, attr3) 
               seznam_premikov.append(faza02_1); del faza02_1

           elif s == '0_21':
               dij = np.array([0,-2,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)        |          # 0_20
                                #(np.pad(attr3,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==grain)      |         # 0_22
                                #(np.pad(attr3,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)    |          # 0_12
                                (np.pad(attr3,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)         |         # 0_10
                                (np.pad(attr3,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)      )       # 0_11       
               else:
                   CG9 = True
               faza0_21 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)& CG9, grain, attr3)
               seznam_premikov.append(faza0_21); del faza0_21

           elif s == '0_2_1':
               dij = np.array([0,-2,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_9:
                   CG9 = (#(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)     |       # 0_20
                                #(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)    |      # 0_2_2
                                #(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |       # 0_1_2
                                (np.pad(attr3,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)      |       # 0_10
                                (np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    )     # 0_1_1             
               else:
                   CG9 = True
               faza0_2_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)& CG9, grain, attr3) 
               seznam_premikov.append(faza0_2_1); del faza0_2_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        

               ''' ----------------------------------------- GROUP 10 ::: [102], [120] >>> 4 sites -------------------------------------------'''
           elif s == '102':
               dij = np.array([1,0,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza102 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==grain)& (
                                                                       #(np.pad(attr3,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==grain)       |       # 101
                                                                       #(np.pad(attr3,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)          |       # 002
                                                                       #(np.pad(attr3,((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==grain)    |       # 112
                                                                       #(np.pad(attr3,((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==grain)     |       # 1_12
                                                                       #(np.pad(attr3,((2,0),(0,0), (2,0)), 'constant')[:-2,:,:-2]==grain)               # 202
                                                True), grain, attr3) 
               seznam_premikov.append(faza102); del faza102

           elif s == '10_2':
               dij = np.array([1,0,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza10_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==grain)& (
                                                                       #(np.pad(attr3,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==grain)       |       # 10_1
                                                                       #(np.pad(attr3,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)          |       # 00_2
                                                                       #(np.pad(attr3,((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==grain)    |       # 11_2
                                                                       #(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)     |       # 1_1_2
                                                                       #(np.pad(attr3,((2,0),(0,0), (0,2)), 'constant')[:-2,:,2:]==grain)               # 20_2
                                                True), grain, attr3) 
               seznam_premikov.append(faza10_2); del faza10_2

           elif s == '120':
               dij = np.array([1,2,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza120 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==grain)& (
                                                                       #(np.pad(attr3,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==grain)       |       # 110
                                                                       #(np.pad(attr3,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)          |       # 020
                                                                       #(np.pad(attr3,((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==grain)    |       # 121
                                                                       #(np.pad(attr3,((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==grain)     |       # 12_1
                                                                       #(np.pad(attr3,((2,0),(2,0), (0,0)), 'constant')[:-2,:-2,:]==grain)               # 220
                                                True), grain, attr3) 
               seznam_premikov.append(faza120); del faza120

           elif s == '1_20':
               dij = np.array([1,-2,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza1_20 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)& (
                                                                       #(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       #(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       #(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       #(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       #(np.pad(attr3,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, attr3) 
               seznam_premikov.append(faza1_20); del faza1_20
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
        
               ''' ----------------------------------------- GROUP 11 ::: [_102], [_120] >>> 4 sites -------------------------------------------'''
           elif s == '_102':
               dij = np.array([-1,0,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_102 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(0,0), (2,0)), 'constant')[1:,:,:-2]==grain)& (
                                                                       #wrong(np.pad(attr3,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==grain)       |       # 101
                                                                       #wrong(np.pad(attr3,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)          |       # 002
                                                                       #wrong(np.pad(attr3,((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==grain)    |       # 112
                                                                       #wrong(np.pad(attr3,((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==grain)     |       # 1_12
                                                                       #wrong(np.pad(attr3,((2,0),(0,0), (2,0)), 'constant')[:-2,:,:-2]==grain)               # 202
                                                True), grain, attr3) 
               seznam_premikov.append(faza_102); del faza_102

           elif s == '_10_2':
               dij = np.array([-1,0,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_10_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(0,0), (0,2)), 'constant')[1:,:,2:]==grain)& (
                                                                       #wrong(np.pad(attr3,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==grain)       |       # 10_1
                                                                       #wrong(np.pad(attr3,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)          |       # 00_2
                                                                       #wrong(np.pad(attr3,((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==grain)    |       # 11_2
                                                                       #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)     |       # 1_1_2
                                                                       #wrong(np.pad(attr3,((2,0),(0,0), (0,2)), 'constant')[:-2,:,2:]==grain)               # 20_2
                                               True), grain, attr3) 
               seznam_premikov.append(faza_10_2); del faza_10_2

           elif s == '_120':
               dij = np.array([-1,2,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_120 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(2,0), (0,0)), 'constant')[1:,:-2,:]==grain)& (
                                                                       #(np.pad(attr3,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==grain)       |       # 110
                                                                       #(np.pad(attr3,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)          |       # 020
                                                                       #(np.pad(attr3,((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==grain)    |       # 121
                                                                       #(np.pad(attr3,((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==grain)     |       # 12_1
                                                                       #(np.pad(attr3,((2,0),(2,0), (0,0)), 'constant')[:-2,:-2,:]==grain)               # 220
                                                True), grain, attr3) 
               seznam_premikov.append(faza_120); del faza_120

           elif s == '_1_20':
               dij = np.array([-1,-2,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_1_20 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(0,2), (0,0)), 'constant')[1:,2:,:]==grain)& (
                                                                       #(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       #(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       #(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       #(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       #(np.pad(attr3,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, attr3) 
               seznam_premikov.append(faza_1_20); del faza_1_20
               

               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 12 ::: [201], [210] >>> 4 sites -------------------------------------------'''
           elif s == '201':
               dij = np.array([2,0,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza201 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,0), (1,0)), 'constant')[:-2,:,:-1]==grain)& (
                                                                       # wrong(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(attr3,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, attr3) 
               seznam_premikov.append(faza201); del faza201

           elif s == '20_1':
               dij = np.array([2,0,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza20_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,0), (0,1)), 'constant')[:-2,:,1:]==grain)& (
                                                                       # wrong(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(attr3,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, attr3) 
               seznam_premikov.append(faza20_1); del faza20_1

           elif s == '210':
               dij = np.array([2,1,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza210 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(1,0), (0,0)), 'constant')[:-2,:-1,:]==grain)& (
                                                                       # wrong(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(attr3,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, attr3) 
               seznam_premikov.append(faza210); del faza210

           elif s == '2_10':
               dij = np.array([2,-1,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza2_10 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,1), (0,0)), 'constant')[:-2,1:,:]==grain)& (
                                                                       # wrong(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(attr3,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, attr3) 
               seznam_premikov.append(faza2_10); del faza2_10

               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 13 ::: [_201], [_210] >>> 4 sites -------------------------------------------'''
           elif s == '_201':
               dij = np.array([-2,0,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_201 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,0), (1,0)), 'constant')[2:,:,:-1]==grain)& (
                                                                       # wrong(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(attr3,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, attr3) 
               seznam_premikov.append(faza_201); del faza_201

           elif s == '_20_1':
               dij = np.array([-2,0,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_20_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,0), (0,1)), 'constant')[2:,:,1:]==grain)& (
                                                                       # wrong(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(attr3,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, attr3) 
               seznam_premikov.append(faza_20_1); del faza_20_1

           elif s == '_210':
               dij = np.array([-2,1,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_210 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(1,0), (0,0)), 'constant')[2:,:-1,:]==grain)& (
                                                                       # wrong(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(attr3,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, attr3) 
               seznam_premikov.append(faza_210); del faza_210

           elif s == '_2_10':
               dij = np.array([-2,-1,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_2_10 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_4th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,1), (0,0)), 'constant')[2:,1:,:]==grain)& (
                                                                       # wrong(np.pad(attr3,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==grain)       |       # 1_10
                                                                       # wrong(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)          |       # 0_20
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)    |       # 1_21
                                                                       # wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)     |       # 1_2_1
                                                                       # wrong(np.pad(attr3,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)               # 2_20
                                                True), grain, attr3) 
               seznam_premikov.append(faza_2_10); del faza_2_10
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 14 ::: [002], [020] >>> 4 sites -------------------------------------------'''
           elif s == '002':
               dij = np.array([0,0,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_14:
                   CG14 = (#(np.pad(attr3,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)      |        # 012           
                                  #(np.pad(attr3,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)       |         # 0_12
                                  (np.pad(attr3,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain)         |         # 001
                                  (np.pad(attr3,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)       |        # 011
                                  (np.pad(attr3,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)        )      # 0_11      
               else:
                   CG14 = True
               faza002 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_5th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==grain)& CG14, grain, attr3) 
               seznam_premikov.append(faza002); del faza002

           elif s == '00_2':
               dij = np.array([0,0,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_14:
                   CG14 = (#(np.pad(attr3,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)    |       # 01_2
                                  #(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)     |       # 0_1_2
                                  (np.pad(attr3,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==grain)       |       # 00_1
                                  (np.pad(attr3,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)     |      # 01_1
                                  (np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)      )    # 0_1_1    
               else:
                   CG14 = True
               faza00_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_5th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==grain)& CG14, grain, attr3)
               seznam_premikov.append(faza00_2); del faza00_2

           elif s == '020':
               dij = np.array([0,2,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_14:
                   CG14 = (#(np.pad(attr3,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)     |         # 021
                                  #(np.pad(attr3,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)      |       # 02_1
                                  (np.pad(attr3,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==grain)        |         # 010 
                                  (np.pad(attr3,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)     |         # 011
                                  (np.pad(attr3,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)      )       # 01_1
               else:
                   CG14 = True
               faza020 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_5th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==grain)& CG14, grain, attr3)
               seznam_premikov.append(faza020); del faza020

           elif s == '0_20':
               dij = np.array([0,-2,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_14:
                   CG14 = (#(np.pad(attr3,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)    |        # 0_21
                                  #(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)     |        # 0_2_1
                                  (np.pad(attr3,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==grain)        |        # 0_10
                                  (np.pad(attr3,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)     |        # 0_11 
                                  (np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)      )      # 0_1_1
               else:
                   CG14 = True
               faza0_20 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_5th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==grain)& CG14, grain, attr3)     
               seznam_premikov.append(faza0_20); del faza0_20
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 15 ::: [200] >>> 1 site -------------------------------------------'''
           elif s == '200':
               dij = np.array([2,0,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza200 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_5th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,0), (0,0)), 'constant')[:-2,:,:]==grain)& (
                                                True), grain, attr3)     
               seznam_premikov.append(faza200); del faza200
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 16 ::: [_200] >>> 1 site -------------------------------------------'''
           elif s == '_200':
               dij = np.array([-2,0,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_200 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_5th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,0), (0,0)), 'constant')[2:,:,:]==grain)& (
                                                True), grain, attr3)     
               seznam_premikov.append(faza_200); del faza_200
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 17 ::: [022] >>> 4 sites -------------------------------------------'''
           elif s == '022':
               dij = np.array([0,2,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_17:
                   CG17 = (#(np.pad(attr3,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)      |       # 012
                                  #(np.pad(attr3,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)      |       # 021
                                  #(np.pad(attr3,((0,0),(1,0), (3,0)), 'constant')[:,:-1,:-3]==grain)      |       # 013
                                  #(np.pad(attr3,((0,0),(3,0), (1,0)), 'constant')[:,:-3,:-1]==grain)      |       # 031
                                  (np.pad(attr3,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==grain)       )    # 011
               else:
                   CG17 = True
               faza022 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_6th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==grain)& CG17, grain, attr3) 
               seznam_premikov.append(faza022); del faza022

           elif s == '02_2':
               dij = np.array([0,2,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_17:
                   CG17 = (#(np.pad(attr3,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)    |        # 01_2
                                  #(np.pad(attr3,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)    |        # 02_1
                                  #(np.pad(attr3,((0,0),(1,0), (0,3)), 'constant')[:,:-1,3:]==grain)    |        # 01_3
                                  #(np.pad(attr3,((0,0),(3,0), (0,1)), 'constant')[:,:-3,1:]==grain)   |        # 03_1
                                  (np.pad(attr3,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==grain)     )      # 01_1
               else:
                   CG17 = True
               faza02_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_6th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==grain)& CG17, grain, attr3) 
               seznam_premikov.append(faza02_2); del faza02_2

           elif s == '0_22':
               dij = np.array([0,-2,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_17:
                   CG17 = (#(np.pad(attr3,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)     |        # 0_12
                                  #(np.pad(attr3,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)     |        # 0_21
                                  #(np.pad(attr3,((0,0),(0,1), (3,0)), 'constant')[:,1:,:-3]==grain)     |        # 0_13
                                  #(np.pad(attr3,((0,0),(0,3), (1,0)), 'constant')[:,3:,:-1]==grain)     |        # 0_31
                                  (np.pad(attr3,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==grain)      )     # 0_11
               else:
                   CG17 = True
               faza0_22 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_6th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==grain)& CG17, grain, attr3) 
               seznam_premikov.append(faza0_22); del faza0_22

           elif s == '0_2_2':
               dij = np.array([0,-2,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               if constrains_of_group_17:
                   CG17 = (#(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)     |      # 0_1_2                
                                  #(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)     |      # 0_2_1
                                  #(np.pad(attr3,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)     |      # 0_1_3
                                  #(np.pad(attr3,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)     |     # 0_3_1
                                  (np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)      )    # 0_1_1
               else:
                   CG17 = True
               faza0_2_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_6th(attr11, attr12, attr15) )&(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)& CG17, grain, attr3)
               seznam_premikov.append(faza0_2_2); del faza0_2_2
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 18 ::: [202], [220] >>> 4 sites -------------------------------------------'''
           elif s == '202':
               dij = np.array([2,0,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza202 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_6th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,0), (2,0)), 'constant')[:-2,:,:-2]==grain)& (                                                 
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(attr3,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, attr3)
               seznam_premikov.append(faza202); del faza202

           elif s == '20_2':
               dij = np.array([2,0,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza20_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_6th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,0), (0,2)), 'constant')[:-2,:,2:]==grain)& (                                                 
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(attr3,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, attr3)
               seznam_premikov.append(faza20_2); del faza20_2

           elif s == '220':
               dij = np.array([2,2,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza220 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_6th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(2,0), (0,0)), 'constant')[:-2,:-2,:]==grain)& (                                                 
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(attr3,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, attr3)
               seznam_premikov.append(faza220); del faza220

           elif s == '2_20':
               dij = np.array([2,-2,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza2_20 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_6th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==grain)& (                                                 
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(attr3,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, attr3)
               seznam_premikov.append(faza2_20); del faza2_20
           

               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 19 ::: [_202], [_220] >>> 4 sites -------------------------------------------'''
           elif s == '_202':
               dij = np.array([-2,0,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_202 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_6th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,0), (2,0)), 'constant')[2:,:,:-2]==grain)& (                                                 
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(attr3,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_202); del faza_202

           elif s == '_20_2':
               dij = np.array([-2,0,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_20_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_6th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,0), (0,2)), 'constant')[2:,:,2:]==grain)& (                                                 
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(attr3,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_20_2); del faza_20_2

           elif s == '_220':
               dij = np.array([-2,2,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_220 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_6th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(2,0), (0,0)), 'constant')[2:,:-2,:]==grain)& (                                                 
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(attr3,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_220); del faza_220

           elif s == '_2_20':
               dij = np.array([-2,-2,0])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_2_20 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_6th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,2), (0,0)), 'constant')[2:,2:,:]==grain)& (                                                 
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)    |      # 0_1_2                
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==grain)    |      # 0_1_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         # 0_2_1
                                                                            #wrong(np.pad(attr3,((0,0),(0,1), (0,3)), 'constant')[:,1:,3:]==grain)    |      # 0_1_3
                                                                            #wrong(np.pad(attr3,((0,0),(0,3), (0,1)), 'constant')[:,3:,1:]==grain)         # 0_3_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_2_20); del faza_2_20
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 20 ::: [112], [121] >>> 8 sites -------------------------------------------'''
           elif s == '112':
               dij = np.array([1,1,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza112 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==grain)& (                                                 
                                                                           #(np.pad(attr3,((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==grain)       |      # 102                
                                                                           #(np.pad(attr3,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)    |      # 111
                                                                           #(np.pad(attr3,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)       |      # 012
                                                                           #(np.pad(attr3,((2,0),(1,0), (2,0)), 'constant')[:-2,:-1,:-2]==grain)    |      # 212
                                                                           #(np.pad(attr3,((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==grain)           # 122
                                               True ), grain, attr3)
               seznam_premikov.append(faza112); del faza112

           elif s == '11_2':
               dij = np.array([1,1,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza11_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==grain)& (                                                 
                                                                           #(np.pad(attr3,((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==grain)       |      # 10_2                
                                                                           #(np.pad(attr3,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)    |      # 11_1
                                                                           #(np.pad(attr3,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)       |      # 01_2
                                                                           #(np.pad(attr3,((2,0),(1,0), (0,2)), 'constant')[:-2,:-1,2:]==grain)    |      # 21_2
                                                                           #(np.pad(attr3,((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==grain)           # 12_2
                                                True), grain, attr3)
               seznam_premikov.append(faza11_2); del faza11_2

           elif s == '1_12':
               dij = np.array([1,-1,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza1_12 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==grain)& (                                                 
                                                                           #(np.pad(attr3,((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==grain)      |      # 102                
                                                                           #(np.pad(attr3,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)    |      # 1_11
                                                                           #(np.pad(attr3,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)       |      # 0_12
                                                                           #(np.pad(attr3,((2,0),(0,1), (2,0)), 'constant')[:-2,1:,:-2]==grain)    |      # 2_12
                                                                           #(np.pad(attr3,((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==grain)           # 1_22
                                                True), grain, attr3)
               seznam_premikov.append(faza1_12); del faza1_12

           elif s == '1_1_2':
               dij = np.array([1,-1,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza1_1_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)& (                                                 
                                                                           #(np.pad(attr3,((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==grain)      |      # 10_2                
                                                                           #(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)    |      # 1_1_1
                                                                           #(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)       |      # 0_1_2
                                                                           #(np.pad(attr3,((2,0),(0,1), (0,2)), 'constant')[:-2,1:,2:]==grain)    |      # 2_1_2
                                                                           #(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)           # 1_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza1_1_2); del faza1_1_2

           elif s == '121':
               dij = np.array([1,2,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza121 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==grain)& (                                                 
                                                                           #(np.pad(attr3,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)      |      # 111                
                                                                           #(np.pad(attr3,((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==grain)         |      # 120
                                                                           #(np.pad(attr3,((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==grain)      |      # 122
                                                                           #(np.pad(attr3,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)         |      # 021
                                                                           #(np.pad(attr3,((2,0),(2,0), (1,0)), 'constant')[:-2,:-2,:-1]==grain)             # 221
                                                True), grain, attr3)
               seznam_premikov.append(faza121); del faza121

           elif s == '12_1':
               dij = np.array([1,2,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza12_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==grain)& (                                                 
                                                                           #(np.pad(attr3,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)      |      # 11_1                
                                                                           #(np.pad(attr3,((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==grain)        |      # 120
                                                                           #(np.pad(attr3,((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==grain)      |      # 12_2
                                                                           #(np.pad(attr3,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)         |      # 02_1
                                                                           #(np.pad(attr3,((2,0),(2,0), (0,1)), 'constant')[:-2,:-2,1:]==grain)             # 22_1
                                                True), grain, attr3)
               seznam_premikov.append(faza12_1); del faza12_1

           elif s == '1_21':
               dij = np.array([1,-2,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza1_21 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)& (                                                 
                                                                           #(np.pad(attr3,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)      |      # 1_11                
                                                                           #(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)         |      # 1_20
                                                                           #(np.pad(attr3,((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==grain)      |      # 1_22
                                                                           #(np.pad(attr3,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)         |      # 0_21
                                                                           #(np.pad(attr3,((2,0),(0,2), (1,0)), 'constant')[:-2,2:,:-1]==grain)             # 2_21
                                                True), grain, attr3)
               seznam_premikov.append(faza1_21); del faza1_21

           elif s == '1_2_1':
               dij = np.array([1,-2,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza1_2_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)& (                                                 
                                                                           #(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #(np.pad(attr3,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, attr3)
               seznam_premikov.append(faza1_2_1); del faza1_2_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 21 ::: [_112], [_121] >>> 8 sites -------------------------------------------'''
           elif s == '_112':
               dij = np.array([-1,1,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_112 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(1,0), (2,0)), 'constant')[1:,:-1,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==grain)       |      # 102                
                                                                           #wrong(np.pad(attr3,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)    |      # 111
                                                                           #wrong(np.pad(attr3,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==grain)       |      # 012
                                                                           #wrong(np.pad(attr3,((2,0),(1,0), (2,0)), 'constant')[:-2,:-1,:-2]==grain)    |      # 212
                                                                           #wrong(np.pad(attr3,((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==grain)           # 122
                                               True ), grain, attr3)
               seznam_premikov.append(faza_112); del faza_112

           elif s == '_11_2':
               dij = np.array([-1,1,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_11_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(1,0), (0,2)), 'constant')[1:,:-1,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==grain)       |      # 10_2                
                                                                           #wrong(np.pad(attr3,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)    |      # 11_1
                                                                           #wrong(np.pad(attr3,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==grain)       |      # 01_2
                                                                           #wrong(np.pad(attr3,((2,0),(1,0), (0,2)), 'constant')[:-2,:-1,2:]==grain)    |      # 21_2
                                                                           #wrong(np.pad(attr3,((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==grain)           # 12_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_11_2); del faza_11_2

           elif s == '_1_12':
               dij = np.array([-1,-1,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_1_12 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(0,1), (2,0)), 'constant')[1:,1:,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==grain)      |      # 102                
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)    |      # 1_11
                                                                           #wrong(np.pad(attr3,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==grain)       |      # 0_12
                                                                           #wrong(np.pad(attr3,((2,0),(0,1), (2,0)), 'constant')[:-2,1:,:-2]==grain)    |      # 2_12
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==grain)           # 1_22
                                                True), grain, attr3)
               seznam_premikov.append(faza_1_12); del faza_1_12

           elif s == '_1_1_2':
               dij = np.array([-1,-1,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_1_1_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(0,1), (0,2)), 'constant')[1:,1:,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==grain)      |      # 10_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)    |      # 1_1_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==grain)       |      # 0_1_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,1), (0,2)), 'constant')[:-2,1:,2:]==grain)    |      # 2_1_2
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)           # 1_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_1_1_2); del faza_1_1_2

           elif s == '_121':
               dij = np.array([-1,2,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_121 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(2,0), (1,0)), 'constant')[1:,:-2,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==grain)      |      # 111                
                                                                           #wrong(np.pad(attr3,((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==grain)         |      # 120
                                                                           #wrong(np.pad(attr3,((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==grain)      |      # 122
                                                                           #wrong(np.pad(attr3,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==grain)         |      # 021
                                                                           #wrong(np.pad(attr3,((2,0),(2,0), (1,0)), 'constant')[:-2,:-2,:-1]==grain)             # 221
                                                True), grain, attr3)
               seznam_premikov.append(faza_121); del faza_121

           elif s == '_12_1':
               dij = np.array([-1,2,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_12_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(2,0), (0,1)), 'constant')[1:,:-2,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==grain)      |      # 11_1                
                                                                           #wrong(np.pad(attr3,((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==grain)        |      # 120
                                                                           #wrong(np.pad(attr3,((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==grain)      |      # 12_2
                                                                           #wrong(np.pad(attr3,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==grain)         |      # 02_1
                                                                           #wrong(np.pad(attr3,((2,0),(2,0), (0,1)), 'constant')[:-2,:-2,1:]==grain)             # 22_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_12_1); del faza_12_1

           elif s == '_1_21':
               dij = np.array([-1,-2,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_1_21 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(0,2), (1,0)), 'constant')[1:,2:,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==grain)      |      # 1_11                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)         |      # 1_20
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==grain)      |      # 1_22
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==grain)         |      # 0_21
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (1,0)), 'constant')[:-2,2:,:-1]==grain)             # 2_21
                                                True), grain, attr3)
               seznam_premikov.append(faza_1_21); del faza_1_21

           elif s == '_1_2_1':
               dij = np.array([-1,-2,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_1_2_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(0,2), (0,1)), 'constant')[1:,2:,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_1_2_1); del faza_1_2_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 22 ::: [211] >>> 4 sites -------------------------------------------'''
           elif s == '211':
               dij = np.array([2,1,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza211 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(1,0), (1,0)), 'constant')[:-2,:-1,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, attr3)
               seznam_premikov.append(faza211); del faza211

           elif s == '21_1':
               dij = np.array([2,1,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza21_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(1,0), (0,1)), 'constant')[:-2,:-1,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, attr3)
               seznam_premikov.append(faza21_1); del faza21_1

           elif s == '2_11':
               dij = np.array([2,-1,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza2_11 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,1), (1,0)), 'constant')[:-2,1:,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, attr3)
               seznam_premikov.append(faza2_11); del faza2_11

           elif s == '2_1_1':
               dij = np.array([2,-1,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza2_1_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,1), (0,1)), 'constant')[:-2,1:,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, attr3)
               seznam_premikov.append(faza2_1_1); del faza2_1_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 23 ::: [_211] >>> 4 sites -------------------------------------------'''
           elif s == '_211':
               dij = np.array([-2,1,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_211 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(1,0), (1,0)), 'constant')[2:,:-1,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_211); del faza_211

           elif s == '_21_1':
               dij = np.array([-2,1,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_21_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(1,0), (0,1)), 'constant')[2:,:-1,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_21_1); del faza_21_1

           elif s == '_2_11':
               dij = np.array([-2,-1,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_2_11 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,1), (1,0)), 'constant')[2:,1:,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_2_11); del faza_2_11

           elif s == '_2_1_1':
               dij = np.array([-2,-1,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_2_1_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_7th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,1), (0,1)), 'constant')[2:,1:,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==grain)      |      # 1_1_1                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==grain)        |      # 1_20
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)      |      # 1_2_2
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==grain)         |      # 0_2_1
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)             # 2_2_1
                                                True), grain, attr3)
               seznam_premikov.append(faza_2_1_1); del faza_2_1_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 24 ::: [122] >>> 4 sites -------------------------------------------'''
           elif s == '122':
               dij = np.array([1,2,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza122 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==grain)& (                                                 
                                                                           #(np.pad(attr3,((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==grain)      |      # 112                
                                                                           #(np.pad(attr3,((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==grain)      |      # 121
                                                                           #(np.pad(attr3,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==grain)         |      # 022
                                                                           #(np.pad(attr3,((2,0),(2,0), (2,0)), 'constant')[:-2,:-2,:-2]==grain)             # 222
                                                True), grain, attr3)
               seznam_premikov.append(faza122); del faza122

           elif s == '12_2':
               dij = np.array([1,2,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza12_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==grain)& (                                                 
                                                                           #(np.pad(attr3,((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==grain)      |      # 11_2                
                                                                           #(np.pad(attr3,((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==grain)      |      # 12_1
                                                                           #(np.pad(attr3,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==grain)         |      # 02_2
                                                                           #(np.pad(attr3,((2,0),(2,0), (0,2)), 'constant')[:-2,:-2,2:]==grain)             # 22_2
                                                True), grain, attr3)
               seznam_premikov.append(faza12_2); del faza12_2

           elif s == '1_22':
               dij = np.array([1,-2,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza1_22 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==grain)& (                                                 
                                                                           #(np.pad(attr3,((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==grain)      |      # 1_12                
                                                                           #(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)      |      # 1_21
                                                                           #(np.pad(attr3,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==grain)         |      # 0_22
                                                                           #(np.pad(attr3,((2,0),(0,2), (2,0)), 'constant')[:-2,2:,:-2]==grain)             # 2_22
                                                True), grain, attr3)
               seznam_premikov.append(faza1_22); del faza1_22

           elif s == '1_2_2':
               dij = np.array([1,-2,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza1_2_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==grain)& (                                                 
                                                                           #(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza1_2_2); del faza1_2_2
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 25 ::: [_122] >>> 4 sites -------------------------------------------'''
           elif s == '_122':
               dij = np.array([-1,2,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_122 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(2,0), (2,0)), 'constant')[1:,:-2,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==grain)      |      # 112                
                                                                           #wrong(np.pad(attr3,((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==grain)      |      # 121
                                                                           #wrong(np.pad(attr3,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==grain)         |      # 022
                                                                           #wrong(np.pad(attr3,((2,0),(2,0), (2,0)), 'constant')[:-2,:-2,:-2]==grain)             # 222
                                                True), grain, attr3)
               seznam_premikov.append(faza_122); del faza_122

           elif s == '_12_2':
               dij = np.array([-1,2,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_12_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(2,0), (0,2)), 'constant')[1:,:-2,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==grain)      |      # 11_2                
                                                                           #wrong(np.pad(attr3,((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==grain)      |      # 12_1
                                                                           #wrong(np.pad(attr3,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==grain)         |      # 02_2
                                                                           #wrong(np.pad(attr3,((2,0),(2,0), (0,2)), 'constant')[:-2,:-2,2:]==grain)             # 22_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_12_2); del faza_12_2

           elif s == '_1_22':
               dij = np.array([-1,-2,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_1_22 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(0,2), (2,0)), 'constant')[1:,2:,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==grain)      |      # 1_12                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==grain)      |      # 1_21
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==grain)         |      # 0_22
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (2,0)), 'constant')[:-2,2:,:-2]==grain)             # 2_22
                                                True), grain, attr3)
               seznam_premikov.append(faza_1_22); del faza_1_22

           elif s == '_1_2_2':
               dij = np.array([-1,-2,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_1_2_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((0,1),(0,2), (0,2)), 'constant')[1:,2:,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_1_2_2); del faza_1_2_2

               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 26 ::: [212], [221] >>> 8 sites -------------------------------------------'''
           elif s == '212':
               dij = np.array([2,1,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza212 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(1,0), (2,0)), 'constant')[:-2,:-1,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza212); del faza212

           elif s == '21_2':
               dij = np.array([2,1,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza21_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(1,0), (0,2)), 'constant')[:-2,:-1,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza21_2); del faza21_2

           elif s == '2_12':
               dij = np.array([2,-1,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza2_12 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,1), (2,0)), 'constant')[:-2,1:,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza2_12); del faza2_12

           elif s == '2_1_2':
               dij = np.array([2,-1,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza2_1_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,1), (0,2)), 'constant')[:-2,1:,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza2_1_2); del faza2_1_2

           elif s == '221':
               dij = np.array([2,2,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza221 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(2,0), (1,0)), 'constant')[:-2,:-2,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza221); del faza221

           elif s == '22_1':
               dij = np.array([2,2,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza22_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(2,0), (0,1)), 'constant')[:-2,:-2,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza22_1); del faza22_1

           elif s == '2_21':
               dij = np.array([2,-2,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza2_21 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,2), (1,0)), 'constant')[:-2,2:,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza2_21); del faza2_21

           elif s == '2_2_1':
               dij = np.array([2,-2,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza2_2_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza2_2_1); del faza2_2_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 27 ::: [_212], [_221] >>> 8 sites -------------------------------------------'''
           elif s == '_212':
               dij = np.array([-2,1,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_212 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(1,0), (2,0)), 'constant')[2:,:-1,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_212); del faza_212

           elif s == '_21_2':
               dij = np.array([-2,1,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_21_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(1,0), (0,2)), 'constant')[2:,:-1,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_21_2); del faza_21_2

           elif s == '_2_12':
               dij = np.array([-2,-1,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_2_12 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,1), (2,0)), 'constant')[2:,1:,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_2_12); del faza_2_12

           elif s == '_2_1_2':
               dij = np.array([-2,-1,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_2_1_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,1), (0,2)), 'constant')[2:,1:,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_2_1_2); del faza_2_1_2

           elif s == '_221':
               dij = np.array([-2,2,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_221 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(2,0), (1,0)), 'constant')[2:,:-2,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_221); del faza_221

           elif s == '_22_1':
               dij = np.array([-2,2,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_22_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(2,0), (0,1)), 'constant')[2:,:-2,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_22_1); del faza_22_1

           elif s == '_2_21':
               dij = np.array([-2,-2,1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_2_21 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,2), (1,0)), 'constant')[2:,2:,:-1]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_2_21); del faza_2_21

           elif s == '_2_2_1':
               dij = np.array([-2,-2,-1])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_2_2_1 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_8th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,2), (0,1)), 'constant')[2:,2:,1:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_2_2_1); del faza_2_2_1
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               
               ''' ----------------------------------------- GROUP 28 ::: [222] >>> 4 sites -------------------------------------------'''
           elif s == '222':
               dij = np.array([2,2,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza222 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_9th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(2,0), (2,0)), 'constant')[:-2,:-2,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza222); del faza222

           elif s == '22_2':
               dij = np.array([2,2,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza22_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_9th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(2,0), (0,2)), 'constant')[:-2,:-2,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza22_2); del faza22_2

           elif s == '2_22':
               dij = np.array([2,-2,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza2_22 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_9th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,2), (2,0)), 'constant')[:-2,2:,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza2_22); del faza2_22

           elif s == '2_2_2':
               dij = np.array([2,-2,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza2_2_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_9th(attr11, attr12, attr15) )&(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza2_2_2); del faza2_2_2
               ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''        
               # *** UNDER CONattr6TRUCTION ! ***
               ''' ----------------------------------------- GROUP 29 ::: [_222] >>> 4 sites -------------------------------------------'''
           elif s == '_222':
               dij = np.array([-2,2,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_222 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_9th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(2,0), (2,0)), 'constant')[2:,:-2,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_222); del faza_222

           elif s == '_22_2':
               dij = np.array([-2,2,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_22_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_9th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(2,0), (0,2)), 'constant')[2:,:-2,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_22_2); del faza_22_2

           elif s == '_2_22':
               dij = np.array([-2,-2,2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_2_22 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_9th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,2), (2,0)), 'constant')[2:,2:,:-2]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_2_22); del faza_2_22

           elif s == '_2_2_2':
               dij = np.array([-2,-2,-2])
               wij = attr13(attr4[grain]['oi'], dij)
               lij = wij*attr5*attr6[s]
               faza_2_2_2 = np.where((attr8 < attr9)&(attr3==0)& (attr7==0)& (attr10==0) &(lij>=Dsr_9th(attr11, attr12, attr15) )&(np.pad(attr3,((0,2),(0,2), (0,2)), 'constant')[2:,2:,2:]==grain)& (                                                 
                                                                           #wrong(np.pad(attr3,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==grain)      |      # 1_1_2                
                                                                           #wrong(np.pad(attr3,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==grain)      |      # 1_2_1
                                                                           #wrong(np.pad(attr3,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==grain)         |      # 0_2_2
                                                                           #wrong(np.pad(attr3,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==grain)             # 2_2_2
                                                True), grain, attr3)
               seznam_premikov.append(faza_2_2_2); del faza_2_2_2

           attr3 = seznam_premikov[0]
           
       attr14 = attr6[s]
       #return attr1, attr3, attr4, attr5, attr6, attr7, attr8, attr9, attr10, attr11, attr12, attr13, attr14, attr15
       return attr3, attr14
       #q.put(attr3, attr14)
        

# attr1= GDSM, attr2= Selection, attr3= faza, attr4= asc, attr5= vg, attr6= S, attr7= taula, attr8= T_next, attr9= T, attr10= likvid, attr11= R, attr12= cell, attr13= W, attr14= cas, attr15= random_distances 
       


