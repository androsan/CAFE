import numpy as np
import math

def NN(func1, attr1):
    for negind in attr1:
        func1(negind)
    
def MM(attr1, attr2, attr3, attr4, attr5, attr6):
    for s in attr1:
        for negind in attr2:
            ''' ----------------------------------------- GROUP 1 ::: [001], [010] >>> 4 sites -------------------------------------'''
            if s == '001':
                if negind == -1:
                    cas001 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas001; del cas001
                else:
                    cas001 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas001; del cas001

            elif s == '00_1':
                if negind == -1:
                    cas00_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas00_1; del cas00_1
                else:
                    cas00_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,0), (0,1)), 'constant')[:,:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas00_1; del cas00_1

            elif s == '010':
                if negind == -1:
                    cas010 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas010; del cas010
                else:
                    cas010 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(1,0), (0,0)), 'constant')[:,:-1,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas010; del cas010

            elif s == '0_10':
                if negind == -1:
                    cas0_10 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas0_10; del cas0_10
                else:
                    cas0_10 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,1), (0,0)), 'constant')[:,1:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas0_10; del cas0_10
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 2 ::: [100] >>> 1 site -------------------------------------------'''
            elif s == '100':
                if negind == -1:
                    cas100 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas100; del cas100
                else:
                    cas100 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,0), (0,0)), 'constant')[:-1,:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas100; del cas100
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 3 ::: [_100] >>> 1 site -------------------------------------------'''
            elif s == '_100':
                if negind == -1:
                    cas_100 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_100; del cas_100
                else:
                    cas_100 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,0), (0,0)), 'constant')[1:,:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_100; del cas_100
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 4 ::: [011] >>> 4 sites -------------------------------------------'''
            elif s == '011':
                if negind == -1:
                    cas011 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas011; del cas011
                else:
                    cas011 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(1,0), (1,0)), 'constant')[:,:-1,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas011; del cas011

            elif s == '01_1':
                if negind == -1:
                    cas01_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas01_1; del cas01_1
                else:
                    cas01_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(1,0), (0,1)), 'constant')[:,:-1,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas01_1; del cas01_1

            elif s == '0_11':
                if negind == -1:
                    cas0_11 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas0_11; del cas0_11
                else:
                    cas0_11 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,1), (1,0)), 'constant')[:,1:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas0_11; del cas0_11

            elif s == '0_1_1':
                if negind == -1:
                    cas0_1_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas0_1_1; del cas0_1_1
                else:
                    cas0_1_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,1), (0,1)), 'constant')[:,1:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas0_1_1; del cas0_1_1
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 5 ::: [101], [110] >>> 4 sites -------------------------------------------'''
            elif s == '101':
                if negind == -1:
                    cas101 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas101; del cas101
                else:
                    cas101 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,0), (1,0)), 'constant')[:-1,:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas101; del cas101

            elif s == '10_1':
                if negind == -1:
                    cas10_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas10_1; del cas10_1
                else:
                    cas10_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,0), (0,1)), 'constant')[:-1,:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas10_1; del cas10_1

            elif s == '110':
                if negind == -1:
                    cas110 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas110; del cas110
                else:
                    cas110 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(1,0), (0,0)), 'constant')[:-1,:-1,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas110; del cas110

            elif s == '1_10':
                if negind == -1:
                    cas1_10 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas1_10; del cas1_10
                else:
                    cas1_10 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,1), (0,0)), 'constant')[:-1,1:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas1_10; del cas1_10
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 6 ::: [_101], [_110] >>> 4 sites -------------------------------------------'''
            elif s == '_101':
                if negind == -1:
                    cas_101 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,0), (1,0)), 'constant')[1:,:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas_101; del cas_101
                else:
                    cas_101 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,0), (1,0)), 'constant')[1:,:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_101; del cas_101

            elif s == '_10_1':
                if negind == -1:
                    cas_10_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,0), (0,1)), 'constant')[1:,:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_10_1; del cas_10_1
                else:
                    cas_10_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,0), (0,1)), 'constant')[1:,:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_10_1; del cas_10_1

            elif s == '_110':
                if negind == -1:
                    cas_110 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(1,0), (0,0)), 'constant')[1:,:-1,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_110; del cas_110
                else:
                    cas_110 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(1,0), (0,0)), 'constant')[1:,:-1,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_110; del cas_110

            elif s == '_1_10':
                if negind == -1:
                    cas_1_10 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,1), (0,0)), 'constant')[1:,1:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_1_10; del cas_1_10
                else:
                    cas_1_10 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,1), (0,0)), 'constant')[1:,1:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_1_10; del cas_1_10
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 7 ::: [111] >>> 4 sites -------------------------------------------'''
            elif s == '111':
                if negind == -1:
                    cas111 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas111; del cas111
                else:
                    cas111 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(1,0), (1,0)), 'constant')[:-1,:-1,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas111; del cas111

            elif s == '11_1':
                if negind == -1:
                    cas11_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas11_1; del cas11_1
                else:
                    cas11_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(1,0), (0,1)), 'constant')[:-1,:-1,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas11_1; del cas11_1

            elif s == '1_11':
                if negind == -1:
                    cas1_11 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas1_11; del cas1_11
                else:
                    cas1_11 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,1), (1,0)), 'constant')[:-1,1:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas1_11; del cas1_11

            elif s == '1_1_1':
                if negind == -1:
                    cas1_1_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas1_1_1; del cas1_1_1
                else:
                    cas1_1_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,1), (0,1)), 'constant')[:-1,1:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas1_1_1; del cas1_1_1
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 8 ::: [_111] >>> 4 sites -------------------------------------------'''
            elif s == '_111':
                if negind == -1:
                    cas_111 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas_111; del cas_111
                else:
                    cas_111 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(1,0), (1,0)), 'constant')[1:,:-1,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_111; del cas_111

            elif s == '_11_1':
                if negind == -1:
                    cas_11_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_11_1; del cas_11_1
                else:
                    cas_11_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(1,0), (0,1)), 'constant')[1:,:-1,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_11_1; del cas_11_1

            elif s == '_1_11':
                if negind == -1:
                    cas_1_11 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas_1_11; del cas_1_11
                else:
                    cas_1_11 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,1), (1,0)), 'constant')[1:,1:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_1_11; del cas_1_11

            elif s == '_1_1_1':
                if negind == -1:
                    cas_1_1_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_1_1_1; del cas_1_1_1
                else:
                    cas_1_1_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,1), (0,1)), 'constant')[1:,1:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_1_1_1; del cas_1_1_1
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 9 ::: [012], [021] >>> 8 sites -------------------------------------------'''
            elif s == '012':
               if negind == -1:
                    cas012 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas012; del cas012
               else:
                    cas012 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(1,0), (2,0)), 'constant')[:,:-1,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas012; del cas012

            elif s == '01_2':
               if negind == -1:
                    cas01_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas01_2; del cas01_2
               else:
                    cas01_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(1,0), (0,2)), 'constant')[:,:-1,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas01_2; del cas01_2

            elif s == '0_12':
               if negind == -1:
                    cas0_12 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas0_12; del cas0_12
               else:
                    cas0_12 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,1), (2,0)), 'constant')[:,1:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas0_12; del cas0_12

            elif s == '0_1_2':
               if negind == -1:
                    cas0_1_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas0_1_2; del cas0_1_2
               else:
                    cas0_1_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,1), (0,2)), 'constant')[:,1:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas0_1_2; del cas0_1_2

            elif s == '021':
               if negind == -1:
                    cas021 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas021; del cas021
               else:
                    cas021 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(2,0), (1,0)), 'constant')[:,:-2,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas021; del cas021

            elif s == '02_1':
               if negind == -1:
                    cas02_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas02_1; del cas02_1
               else:
                    cas02_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(2,0), (0,1)), 'constant')[:,:-2,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas02_1; del cas02_1

            elif s == '0_21':
               if negind == -1:
                    cas0_21 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas0_21; del cas0_21
               else:
                    cas0_21 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,2), (1,0)), 'constant')[:,2:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas0_21; del cas0_21

            elif s == '0_2_1':
               if negind == -1:
                    cas0_2_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas0_2_1; del cas0_2_1
               else:
                    cas0_2_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,2), (0,1)), 'constant')[:,2:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas0_2_1; del cas0_2_1
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 10 ::: [102], [120] >>> 4 sites -------------------------------------------'''
            elif s == '102':
               if negind == -1:
                    cas102 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas102; del cas102
               else:
                    cas102 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,0), (2,0)), 'constant')[:-1,:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas102; del cas102

            elif s == '10_2':
               if negind == -1:
                    cas10_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas10_2; del cas10_2
               else:
                    cas10_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,0), (0,2)), 'constant')[:-1,:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas10_2; del cas10_2

            elif s == '120':
               if negind == -1:
                    cas120 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas120; del cas120
               else:
                    cas120 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(2,0), (0,0)), 'constant')[:-1,:-2,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas120; del cas120

            elif s == '1_20':
               if negind == -1:
                    cas1_20 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas1_20; del cas1_20
               else:
                    cas1_20 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,2), (0,0)), 'constant')[:-1,2:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas1_20; del cas1_20
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 11 ::: [_102], [_120] >>> 4 sites -------------------------------------------'''
            elif s == '_102':
               if negind == -1:
                    cas_102 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,0), (2,0)), 'constant')[1:,:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas_102; del cas_102
               else:
                    cas_102 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,0), (2,0)), 'constant')[1:,:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_102; del cas_102

            elif s == '_10_2':
               if negind == -1:
                    cas_10_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,0), (0,2)), 'constant')[1:,:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_10_2; del cas_10_2
               else:
                    cas_10_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,0), (0,2)), 'constant')[1:,:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_10_2; del cas_10_2

            elif s == '_120':
               if negind == -1:
                    cas_120 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(2,0), (0,0)), 'constant')[1:,:-2,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_120; del cas_120
               else:
                    cas_120 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(2,0), (0,0)), 'constant')[1:,:-2,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_120; del cas_120

            elif s == '_1_20':
               if negind == -1:
                    cas_1_20 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,2), (0,0)), 'constant')[1:,2:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_1_20; del cas_1_20
               else:
                    cas_1_20 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,2), (0,0)), 'constant')[1:,2:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_1_20; del cas_1_20
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 12 ::: [201], [210] >>> 4 sites -------------------------------------------'''
            elif s == '201':
               if negind == -1:
                    cas201 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,0), (1,0)), 'constant')[:-2,:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas201; del cas201
               else:
                    cas201 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,0), (1,0)), 'constant')[:-2,:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas201; del cas201

            elif s == '20_1':
               if negind == -1:
                    cas20_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,0), (0,1)), 'constant')[:-2,:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas20_1; del cas20_1
               else:
                    cas20_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,0), (0,1)), 'constant')[:-2,:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas20_1; del cas20_1

            elif s == '210':
               if negind == -1:
                    cas210 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(1,0), (0,0)), 'constant')[:-2,:-1,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas210; del cas210
               else:
                    cas210 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(1,0), (0,0)), 'constant')[:-2,:-1,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas210; del cas210

            elif s == '2_10':
               if negind == -1:
                    cas2_10 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,1), (0,0)), 'constant')[:-2,1:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas2_10; del cas2_10
               else:
                    cas2_10 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,1), (0,0)), 'constant')[:-2,1:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas2_10; del cas2_10
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 13 ::: [_201], [_210] >>> 4 sites -------------------------------------------'''
            elif s == '_201':
               if negind == -1:
                    cas_201 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,0), (1,0)), 'constant')[2:,:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas_201; del cas_201
               else:
                    cas_201 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,0), (1,0)), 'constant')[2:,:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_201; del cas_201

            elif s == '_20_1':
               if negind == -1:
                    cas_20_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,0), (0,1)), 'constant')[2:,:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_20_1; del cas_20_1
               else:
                    cas_20_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,0), (0,1)), 'constant')[2:,:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_20_1; del cas_20_1

            elif s == '_210':
               if negind == -1:
                    cas_210 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(1,0), (0,0)), 'constant')[2:,:-1,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_210; del cas_210
               else:
                    cas_210 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(1,0), (0,0)), 'constant')[2:,:-1,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_210; del cas_210

            elif s == '_2_10':
               if negind == -1:
                    cas_2_10 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,1), (0,0)), 'constant')[2:,1:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_2_10; del cas_2_10
               else:
                    cas_2_10 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,1), (0,0)), 'constant')[2:,1:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_2_10; del cas_2_10
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 14 ::: [002], [020] >>> 4 sites -------------------------------------------'''
            elif s == '002':
                if negind == -1:
                    cas002 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas002; del cas002
                else:
                    cas002 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,0), (2,0)), 'constant')[:,:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas002; del cas002

            elif s == '00_2':
                if negind == -1:
                    cas00_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas00_2; del cas00_2
                else:
                    cas00_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,0), (0,2)), 'constant')[:,:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas00_2; del cas00_2

            elif s == '020':
                if negind == -1:
                    cas020 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas020; del cas020
                else:
                    cas020 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(2,0), (0,0)), 'constant')[:,:-2,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas020; del cas020

            elif s == '0_20':
                if negind == -1:
                    cas0_20 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas0_20; del cas0_20
                else:
                    cas0_20 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,2), (0,0)), 'constant')[:,2:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas0_20; del cas0_20
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 15 ::: [200] >>> 1 site -------------------------------------------'''
            elif s == '200':
                if negind == -1:
                    cas200 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,0), (0,0)), 'constant')[:-2,:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas200; del cas200
                else:
                    cas200 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,0), (0,0)), 'constant')[:-2,:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas200; del cas200
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 16 ::: [_200] >>> 1 site -------------------------------------------'''
            elif s == '_200':
                if negind == -1:
                    cas_200 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,0), (0,0)), 'constant')[2:,:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_200; del cas_200
                else:
                    cas_200 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,0), (0,0)), 'constant')[2:,:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_200; del cas_200
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 17 ::: [022] >>> 4 sites -------------------------------------------'''
            elif s == '022':
               if negind == -1:
                    cas022 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas022; del cas022
               else:
                    cas022 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(2,0), (2,0)), 'constant')[:,:-2,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas022; del cas022

            elif s == '02_2':
               if negind == -1:
                    cas02_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas02_2; del cas02_2
               else:
                    cas02_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(2,0), (0,2)), 'constant')[:,:-2,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas02_2; del cas02_2

            elif s == '0_22':
               if negind == -1:
                    cas0_22 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas0_22; del cas0_22
               else:
                    cas0_22 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,2), (2,0)), 'constant')[:,2:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas0_22; del cas0_22

            elif s == '0_2_2':
               if negind == -1:
                    cas0_2_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas0_2_2; del cas0_2_2
               else:
                    cas0_2_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,0),(0,2), (0,2)), 'constant')[:,2:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas0_2_2; del cas0_2_2
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 18 ::: [202], [220] >>> 4 sites -------------------------------------------'''
            elif s == '202':
               if negind == -1:
                    cas202 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,0), (2,0)), 'constant')[:-2,:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas202; del cas202
               else:
                    cas202 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,0), (2,0)), 'constant')[:-2,:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas202; del cas202

            elif s == '20_2':
               if negind == -1:
                    cas20_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,0), (0,2)), 'constant')[:-2,:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas20_2; del cas20_2
               else:
                    cas20_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,0), (0,2)), 'constant')[:-2,:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas20_2; del cas20_2

            elif s == '220':
               if negind == -1:
                    cas220 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(2,0), (0,0)), 'constant')[:-2,:-2,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas220; del cas220
               else:
                    cas220 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(2,0), (0,0)), 'constant')[:-2,:-2,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas220; del cas220

            elif s == '2_20':
               if negind == -1:
                    cas2_20 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas2_20; del cas2_20
               else:
                    cas2_20 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,2), (0,0)), 'constant')[:-2,2:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas2_20; del cas2_20
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 19 ::: [_202], [_220] >>> 4 sites -------------------------------------------'''
            elif s == '_202':
               if negind == -1:
                    cas_202 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,0), (2,0)), 'constant')[2:,:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas_202; del cas_202
               else:
                    cas_202 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,0), (2,0)), 'constant')[2:,:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_202; del cas_202

            elif s == '_20_2':
               if negind == -1:
                    cas_20_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,0), (0,2)), 'constant')[2:,:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_20_2; del cas_20_2
               else:
                    cas_20_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,0), (0,2)), 'constant')[2:,:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_20_2; del cas_20_2

            elif s == '_220':
               if negind == -1:
                    cas_220 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(2,0), (0,0)), 'constant')[2:,:-2,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_220; del cas_220
               else:
                    cas_220 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(2,0), (0,0)), 'constant')[2:,:-2,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_220; del cas_220

            elif s == '_2_20':
               if negind == -1:
                    cas_2_20 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,2), (0,0)), 'constant')[2:,2:,:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_2_20; del cas_2_20
               else:
                    cas_2_20 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,2), (0,0)), 'constant')[2:,2:,:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_2_20; del cas_2_20
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 20 ::: [112], [121] >>> 8 sites -------------------------------------------'''
            elif s == '112':
               if negind == -1:
                    cas112 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas112; del cas112
               else:
                    cas112 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(1,0), (2,0)), 'constant')[:-1,:-1,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas112; del cas112

            elif s == '11_2':
               if negind == -1:
                    cas11_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas11_2; del cas11_2
               else:
                    cas11_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(1,0), (0,2)), 'constant')[:-1,:-1,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas11_2; del cas11_2

            elif s == '1_12':
               if negind == -1:
                    cas1_12 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas1_12; del cas1_12
               else:
                    cas1_12 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,1), (2,0)), 'constant')[:-1,1:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas1_12; del cas1_12

            elif s == '1_1_2':
               if negind == -1:
                    cas1_1_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas1_1_2; del cas1_1_2
               else:
                    cas1_1_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,1), (0,2)), 'constant')[:-1,1:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas1_1_2; del cas1_1_2

            elif s == '121':
               if negind == -1:
                    cas121 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas121; del cas121
               else:
                    cas121 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(2,0), (1,0)), 'constant')[:-1,:-2,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas121; del cas121

            elif s == '12_1':
               if negind == -1:
                    cas12_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas12_1; del cas12_1
               else:
                    cas12_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(2,0), (0,1)), 'constant')[:-1,:-2,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas12_1; del cas12_1

            elif s == '1_21':
               if negind == -1:
                    cas1_21 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas1_21; del cas1_21
               else:
                    cas1_21 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,2), (1,0)), 'constant')[:-1,2:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas1_21; del cas1_21

            elif s == '1_2_1':
               if negind == -1:
                    cas1_2_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas1_2_1; del cas1_2_1
               else:
                    cas1_2_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,2), (0,1)), 'constant')[:-1,2:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas1_2_1; del cas1_2_1
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 21 ::: [_112], [_121] >>> 8 sites -------------------------------------------'''
            elif s == '_112':
               if negind == -1:
                    cas_112 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(1,0), (2,0)), 'constant')[1:,:-1,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas_112; del cas_112
               else:
                    cas_112 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(1,0), (2,0)), 'constant')[1:,:-1,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_112; del cas_112

            elif s == '_11_2':
               if negind == -1:
                    cas_11_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(1,0), (0,2)), 'constant')[1:,:-1,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_11_2; del cas_11_2
               else:
                    cas_11_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(1,0), (0,2)), 'constant')[1:,:-1,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_11_2; del cas_11_2

            elif s == '_1_12':
               if negind == -1:
                    cas_1_12 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,1), (2,0)), 'constant')[1:,1:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas_1_12; del cas_1_12
               else:
                    cas_1_12 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,1), (2,0)), 'constant')[1:,1:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_1_12; del cas_1_12

            elif s == '_1_1_2':
               if negind == -1:
                    cas_1_1_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,1), (0,2)), 'constant')[1:,1:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_1_1_2; del cas_1_1_2
               else:
                    cas_1_1_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,1), (0,2)), 'constant')[1:,1:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_1_1_2; del cas_1_1_2

            elif s == '_121':
               if negind == -1:
                    cas_121 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(2,0), (1,0)), 'constant')[1:,:-2,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas_121; del cas_121
               else:
                    cas_121 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(2,0), (1,0)), 'constant')[1:,:-2,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_121; del cas_121

            elif s == '_12_1':
               if negind == -1:
                    cas_12_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(2,0), (0,1)), 'constant')[1:,:-2,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_12_1; del cas_12_1
               else:
                    cas_12_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(2,0), (0,1)), 'constant')[1:,:-2,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_12_1; del cas_12_1

            elif s == '_1_21':
               if negind == -1:
                    cas_1_21 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,2), (1,0)), 'constant')[1:,2:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas_1_21; del cas_1_21
               else:
                    cas_1_21 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,2), (1,0)), 'constant')[1:,2:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_1_21; del cas_1_21

            elif s == '_1_2_1':
               if negind == -1:
                    cas_1_2_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,2), (0,1)), 'constant')[1:,2:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_1_2_1; del cas_1_2_1
               else:
                    cas_1_2_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,2), (0,1)), 'constant')[1:,2:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_1_2_1; del cas_1_2_1
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 22 ::: [211] >>> 4 sites -------------------------------------------'''
            elif s == '211':
               if negind == -1:
                    cas211 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(1,0), (1,0)), 'constant')[:-2,:-1,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas211; del cas211
               else:
                    cas211 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(1,0), (1,0)), 'constant')[:-2,:-1,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas211; del cas211

            elif s == '21_1':
               if negind == -1:
                    cas21_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(1,0), (0,1)), 'constant')[:-2,:-1,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas21_1; del cas21_1
               else:
                    cas21_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(1,0), (0,1)), 'constant')[:-2,:-1,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas21_1; del cas21_1

            elif s == '2_11':
               if negind == -1:
                    cas2_11 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,1), (1,0)), 'constant')[:-2,1:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas2_11; del cas2_11
               else:
                    cas2_11 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,1), (1,0)), 'constant')[:-2,1:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas2_11; del cas2_11

            elif s == '2_1_1':
               if negind == -1:
                    cas2_1_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,1), (0,1)), 'constant')[:-2,1:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas2_1_1; del cas2_1_1
               else:
                    cas2_1_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,1), (0,1)), 'constant')[:-2,1:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas2_1_1; del cas2_1_1
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 23 ::: [_211] >>> 4 sites -------------------------------------------'''
            elif s == '_211':
               if negind == -1:
                    cas_211 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(1,0), (1,0)), 'constant')[2:,:-1,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas_211; del cas_211
               else:
                    cas_211 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(1,0), (1,0)), 'constant')[2:,:-1,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_211; del cas_211

            elif s == '_21_1':
               if negind == -1:
                    cas_21_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(1,0), (0,1)), 'constant')[2:,:-1,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_21_1; del cas_21_1
               else:
                    cas_21_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(1,0), (0,1)), 'constant')[2:,:-1,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_21_1; del cas_21_1

            elif s == '_2_11':
               if negind == -1:
                    cas_2_11 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,1), (1,0)), 'constant')[2:,1:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas_2_11; del cas_2_11
               else:
                    cas_2_11 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,1), (1,0)), 'constant')[2:,1:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_2_11; del cas_2_11

            elif s == '_2_1_1':
               if negind == -1:
                    cas_2_1_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,1), (0,1)), 'constant')[2:,1:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_2_1_1; del cas_2_1_1
               else:
                    cas_2_1_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,1), (0,1)), 'constant')[2:,1:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_2_1_1; del cas_2_1_1
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 24 ::: [122] >>> 4 sites -------------------------------------------'''
            elif s == '122':
               if negind == -1:
                    cas122 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas122; del cas122
               else:
                    cas122 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(2,0), (2,0)), 'constant')[:-1,:-2,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas122; del cas122

            elif s == '12_2':
               if negind == -1:
                    cas12_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas12_2; del cas12_2
               else:
                    cas12_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(2,0), (0,2)), 'constant')[:-1,:-2,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas12_2; del cas12_2

            elif s == '1_22':
               if negind == -1:
                    cas1_22 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas1_22; del cas1_22
               else:
                    cas1_22 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,2), (2,0)), 'constant')[:-1,2:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas1_22; del cas1_22

            elif s == '1_2_2':
               if negind == -1:
                    cas1_2_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas1_2_2; del cas1_2_2
               else:
                    cas1_2_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((1,0),(0,2), (0,2)), 'constant')[:-1,2:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas1_2_2; del cas1_2_2
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 25 ::: [_122] >>> 4 sites -------------------------------------------'''
            elif s == '_122':
               if negind == -1:
                    cas_122 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(2,0), (2,0)), 'constant')[1:,:-2,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas_122; del cas_122
               else:
                    cas_122 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(2,0), (2,0)), 'constant')[1:,:-2,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_122; del cas_122

            elif s == '_12_2':
               if negind == -1:
                    cas_12_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(2,0), (0,2)), 'constant')[1:,:-2,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_12_2; del cas_12_2
               else:
                    cas_12_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(2,0), (0,2)), 'constant')[1:,:-2,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_12_2; del cas_12_2

            elif s == '_1_22':
               if negind == -1:
                    cas_1_22 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,2), (2,0)), 'constant')[1:,2:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas_1_22; del cas_1_22
               else:
                    cas_1_22 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,2), (2,0)), 'constant')[1:,2:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_1_22; del cas_1_22

            elif s == '_1_2_2':
               if negind == -1:
                    cas_1_2_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,1),(0,2), (0,2)), 'constant')[1:,2:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_1_2_2; del cas_1_2_2
               else:
                    cas_1_2_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,1),(0,2), (0,2)), 'constant')[1:,2:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_1_2_2; del cas_1_2_2
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 26 ::: [212], [221] >>> 8 sites -------------------------------------------'''
            elif s == '212':
               if negind == -1:
                    cas212 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(1,0), (2,0)), 'constant')[:-2,:-1,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas212; del cas212
               else:
                    cas212 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(1,0), (2,0)), 'constant')[:-2,:-1,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas212; del cas212

            elif s == '21_2':
               if negind == -1:
                    cas21_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(1,0), (0,2)), 'constant')[:-2,:-1,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas21_2; del cas21_2
               else:
                    cas21_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(1,0), (0,2)), 'constant')[:-2,:-1,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas21_2; del cas21_2

            elif s == '2_12':
               if negind == -1:
                    cas2_12 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,1), (2,0)), 'constant')[:-2,1:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas2_12; del cas2_12
               else:
                    cas2_12 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,1), (2,0)), 'constant')[:-2,1:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas2_12; del cas2_12

            elif s == '2_1_2':
               if negind == -1:
                    cas2_1_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,1), (0,2)), 'constant')[:-2,1:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas2_1_2; del cas2_1_2
               else:
                    cas2_1_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,1), (0,2)), 'constant')[:-2,1:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas2_1_2; del cas2_1_2

            elif s == '221':
               if negind == -1:
                    cas221 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(2,0), (1,0)), 'constant')[:-2,:-2,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas221; del cas221
               else:
                    cas221 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(2,0), (1,0)), 'constant')[:-2,:-2,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas221; del cas221

            elif s == '22_1':
               if negind == -1:
                    cas22_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(2,0), (0,1)), 'constant')[:-2,:-2,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas22_1; del cas22_1
               else:
                    cas22_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(2,0), (0,1)), 'constant')[:-2,:-2,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas22_1; del cas22_1

            elif s == '2_21':
               if negind == -1:
                    cas2_21 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,2), (1,0)), 'constant')[:-2,2:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas2_21; del cas2_21
               else:
                    cas2_21 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,2), (1,0)), 'constant')[:-2,2:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas2_21; del cas2_21

            elif s == '2_2_1':
               if negind == -1:
                    cas2_2_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas2_2_1; del cas2_2_1
               else:
                    cas2_2_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,2), (0,1)), 'constant')[:-2,2:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas2_2_1; del cas2_2_1
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 27 ::: [_212], [_221] >>> 8 sites -------------------------------------------'''
            elif s == '_212':
               if negind == -1:
                    cas_212 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(1,0), (2,0)), 'constant')[2:,:-1,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas_212; del cas_212
               else:
                    cas_212 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(1,0), (2,0)), 'constant')[2:,:-1,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_212; del cas_212

            elif s == '_21_2':
               if negind == -1:
                    cas_21_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(1,0), (0,2)), 'constant')[2:,:-1,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_21_2; del cas_21_2
               else:
                    cas_21_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(1,0), (0,2)), 'constant')[2:,:-1,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_21_2; del cas_21_2

            elif s == '_2_12':
               if negind == -1:
                    cas_2_12 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,1), (2,0)), 'constant')[2:,1:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas_2_12; del cas_2_12
               else:
                    cas_2_12 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,1), (2,0)), 'constant')[2:,1:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_2_12; del cas_2_12

            elif s == '_2_1_2':
               if negind == -1:
                    cas_2_1_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,1), (0,2)), 'constant')[2:,1:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_2_1_2; del cas_2_1_2
               else:
                    cas_2_1_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,1), (0,2)), 'constant')[2:,1:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_2_1_2; del cas_2_1_2

            elif s == '_221':
               if negind == -1:
                    cas_221 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(2,0), (1,0)), 'constant')[2:,:-2,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas_221; del cas_221
               else:
                    cas_221 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(2,0), (1,0)), 'constant')[2:,:-2,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_221; del cas_221

            elif s == '_22_1':
               if negind == -1:
                    cas_22_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(2,0), (0,1)), 'constant')[2:,:-2,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_22_1; del cas_22_1
               else:
                    cas_22_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(2,0), (0,1)), 'constant')[2:,:-2,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_22_1; del cas_22_1

            elif s == '_2_21':
               if negind == -1:
                    cas_2_21 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,2), (1,0)), 'constant')[2:,2:,:-1]==negind), attr2[negind], attr5)
                    attr6[s]=cas_2_21; del cas_2_21
               else:
                    cas_2_21 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,2), (1,0)), 'constant')[2:,2:,:-1]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_2_21; del cas_2_21

            elif s == '_2_2_1':
               if negind == -1:
                    cas_2_2_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,2), (0,1)), 'constant')[2:,2:,1:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_2_2_1; del cas_2_2_1
               else:
                    cas_2_2_1 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,2), (0,1)), 'constant')[2:,2:,1:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_2_2_1; del cas_2_2_1
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    
                    ''' ----------------------------------------- GROUP 28 ::: [222] >>> 4 sites -------------------------------------------'''
            elif s == '222':
               if negind == -1:
                    cas222 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(2,0), (2,0)), 'constant')[:-2,:-2,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas222; del cas222
               else:
                    cas222 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(2,0), (2,0)), 'constant')[:-2,:-2,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas222; del cas222

            elif s == '22_2':
               if negind == -1:
                    cas22_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(2,0), (0,2)), 'constant')[:-2,:-2,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas22_2; del cas22_2
               else:
                    cas22_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(2,0), (0,2)), 'constant')[:-2,:-2,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas22_2; del cas22_2

            elif s == '2_22':
               if negind == -1:
                    cas2_22 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,2), (2,0)), 'constant')[:-2,2:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas2_22; del cas2_22
               else:
                    cas2_22 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,2), (2,0)), 'constant')[:-2,2:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas2_22; del cas2_22

            elif s == '2_2_2':
               if negind == -1:
                    cas2_2_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas2_2_2; del cas2_2_2
               else:
                    cas2_2_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((2,0),(0,2), (0,2)), 'constant')[:-2,2:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas2_2_2; del cas2_2_2
                    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
                    # under construction
                    ''' ----------------------------------------- GROUP 29 ::: [_222] >>> 4 sites -------------------------------------------'''
            elif s == '_222':
               if negind == -1:
                    cas_222 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(2,0), (2,0)), 'constant')[2:,:-2,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas_222; del cas_222
               else:
                    cas_222 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(2,0), (2,0)), 'constant')[2:,:-2,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_222; del cas_222

            elif s == '_22_2':
               if negind == -1:
                    cas_22_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(2,0), (0,2)), 'constant')[2:,:-2,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_22_2; del cas_22_2
               else:
                    cas_22_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(2,0), (0,2)), 'constant')[2:,:-2,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_22_2; del cas_22_2

            elif s == '_2_22':
               if negind == -1:
                    cas_2_22 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,2), (2,0)), 'constant')[2:,2:,:-2]==negind), attr2[negind], attr5)
                    attr6[s]=cas_2_22; del cas_2_22
               else:
                    cas_2_22 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,2), (2,0)), 'constant')[2:,2:,:-2]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_2_22; del cas_2_22

            elif s == '_2_2_2':
               if negind == -1:
                    cas_2_2_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr5,((0,2),(0,2), (0,2)), 'constant')[2:,2:,2:]==negind), attr2[negind], attr5)
                    attr6[s]=cas_2_2_2; del cas_2_2_2
               else:
                    cas_2_2_2 = np.where((np.isin(attr3, attr4, invert=True)) & (np.pad(attr6[s],((0,2),(0,2), (0,2)), 'constant')[2:,2:,2:]==negind), attr2[negind], attr6[s])
                    attr6[s]=cas_2_2_2; del cas_2_2_2


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


def GG(attr1, attr2, attr3, attr4, attr5, attr6, attr7, attr8, attr9, attr10, attr11, attr12, attr13, attr14, attr15):
    # attr1= GDSM, attr2= Selection, attr3= faza, attr4= asc, attr5= vg, attr6= S, attr7= taula, attr8= T_next, attr9= T, attr10= likvid, attr11= R, attr12= cell, attr13= W, attr14= cas, attr15= random_distances 
    for choice in attr2:
       if attr1 == 'new':
           seznam_premikov= []           # NEW
       grain=choice[0]
       s=choice[1]
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

    return attr3, attr14
        
    


