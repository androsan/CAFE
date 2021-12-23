import os
import sys
import shutil
import numpy as np
from scipy.ndimage.interpolation import map_coordinates
import matplotlib.pyplot as plt
import random, math, time
plt.ion()

#---------------------------------------------------------------------------------------------  F U N C T I O N S ---------------------------------------------#
def m(pot, stevilo):
   q=np.load((pot+'salome_{}.npy').format(stevilo))
   return q

def N_Int_Files(sf, tf1, tf2, first):
   if first:
      SF_files = sf-1
      TF1_files = tf1 + (SF_files-2)*(tf1-1)
      TF2_files = tf2 + (TF1_files-2)*(tf2-1)
   else:
      SF_files = sf
      TF1_files = tf1 + (SF_files-2)*(tf1-1)
      TF2_files = tf2 + (TF1_files-1)*(tf2-1)
   return SF_files, TF1_files, TF2_files
   
def Get_Sorted_File_List():
   _path_, dirs, files = next(os.walk(ath))
   nums = sorted([int(a[7:].split('.')[0])for a in files])
   strings = [('salome_{}.npy').format(b) for b in nums]
   len_strings = len(files)
   return len_strings, strings

def Merge_2_Space_Interpolated_NPYs(YP, TM, file_0, file_1):
   global Z,X,Y
   Snpy_0 = np.load(ath+file_0);Snpy_1 = np.load(ath+file_1)
   mat_4D = np.array([Snpy_0, Snpy_1]); del Snpy_0; del Snpy_1
   #file =    ('salome_({0},{1}).npy').format(N, M)     #  two .npy files merged
   #np.save(ath+file, mat_4D)                                              #  saving merged .npy files
   Z = mat_4D.shape[1]; X = mat_4D.shape[2]; Y = mat_4D.shape[3]
   return mat_4D

def Matrix_Interpolation(A):
   B_tuple = (time_factor, Z, X, Y)
   new_dims = []
   for original_length, new_length in zip(A.shape, B_tuple):
      new_dims.append(np.linspace(0, original_length-1, new_length))
   coords = np.meshgrid(*new_dims, indexing='ij')
   B = map_coordinates(A, coords, order=1)
   return B

def Matrix_Time_Interpolation(YP, TM, cti, mapa_super, stp2):     # time interpolation in given ../YP/TM folder over all space interpolated npy files
   global ath, TSC, counter, file_count, files
   ath = PATH+YP+'/'+TM+'/'+mapa_super+'/'
   mapa = 'time_factor_'+str(time_factor)
   print('Working in folder ../'+YP+'/'+TM+'/'+mapa_super+'/'+mapa)
   if not os.path.isdir(ath+mapa):
      os.mkdir(ath+mapa)
   file_count, files = Get_Sorted_File_List()
   TSC = []; tsc_counter = -1
   
   for i in range(file_count - 1):
     if i==0:
        tsc = [j for j in range(last, last + time_factor)] 
        TSC.append(tsc)
     else:
        tsc = [j for j in range(tsc[-1]+1, tsc[-1]+time_factor)] 
        TSC.append(tsc)
        tsc_counter -= 1

   start_total = time.time() 
   for ii in range(file_count - 1):
      start = time.time(); counter = 0
      Matrix_4D = Merge_2_Space_Interpolated_NPYs(YP, TM, files[ii], files[ii+1])
      inp_mat = Matrix_4D[:2, :, :, :]    
      out_mat = Matrix_Interpolation(inp_mat); del inp_mat
      
      if ii==0:
          od=None
      else:
          od=1
      for j in out_mat[od:]:
          if tm_index >= cti:
             np.save(ath+mapa+'/salome_'+str(TSC[ii][counter])+'.npy', j)
          counter+=1
      end = time.time()
 
   # for connecting indivudal TR folders, i.e. last (pre-last, etc.) .npy file from working directory (say TR0) puts into next (TR1), which is created if not yet..
   try: 
       ath_next = PATH+YP+'/'+next_tm+'/'+mapa_super+'/'   
       np.save(ath_next+mapa+'/salome_'+str(TSC[ii][counter-1])+'.npy', out_mat[-1])# last .npy file
       if stp2:
          np.save(ath_next+mapa+'/salome_'+str(TSC[ii][counter-2])+'.npy', out_mat[-2])# pre-last .npy file
   except FileNotFoundError:
       os.mkdir(ath_next+mapa)
       np.save(ath_next+mapa+'/salome_'+str(TSC[ii][counter-1])+'.npy', out_mat[-1])# last .npy file
       if stp2:
          np.save(ath_next+mapa+'/salome_'+str(TSC[ii][counter-2])+'.npy', out_mat[-2])# pre-last .npy file
   print('Total computing time =  ',round(time.time()-start_total, 3),'  seconds.'); print()


def Delete_and_Move_files(YP, TM, F1, F2):
   print('Deleting and moving the files ..'); print()
   # 1st STEP:: Deleting only space interpolated .npyfiles from ..YP/TM folder
   mydir = PATH+YP+'/'+TM+'/'
   filelist = [f for f in os.listdir(mydir) if f.endswith(".npy") ]
   for f in filelist:
       os.remove(os.path.join(mydir, f))

   # 2nd STEP:: Moving final time interpolated files from ..YP/TM/../../../time_factor_.. into ..YP/TM/
   source_dir = mydir+F1+'/'+F2
   file_names = os.listdir(source_dir)
   for file_name in file_names[1:]:
       shutil.move(os.path.join(source_dir, file_name), mydir)

   # 3rd STEP:: Deleting ..YP/TM/time_factor_.. folder
   shutil.rmtree(mydir+F1)

   

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------#

case =         'SLM_2D_Source'
subcase =    '0002/'
subpath =    'INTER  time=1, space=8  Z[0-9], X[15-27], Y[12-97], 1500°C, N=12/'
#subpath =     'Testing_the_Toycode/'

PATH =       'C:/sm-2018-w64-0-3/WORK/'+case+'_Files/post processing database/'+subcase+subpath


# Creation of YPTM dictionary of YP folders and TM subfolders
yp_dirs = []
for subdir in os.walk(PATH):
   yp_dirs.append(subdir)
YP_dirlist = yp_dirs[0][1]

YPTM = {}
for ypdir in YP_dirlist:
   tm_dirs = []
   for tmdir in os.walk(PATH+'/'+ypdir):
      tm_dirs.append(tmdir)
   YPTM[ypdir] = tm_dirs[0][1]

# i.e.  YPTM ={'YP0  [12,20]': ['TR1  [1,2]', 'TR2  [2,3]'], 'YP1  [19,27]': ['TR1  [1,2]', 'TR2  [2,3]']}


# H contents :::  'SF' (special folder of interpolated files from previous step), 'TF' (time factor of current step),
#                       'Last' (dictionary of starting indices of salome_{}.npy files of given time folders)                        


H = {1:{'SF':'', 'TF':24, 'Last': {'TR0': None,        'TR1': 0,            'TR2': 139,           'TR3': 301,
                                                  'TR4': 463,          'TR5': 625,         'TR6': 787,           'TR7': 949,
                                                  'TR8': 1111,        'TR9': 1273,        'TR10': 1435,       'TR11': 1597,
                                                  'TR12': 1759,      'TR13': 1921,       'TR14': 2083,       'TR15': 2245,
                                   
                                                 }
              },  #  STEP 1 time interpolation
     
         2: {'SF':'time_factor_24', 'TF':3, 'Last': {'TR0': None,          'TR1': 0,               'TR2': 277,             'TR3': 602,
                                                                         'TR4': 927,            'TR5': 1252,          'TR6': 1577,           'TR7': 1902,
                                                                         'TR8': 2227,          'TR9': 2552,          'TR10': 2877,         'TR11': 3202,
                                                                         'TR12': 3527,         'TR13': 3852,        'TR14': 4177,         'TR15': 4502,
                                                                        }

              }   #  STEP 2 time interpolation
          } 

segments = ['YP0', 'YP1', 'YP2', 'YP3', 'YP4', 'YP5', 'YP6', 'YP7', 'YP8', 'YP9', 'YP10', 'YP11']


step_1 = True
step_2 = True

critical_time_index = 0


delete =   False

''' ................................................................. E  X  E  C  U  T  I  O  N .................................................'''
''' Iteration over given YP folders
   and TM subfolders in YPTM database '''

for yp in YPTM:
   
   if yp[:3] in segments[ 7 : 8 ]:   # PICK  YP  segments from segments list
      tm_list = sorted(YPTM[yp], key=lambda x: int(x.split('  ')[0][2:]))
      
      for tm in tm_list  [ 7 : 13 ]:   # PICK  TR time frames from corresponding YP folder  
         tm_index = tm_list.index(tm)
         if step_1:
            
            STEP = 1
            print(20*'*'+' STEP = '+str(STEP)+' '+20*'*'); print()
            super_folder = H[STEP]['SF']; time_factor = H[STEP]['TF']
            last = H[STEP]['Last'][tm.split('  ')[0]]; next_tm = tm_list[tm_index+1]
            Matrix_Time_Interpolation(yp, tm, critical_time_index, super_folder, False)
            print()
            #................... sys.exit()#   >>>  S . T . O . P .   line  STOP the execution <<< EXIT() ----- !! !  !   !
         if step_2:
            STEP = 2
            print(20*'*'+' STEP = '+str(STEP)+' '+20*'*'); print()
            super_folder = H[STEP]['SF']; time_factor = H[STEP]['TF']
            last = H[STEP]['Last'][tm.split('  ')[0]]; next_tm = tm_list[tm_index+1]
            Matrix_Time_Interpolation(yp, tm, critical_time_index, super_folder, True)

         if delete:
            Delete_and_Move_files(yp, tm, 'time_factor_10', 'time_factor_6')
         
         print(30*'-'+' DONE ! '+30*'-'); print()











