import os
import numpy as np
from scipy.ndimage.interpolation import map_coordinates
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random, math, time
plt.ion()

# FE analysis result from SALOME

case =      'SLM_2D_Source'
subcase = '0002'

PATH = 'C:/sm-2018-w64-0-3/WORK/'+case+'_Files/post processing database/'+subcase+'/'


file =      case+'_'+subcase+'.npy'  #  .npy filename with FE analysis results (after meshio_postprocessing.py) 
Matrix_4D  =  np.load(PATH+file)
               
''' ****************************************************************************************************************************************'''
space_factor =     8                               # Factor of spatial resolution increase

FEM_time_factor =      1                      # FEM time factor to make FEM analysis shorter ---> for this factor , i.e. diminish number of time steps
extra_time_factor =     1                       # ::: special time factor ::: to catch the mesh dependency effect !! !  !


time_factor = FEM_time_factor * extra_time_factor
increase_tuple  = (time_factor * space_factor, space_factor, space_factor, space_factor)    # resolution increase factor by linear interpolation 
''' ****************************************************************************************************************************************'''
Time_Range =  ( 1 , 14   )

long_track =     True
N =                     12                 # Number of equally sized partitions of the track
yp =                  'YP7'               # Name of Y partition

total_time_range =(1, 16)    #  total time range (number of FEM time steps)to set the domain constrains

auto_reduce_domain_size = False           #  True to find the constrains and create mapa;    False to make interpolations of individual time steps
interpolation = True                                    #  False to make .npy file of each Matrix_4D time step;    True to make time & space interpolation (scale up) of individual time step


Tmelt_Celsius = 1500                                                                     #   Melting point; unit:  degrees C [deg. C] - for domain size reduction
Tmelt= Tmelt_Celsius + 273                                                            #   Melting point; unit:  KELVIN [K]

TR_list = [i for i in range(Time_Range[0], Time_Range[1]+1)]

''' domain limits ----> reduction of FE domain to size a bit larger than the molten track, where the crystallization occurs .. '''

def Domain_Size_Reduction(domain, threshold):
   global z_min, z_max, x_min, x_max, y_min, y_max
   ax2=np.where(np.any(domain > threshold, axis=2))
   ax1=np.where(np.any(domain > threshold, axis=1))
   z_min = np.min(ax2[0]);  z_max = np.max(ax2[0])+1
   x_min = np.min(ax2[1]);  x_max = np.max(ax2[1])+1
   y_min = np.min(ax1[1]);  y_max = np.max(ax1[1])+1
   
   reduced_domain = domain[z_min:z_max, x_min:x_max, y_min:y_max]
   return reduced_domain


if auto_reduce_domain_size:
   limits = [ [], [], [], [], [], [] ]
   for i in range(total_time_range[0], total_time_range[1]+1):
      Domain_Size_Reduction(Matrix_4D[i], Tmelt)
      limits[0].append(z_min); limits[1].append(z_max)
      limits[2].append(x_min); limits[3].append(x_max)
      limits[4].append(y_min); limits[5].append(y_max)

   z_min, z_max, x_min, x_max, y_min, y_max = min(limits[0]), max(limits[1]), min(limits[2]), max(limits[3]), min(limits[4]), max(limits[5])
   z_min = 7
   z_max = 8

else:
    z_min =                           0        #3
    z_max =                          9        #9

    x_min =                           15       # 15 --- for salome_4D
    x_max =                          27       # 27 --- for salome_4D

    y_min =                           12       # 12 --- for salome_4D
    y_max =                          97       # 97 --- for salome_4D


mapa =   'INTER  time='+str(time_factor)+', space='+str(space_factor)+'  Z['+str(z_min)+'-'+str(z_max)+'], X['+str(x_min)+'-'+str(x_max)+'], Y['+str(y_min)+'-'+str(y_max)+'], '+str(Tmelt_Celsius)+'°C'+', N='+str(N)+'/'                  

if not os.path.isdir(PATH+mapa):
   os.mkdir(PATH+mapa)


def f(x, tresh):
    T = np.load(PATH+mapa+subfolder+'/salome_'+str(x)+'.npy')
    P = np.zeros(T.shape)
    P[T >(tresh+273)]=1
    return T, P

def p(t_field, tresh):
	p=np.zeros(t_field.shape)
	p[t_field>(tresh+273)]=1
	plt.imshow(p)

def Matrix_Interpolation(A, it):
   if auto_reduce_domain_size:
      B_tuple = (it[0], (z_max-z_min)*it[1], (x_max-x_min)*it[2], (y_max-y_min)*it[3])
   else:
      B_tuple = (it[0], A.shape[1]*it[1], A.shape[2]*it[2], A.shape[3]*it[3])
   new_dims = []
   for original_length, new_length in zip(A.shape, B_tuple):
      new_dims.append(np.linspace(0, original_length-1, new_length))

   coords = np.meshgrid(*new_dims, indexing='ij')
   B = map_coordinates(A, coords, order=1)
   return B

def TimeScaleCounter(tf,sf,tr):
   tsc = [i for i in range(tr[0]*tf*sf, tr[1]*tf*sf)]
   return tsc

def Make_Y_Partitions(Ymin, Ymax, slices):
    dDY = (Ymax-Ymin)/slices
    a=[Ymin+i*dDY for i in range(slices)]
    b=[Ymin+i*dDY+1 for i in range(1,slices+1)]
    YP={}
    for num,(i,j) in enumerate(zip(a,b)):
        YP['YP'+str(num)] = (int(i), int(j))
    return YP


special_case = True

''' ..................... LONG track ............................................................................................................................................ '''
if long_track and not auto_reduce_domain_size:
    start_total = time.time()
    
    for time_snap in range(len(TR_list)-1):
       time_range = (TR_list[time_snap], TR_list[time_snap+1])
       TSC = TimeScaleCounter(time_factor, space_factor, time_range)
       g = time_range[0]
       counter = 0

       YP = Make_Y_Partitions(y_min, y_max-1, N)
       
       for n in range(time_range[0], time_range[1]):
           yp_mapa = yp+'  ['+str(YP[yp][0])+','+str(YP[yp][1])+']'
           if not os.path.isdir(PATH+mapa+yp_mapa):
              os.mkdir(PATH+mapa+yp_mapa+'/')
           print(35*'~')
           print('Processing time step ',n,'/',Time_Range[1]-1,' ..')
           start = time.time()
           if interpolation:
              #inp_mat = Matrix_4D[n:n+2, z_min:z_max, x_min:x_max, y_min:y_max]                      # different sizes of y ranges ---> y_max - y_min set by Domain_Size_Reduction()for given time_ranges
              inp_mat = Matrix_4D[n:n+2, z_min:z_max, x_min:x_max,YP[yp][0]:YP[yp][1]]                  # same sizes of y ranges, i.e. y_max-y_min,  as defined in YP dictionary 
              out_mat = Matrix_Interpolation(inp_mat,increase_tuple)
              if counter==0:
                 od=None
              else:
                 od=1
              del inp_mat
              
           else:
              out_mat = Matrix_4D[ : , z_min:z_max, x_min:x_max, y_min:y_max] 
              od=None

           for i in out_mat[1:]:
               try:
                   subfolder = '/TR'+str(n)+'  ['+str(n)+','+str(n+1)+']'
                   np.save(PATH+mapa+yp_mapa+subfolder+'/salome_'+str(TSC[counter]-g)+'.npy', i)
               except FileNotFoundError:
                   os.mkdir(PATH+mapa+yp_mapa+subfolder)
                   np.save(PATH+mapa+yp_mapa+subfolder+'/salome_'+str(TSC[counter]-g)+'.npy', i)
               counter+=1

           if special_case:      # for connecting indivudal TR folders, i.e. last .npy file from working directory (say TR0) puts into next (TR1), which is created if not yet..
               try: 
                   subf = '/TR'+str(n+1)+'  ['+str(n+1)+','+str(n+2)+']'
                   np.save(PATH+mapa+yp_mapa+subf+'/salome_'+str(TSC[counter-1]-g)+'.npy', out_mat[-1])
               except FileNotFoundError:
                   os.mkdir(PATH+mapa+yp_mapa+subf)
                   np.save(PATH+mapa+yp_mapa+subf+'/salome_'+str(TSC[counter-1]-g)+'.npy', out_mat[-1])

           end = time.time()
           print('Computing time: ',round(end-start, 3),'  sec.')
           print(35*'~'); print()
        
    print('Total computing time =  ',round(time.time()-start_total, 3),'  seconds.')




''' ................. SHORT track ................................................................................................................................................'''
if not long_track and not auto_reduce_domain_size:
    start_total = time.time()
    TSC = TimeScaleCounter(time_factor, space_factor, Time_Range)
    counter = 0
    for n in range(Time_Range[0], Time_Range[1]):
        print(35*'~')
        print('Processing time step ',n,'/',Time_Range[1]-1,' ..')
        start = time.time()
        inp_mat = Matrix_4D[n:n+2, z_min:z_max, x_min:x_max, y_min:y_max]    
        out_mat = Matrix_Interpolation(inp_mat,increase_tuple); del inp_mat

        if counter==0:
            od=None
        else:
            od=1
               
        for i in out_mat[od:]:
            np.save(PATH+mapa+'/salome_'+str(TSC[counter])+'.npy', i)
            counter+=1
        end = time.time()
        print('Computing time: ',round(end-start, 3),'  sec.')
        print(35*'~'); print()
        
    print('Total computing time =  ',round(time.time()-start_total, 3),'  seconds.')



''' ............................ checking the results ............................... '''
'''
maxi = np.where(np.max(mat_4d)== mat_4d)
u=[4*i for i in range(time_range[0], time_range[1]+1)]
plt.plot(u,mat_4d[:, maxi[1][0], maxi[2][0], 25,], marker='s', color='green', markersize=12)


#plt.figure()
v=[]; w=[]
for i in range(21):
   w.append(i)
   q=np.load(PATH+mapa+'/salome_'+str(i)+'.npy')
   v.append(q[maxi[1][0], maxi[2][0], 25,])


plt.plot(w,v, marker='o', color='red')
'''








   
