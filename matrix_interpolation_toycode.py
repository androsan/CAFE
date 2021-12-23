import os
import numpy as np
from scipy.ndimage.interpolation import map_coordinates
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random, math, time
plt.ion()

# FE analysis result from SALOME
#PATH = 'C:/sm-2018-w64-0-3/WORK/kode za Salome/Cellular Automata/3D_ss800_f5R8c4f8b23d2N30_10W_t5/'
PATH = 'C:/sm-2018-w64-0-3/WORK/kode za Salome/Cellular Automata/matrix_interpolation_toycode/'
#file =      '3D_ss800_f5R8c4f8b23d2N30_10W_t5.npy'                                     #  .npy filename with FE analysis results (after meshio_postprocessing.py)
mapa =   'CA_1_1'       #  Folder name
#Matrix_4D  =  np.load(PATH+file)


Matrix_4D = np.zeros((6,3,3,3)); m_values = [1 ,6, 11, 16, 21, 26]
for i in range(6):
   Matrix_4D[i] = m_values[i]

                    
''' ****************************************************************************************************************************************'''

FEM_time_factor = 1
scale_up = 2
increase_tuple  =                  (FEM_time_factor * scale_up, scale_up, scale_up, scale_up)                  # resolution increase factor by linear interpolation 
''' ****************************************************************************************************************************************'''

time_range = (1, Matrix_4D.shape[0]-1)
#time_range = (0,31)

auto_reduce_domain_size = False
Tmelt_Celsius = 1000 - 273                                                           #   Melting point; unit:  degrees C [deg. C] - for domain size reduction
Tmelt= Tmelt_Celsius + 273                                                          #   Melting point; unit:  KELVIN [K]

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
   for i in time_range:
      Domain_Size_Reduction(Matrix_4D[i], Tmelt)
      limits[0].append(z_min); limits[1].append(z_max)
      limits[2].append(x_min); limits[3].append(x_max)
      limits[4].append(y_min); limits[5].append(y_max)

   z_min, z_max, x_min, x_max, y_min, y_max = min(limits[0]), max(limits[1]), min(limits[2]), max(limits[3]), min(limits[4]), max(limits[5])

else:
    z_min =                           None    #3
    z_max =                          None    #9

    x_min =                           None    #13
    x_max =                          None    #27

    y_min =                           None    #10
    y_max =                          None    #50

mat_4d = Matrix_4D[:, z_min:z_max, x_min:x_max, y_min:y_max]              #  Reduced 4D matrix

def f(x, tresh):
    T = np.load(PATH+mapa+'/salome_'+str(x)+'.npy')
    P = np.zeros(T.shape)
    P[T >(tresh+273)]=1
    return T, P

def Matrix_Interpolation(A, it):
   if auto_reduce_domain_size:
      B_tuple = (it[0], (z_max-z_min)*it[1], (x_max-x_min)*it[2], (y_max-y_min)*it[3])
   else:
      B_tuple = (A.shape[0]*it[0], A.shape[1]*it[1], A.shape[2]*it[2], A.shape[3]*it[3])
   new_dims = []
   for original_length, new_length in zip(A.shape, B_tuple):
      new_dims.append(np.linspace(0, original_length-1, new_length))

   coords = np.meshgrid(*new_dims, indexing='ij')
   B = map_coordinates(A, coords, order=1)
   return B


'''....................... toycode .......................................................................................... '''

a=np.array([0,1]).astype('float')
#a=np.arange(2).astype('float')
IT = (100,)

def toy_matrix_interpolation(A, it):
   global B_tuple, new_dims, coords, original_length, new_length
   B_tuple = (A.shape[0]*it[0], )
   new_dims = []
   for original_length, new_length in zip(A.shape, B_tuple):
      new_dims.append(np.linspace(0, original_length-1, int(new_length)))
   coords = np.meshgrid(*new_dims, indexing='ij')
   B = map_coordinates(A, coords, order=1)
   return B

b=toy_matrix_interpolation(a, IT); print(np.around(b,2))

plt.plot(a,marker='s', color='grey'); plt.title('original')

plt.figure();plt.title('Interpolated')
plt.plot(b,marker='s', color='magenta')
   
''' ................................................................................................................................... '''


__name__='comment this line, if you want to execute matrix_interpolation'
if __name__=='__main__':

    counter = time_range[0]*FEM_time_factor-1
    start_total = time.time()


    for n in range(time_range[0], time_range[1] + 1):
        print(35*'~')
        print('Processing time step ',n,'/',time_range[1]+1,' ..')
        start = time.time()
        inp_mat = mat_4d[n:n+3]
        out_mat = Matrix_Interpolation(inp_mat,increase_tuple)
        del inp_mat
        for i in out_mat:
            counter+=1
            try:
                np.save(PATH+mapa+'/'+'salome_'+str(counter)+'.npy', i)
            except FileNotFoundError:
                os.mkdir(PATH+mapa)
                np.save(PATH+mapa+'/'+'salome_'+str(counter)+'.npy', i)
        end = time.time()
        print('Computing time: ',round(end-start, 3),'  sec.')
        print(35*'~'); print()
        
    print('Total computing time =  ',round(time.time()-start_total, 3),'  seconds.')


''' ............................ checking the results ............................... '''

'''
u=[i for i in range(time_range[0], time_range[1]+1)]

maxi = np.where(np.max(mat_4d)== mat_4d)

plt.plot(u,mat_4d[time_range[0]:time_range[1]+1, maxi[1][0], maxi[2][0], maxi[3][0],], marker='s', color='green',)


plt.figure()
v=[]; w=[]
for i in range(0, 11):
   w.append(i)
   q=np.load(PATH+mapa+'/salome_'+str(i)+'.npy')
   v.append(q[maxi[1][0], maxi[2][0], maxi[3][0],])

plt.plot(w,v, marker='o', color='red')
'''








   
