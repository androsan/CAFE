''' Program to help set the right parameters of nucleation and growth of grains during SLM '''

import numpy as np
import random, math
import matplotlib.pyplot as plt
plt.ion()

""" =====================================  F  U  N  C  T  I  O  N  S ================================================================"""

def nukleacija_povrsina(pp):                                                  # heterogeneous nucleation (at the liquid-solid interface ---> melt-pool border)
    ''' parameters of heterogeneous nucleation'''
    ns_max =                     5e10                                                #  Number of new nuclei at optimal dT; unit: per SQUARE METER  [m^-2]
    dTs_max0 =                       2                                                 #   Optimal undercooling dT, where the formation of nuclei is the highest; unit:  KELVIN  [K]
    dTs_sigma =                    0.5                                                 #   Standard deviation of Gaussian nuclei distribution [K]
    Ns_Nas =  ns_max/(math.sqrt(2*math.pi))*math.e**(-(pp-dTs_max0)**2/(2*dTs_sigma**2))
    return Ns_Nas

def nukleacija_volumen(vv):                                                   # homogeneous nucleation
    ''' parameters of homogeneous nucleation '''
    nv_max =                      5e14                                             
    dTv_max0 =                       2
    dTv_sigma =                    0.5 
    Nv_Nav =  nv_max/(math.sqrt(2*math.pi))*math.e**(-(vv-dTv_max0)**2/(2*dTv_sigma**2))
    return Nv_Nav

def nakljucje(msm):
   rand = np.random.random_sample(msm.shape)
   #rand = np.random.randint(0,100,(msm.shape))
   #rand[taula==-1]=-1  
   return rand

def taljenje(u, Tm):     #  'u' is temperature field,  'Tm' is melting point (K)
   global taula
   plt.clf()
   taula=np.zeros((Z,X,Y))
   taula[u<Tm] = -1
   return taula

def liquidus(u, liquidus_temp):    #  'u' is temperature field, 'liquidus_temp is'  'Tm + dTliquidus', i.e.
   global likvid
   likvid=np.zeros((Z,X,Y))
   likvid[u>liquidus_temp]=-1
   return likvid

# Loading temperature field (unit: KELVIN)in corresponding space (YP) and time (TR) related folder
def Load_NPY(yp, tr, x):
    npy = np.load(PATH+mapa+yp+'/'+tr+time_factor_folder+'/salome_'+str(x)+'.npy')
    return npy

def Stack_2 (yp1st, yp2nd, tm, cif):
    yp_1st =  Load_NPY(yp1st, tm, cif)
    yp_2nd = Load_NPY(yp2nd, tm, cif)
    yp_stacked = np.dstack((yp_1st,yp_2nd ))
    return yp_stacked

def growth_speed(dT_field):
    vg = 2.03e-4 * dT_field**2 - 0.544e-4 * dT_field
    return vg

""" ====================================================================================================="""
""" ====================================  F  I  L  E  S    &   F  O  L  D  E  R  S ========================================================"""

case =          'SLM_2D_Source'   ;  subcase = '0002'
PATH =        'C:/sm-2018-w64-0-3/WORK/'+case+'_Files/post processing database/'+subcase+'/'

#mapa    =   'INTER  time=40, space=8  Z[7-8], X[15-27], Y[12-97], 1500°C, N=12/'
mapa      =   'INTER  time=1, space=8  Z[0-9], X[15-27], Y[12-97], 1500°C, N=12/'

time_factor_folder = '/time_factor_24/time_factor_3'
""" ====================================================================================================="""

YP0 =                   'YP0  [12,20]'
YP1 =                   'YP1  [19,27]'

TM =                       'TR1  [1,2]'

'''******************'''
TIME_step = 10
'''******************'''

X =       96
Y =       128
z_min =          71
z_max =         72

Z = z_max - z_min

faza = np.zeros((Z,X,Y))                               
NP = np.vectorize(nukleacija_povrsina)
NV = np.vectorize(nukleacija_volumen)
VG = np.vectorize(growth_speed)                      # matrix of growth velocity (vg)

T = Stack_2(YP0, YP1, TM, TIME_step)[z_min:z_max]
time_shift = 1
T_next = Stack_2(YP0, YP1, TM, TIME_step+time_shift)[z_min:z_max]


Tmelt_Celsius = 1507                          ;Tmelt = Tmelt_Celsius + 273   # melting point (Tmelt is in Kelvin)
dTt = T - Tmelt                                     # undecooling [K]

''' Absolute liquidus temperature '''
dTliquidus =    50   

taljenje(T, Tmelt)                              # condition for melting  ----> taula matrix of cell values;  if value -1 then solid (or powder)elif value 0 then liquid
liquidus(T,Tmelt+dTliquidus)              # absolute liquidus line

interface = NP(dTt)                                              
bulk =        NV(dTt)
live =         nakljucje(taula)
vg =          VG(dTt)

ß = 0


""" ====================================>SHOWING the DATA<========================================================"""

X_CUT =      6    ;'''............... X coordinate of the domain for 1D plotting '''
Y_CUT =      54   ;'''............... X coordinate of the domain for 1D plotting '''

# INTERFACE nuclei concentration PLOT(i.e. surface, heterogeneous nucleation..) vs. undercooling (dTt)in given range (delta_T)

delta_T_min = 0
delta_T_max = 4
resolution = 200
delta_T = [i/resolution for i in range(int(delta_T_min*resolution), int(delta_T_max*resolution))]
interface_data = [NP(i) for i in delta_T]
plt.plot(delta_T, interface_data, marker='.', color='magenta'); plt.xlabel(u"\u0394T ("+u"\u00B0C)", fontsize=16); plt.ylabel('Interface Nuclei Concentration', fontsize=16)


# VG PLOT(i.e. vg)vs. undercooling (dTt)in given range (delta_T)
'''
plt.figure()
delta_T_min = -100
delta_T_max = 0
resolution = 200
delta_T = [i/resolution for i in range(int(delta_T_min*resolution), int(delta_T_max*resolution))]
vg_data = [VG(i) for i in delta_T]
plt.plot(delta_T, vg_data, marker='.', color='#8000FF'); plt.xlabel(u"\u0394T ("+u"\u00B0C)", fontsize=16); plt.ylabel('Growth Velocity (m/s)', fontsize=16)
'''


show_TEMPERATURE =            True
show_MELT_POOL =                False
show_UNDERCOOLING =           False
# Nucleation
show_INTERFACE =                   False
show_BULK =                             False
show_LIVE =                              False
# Growth
show_GROWTH_VELOCITY =    False
show_LIQUIDUS =                       False

'" Showing the TEMPERATURE field '''
if show_TEMPERATURE:
    plt.figure()
    plt.imshow(T[0],cmap='hot')
    '''
    plt.plot(T[0,X_CUT,:], marker='.', color='green', label='time step: '+str(TIME_step))
    plt.plot(T_next[0,X_CUT,:], marker='.', color='lime', label='time step: '+str(TIME_step+time_shift))
    plt.title('TEMPERATURE,     TIME step: '+str(TIME_step))
    plt.xlabel('Y-axis'); plt.ylabel('Temperature (K)')
    plt.legend()
    '''

    '" Showing the MELT POOL '''
if show_MELT_POOL:
    plt.figure()
    #plt.imshow(taula[0], cmap='hot')
    plt.plot(taula[0,X_CUT,:], marker='.', color='red')
    plt.title('MELT POOL,     TIME step: '+str(TIME_step))

'" Showing the UNDERCOOLING field '''
if show_UNDERCOOLING:
    plt.figure()
    #dTt[dTt >= 2] = 0 # threshold
    #plt.imshow(dTt[0], cmap='jet')
    plt.plot(dTt[0,X_CUT,:], marker='.', color='blue')
    plt.title('UNDERCOOLING,     TIME step: '+str(TIME_step))

'" Showing the INTERFACE field '''
if show_INTERFACE:
    plt.figure()
    #plt.imshow(interface[0], cmap='jet')
    plt.plot(interface[0,X_CUT,:], marker='.', color='magenta')
    plt.title('INTERFACE,     TIME step: '+str(TIME_step))

'" Showing the BULK field '''
if show_BULK:
    plt.figure()
    #plt.imshow(bulk[0], cmap='jet')
    plt.plot(bulk[0,X_CUT,:], marker='.', color='orange')
    plt.title('BULK,     TIME step: '+str(TIME_step))

'" Showing the LIVE field '''
if show_LIVE:
    plt.figure()
    #plt.imshow(live[0], cmap='jet')
    plt.plot(live[0,X_CUT,:], marker='.', color='black')
    plt.title('LIVE,     TIME step: '+str(TIME_step))


'" Showing the GROWTH VELOCITY field '''
if show_GROWTH_VELOCITY:
    plt.figure()
    #plt.imshow(vg[0], cmap='jet')
    plt.plot(vg[0,X_CUT,:], marker='.', color='#8000FF')
    plt.title('GROWTH VELOCITY,     TIME step: '+str(TIME_step))

    '" Showing the LIQUIDUS field '''
if show_LIQUIDUS:
    plt.figure()
    #plt.imshow(likvid[0], cmap='hot')
    plt.plot(likvid[0,X_CUT,:], marker='.', color='red')
    plt.title('LIQUIDUS,     TIME step: '+str(TIME_step))



PT = (0, X_CUT, Y_CUT)  # Nucleation Point (NP)

print(70*'*')
print(15*' ','Nucleation point data at TIME step ',TIME_step);print()
print('   z = ',PT[0]);  print('   x = ',PT[1]);  print('   y = ',PT[2]); print('   ß = ',ß); print('   M.P.: ',Tmelt,' K   (',Tmelt_Celsius,u"\u00B0C)"); print()
temp = round(T[PT[0],PT[1],PT[2],], 1)
temp_next = round(T_next[PT[0],PT[1],PT[2],], 1)
print('temperature:  ', temp,' K   (',temp-273,u"\u00B0C)")
print('temperature next:  ', temp_next,' K   (',temp_next-273,u"\u00B0C)")
melt_pool_int = taula[PT[0],PT[1],PT[2],]
melt_pool_bool = True if melt_pool_int == 0 else False
print('melt pool :  ',melt_pool_bool,'    (taula: ',melt_pool_int,')')
print('undercooling :  ', round(dTt[PT[0],PT[1],PT[2],], 3),' K')
print('interface :  ', round(interface[PT[0],PT[1],PT[2],], 4))
print('live :  ', round(live[PT[0],PT[1],PT[2],], 4))
print()
print('vg :  ', round(1e3*vg[PT[0],PT[1],PT[2],], 4),' mm/s')
print()
print(40*'-')
print('Booleans'); print(); booleans=[]
print('taula: ',melt_pool_bool); booleans.append(melt_pool_bool)
faza_bool = True if faza[PT[0],PT[1],PT[2],]==0 else False ; print('faza: ',faza_bool); booleans.append(faza_bool)
live_inter_bool = True if live[PT[0],PT[1],PT[2],]<interface[PT[0],PT[1],PT[2],] else False; print('live < interface: ',live_inter_bool); booleans.append(live_inter_bool)
inter_beta_bool = True if interface[PT[0],PT[1],PT[2],]>ß else False; print('interface > ß: ',inter_beta_bool); booleans.append(inter_beta_bool)
T_bool = True if T_next[PT[0],PT[1],PT[2],]<T[PT[0],PT[1],PT[2],] else False ; print('T_next < T: ',T_bool); booleans.append(T_bool)
print()
likvid_bool = True if likvid[PT[0],PT[1],PT[2],]==0 else False; print('likvid: ',likvid_bool)
print(78*'-'); print()
print(20*' ','NUCLEUS is formed :   ','YES !' if all(booleans)else 'NO..')
print(70*'*')

#   Nucleation Condition:
#   taula[k][ii][j]==0 and faza[k][ii][j]==0 and ((live[k][ii][j]/1<interface[k][ii][j] and interface[k][ii][j]>ß) or (live[k][ii][j]/1<bulk[k][ii][j] and bulk[k][ii][j]>ß))and T_next[k][ii][j] < T[k][ii][j]:              

#   Growth Condition:       
#   (faza==0)  &  (taula==0)& (likvid==0) &(lij>=Dsr_1st(R) )&(np.pad(faza,((0,0),(0,0), (1,0)), 'constant')[:,:,:-1]==grain), grain, faza)
           

