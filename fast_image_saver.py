"""++++++++++++++++++++++++++++ FAST creation of .png animation figures by multiprocessing +++++++++++++++++++++++++++++++++++++++++++++++"""
import os
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import matplotlib
matplotlib.use('Agg')

from mpl_toolkits.axes_grid1 import make_axes_locatable

import multiprocessing
import concurrent.futures

first =        10
last =         23
picker =     1
g= 9
smer = 2

PATH = 'C:/sm-2018-w64-0-3/WORK/cooling_plate_Files/' 

#mapa =              PATH+'CA_1'                                       #  Folder name
mapa =                          'CA_4'
flashies =              'flashies'
#flashies =            'parallel_flashies'

animation_figures = 'animation_figures_MP'

cmap = 'jet'
Tmelt = 1200                                                                   #  Melting Point in degrees Celsius
cell_size = 50 / 4                                                            #  FE cell size; unit: MICROMETER [um]
time_step = 5e-05 / 8                                                      #  FE time step; unit: SECOND [s]



# Get time matrix
def Load_Time_Step(num):
    global micro, Z,Y,X, field
    """ Microstructure """
    flash = flashies+'/flashy_snap_'+str(num)+'.npy'
    micro=np.load(mapa+'/'+flash)
    """ Temperature """
    polje = 'salome_'+str(num)+'.npy'
    field_Kelvin=np.load(mapa+'/'+polje)
    #c=np.load(mapa+'/domain_constraints.npy')
    c=np.array([  0,  10,  12,  52,  17, 155])
    field_Kelvin = field_Kelvin[c[0]:c[1], c[2]:c[3], c[4]:c[5]]                               # Constrained temperature field in Kelvins
    field = field_Kelvin - 273                                                                            # Constrained temperature field in degrees Celsius
    Z=micro.shape[0] ; Y=micro.shape[1];  X=micro.shape[2]                               # Domain (microstructure, phase and

def fast_image_saver(im_counter, depth, axis):
    Load_Time_Step(im_counter)
    fig, axs = plt.subplots(nrows=3, sharex=True, sharey=False, figsize=(5,8))
    axs[0].set_title('microstructure (grain orientation-colored)', fontsize=11)
    phase=np.full(micro.shape, (255, 255, 153))
    phase[field<=Tmelt]=(92, 92, 138)
    e=np.zeros((phase.shape[0],phase.shape[1],phase.shape[2],))
    p=phase[:,:,:,0]
    e[p==92]=0; e[p==255]=1
    
    if axis == 2:  # Z_cut
        plt.ylabel('X'); plt.xlabel('Y')
        fig.suptitle('Projekcija v smeri Z, {0}% od začetka domene, TIME: {1} msec. (step # {2})'.format(round((depth+1)*100/Z, 1), "%6.3f" % (round(im_counter*time_step*1000, 3)), im_counter), fontsize=14)
        axs[0].imshow(micro[depth,:,:], cmap=cmap)
        axs[1].imshow(phase[depth,:,:], cmap=cmap)
        x_dim=np.max(np.sum(e[depth,:,:],axis=0))*cell_size ;  y_dim=np.max(np.sum(e[depth,:,:],axis=1) )*cell_size
        temp=field[depth,:,:] 
   
    try:
        plt.savefig(mapa+'/'+animation_figures+'/figure_'+str(im_counter)+'.png', bbox_inches='tight')
    except FileNotFoundError:
        os.mkdir(mapa+'/'+animation_figures)
        plt.savefig(mapa+'/'+animation_figures+'/figure_'+str(im_counter)+'.png', bbox_inches='tight')

    return im_counter


if __name__ == '__main__':

    with concurrent.futures.ThreadPoolExecutor() as executor:
        rezultati = []
        for slika in range(first,last,picker):
            rezultati.append(executor.submit(fast_image_saver, slika, g, smer))
            #fast_image_saver(slika, g, smer)

        out=[i.result()for i in rezultati]
            
        
"""-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""    
