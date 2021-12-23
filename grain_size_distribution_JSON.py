import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.ion()

case      =   'SLM_2D_Source'     ; subcase =   '0002'
PATH   =    'C:/sm-2018-w64-0-3/WORK/'+case+'_Files/post processing database/'+subcase+'/'

mapa       =     'INTER  time=40, space=8  Z[7-8], X[15-27], Y[12-97], 1500°C, N=12/'

JSON_folder =           'Tracks (0,1,2) zipped/'
JSON_root =               'nuclei_data_tracks_(0,1,2)'

pixel_size =     0.625                      # unit: micrometer
pixel_area =    pixel_size**2          # unit: square micrometer


with open(PATH+mapa+'Track_0/nuclei_data.json')as grains:
    g = json.load(grains)
df0 = pd.DataFrame.from_dict(g, orient='index') 

with open(PATH+mapa+'Track_1/nuclei_data_tracks_(0,1).json')as grains:
    g = json.load(grains)
df1 = pd.DataFrame.from_dict(g, orient='index') 

with open(PATH+mapa+JSON_folder+JSON_root+'.json')as grains:
    g = json.load(grains)
df2 = pd.DataFrame.from_dict(g, orient='index')   

G0 = df0['grain size']
G1 = df1['grain size']
G2 = df2['grain size']

G = np.concatenate((G0,G1,G2))

#G_pix = df['grain size']                   # grain size - number of pixels
#G_area = G_pix * pixel_area        # grain size - area

#a,b = np.histogram(G_area)
'''
neg_xrang=1
xrang_poz=10
bin_sajz = 0.1
x_manual = [i*bin_sajz for i in range(int(neg_xrang/bin_sajz), int(xrang_poz/bin_sajz))]
plt.hist(G_area, bins=x_manual )
'''
#x=[i/10 for i in range(0, 535*10, 10)]
#plt.plot(x, G_area, marker = 'o')


'''
x_manual = [i*bin_sajz for i in range(int(neg_xrang/bin_sajz), int(xrang_poz/bin_sajz))]
plt.figure(figsize=(12,8))
a,b = np.histogram(podatki, x_manual)
max_bin_index = np.where(a == a.max())
highest = round(b[max_bin_index[0][0]], 2)
plt.hist(podatki, bins=x_manual,edgecolor='white', label='Povprečje:  '+str(highest)+' '+enota)
plt.legend(loc='upper right', fontsize=20)
plt.title(naslov_grafa+'\n'+file)
plt.xlabel(ime_osi_X, fontsize=18); plt.ylabel(ime_osi_Y, fontsize=18)
plt.xticks(fontsize=14); plt.yticks(fontsize=14)
if save_graph:
    new_file = ime_grafa+'_'+file[:-4]+'.png'
    try:
        plt.gcf().savefig(direktorij+'/Histogrami'+'/'+new_file)
    except FileNotFoundError:
        os.mkdir(direktorij+'/Histogrami')
        plt.gcf().savefig(direktorij+'/Histogrami'+'/'+new_file)
'''
            
