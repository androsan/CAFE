import sys
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
plt.ion()

case      =   'SLM_2D_Source'     ; subcase =   '0002'

#mapa       =        'INTER  time=40, space=8  Z[7-8], X[15-27], Y[12-97], 1500°C, N=12/'
mapa         =        'INTER  time=1, space=8  Z[0-9], X[15-27], Y[12-97], 1500°C, N=12/'


track =                         '2D 1st order Moore, real field/'               # for zipping False only  (combining cut figures of individual track..)
#track =                         '3D 2nd order Moore (2 segments)/'

'''**********************'''
zipping  =    False
'''**********************'''

RGB_mode =     True
bessel =              True

perspektiva =               'Z'        ;  plast_range =             (0,8)
cut_range =                (0,2)

cutoff_limit = 64

cell = 0.625                  # pixel (cell) size in micrometers (for marker in microstructure figure)
marker_size = 10         # marker size in micrometers (to calculate number of pixels in the marker)

X = 96

cuts_RGB =                track+'cuts_RGB/'                               #  folder
cuts_faza =                 track+'cuts_faza/'

cut_file_RGB =          'cut_RGB_'                                        #  fileroot
cut_file_faza =           'cut_faza_'

PATH   =    'C:/sm-2018-w64-0-3/WORK/'+case+'_Files/post processing database/'+subcase+'/'

# Cuts pictures to be without edges. How cool is that! :)
def Cut_Pictures(folder, fileroot, perspective, layer):
    global names
    names =    [fileroot+'{}'.format(i) for i in range (cut_range[0], cut_range[1])]
    if bessel:
        inter_mode = 'bessel' ; bes_txt = 'bessel'
    else:
        inter_mode = None ;  bes_txt = ''
    for nam in names:
        npy = np.load(PATH+mapa+folder+nam+'.npy'); print(nam+'.npy.shape: ', npy.shape)
        if perspective == 'Z':
            if RGB_mode:
                plt.imshow(npy[layer,:,:cutoff_limit,:], interpolation=inter_mode)
            else:
                plt.imshow(npy[layer,:,:cutoff_limit], interpolation=inter_mode)
        elif perspective == 'X':
            plt.imshow(npy[:,layer,:,:], origin='lower', interpolation=inter_mode)
        elif perspective == 'Y':
            plt.imshow(npy[:,:,layer,:], origin='lower', interpolation=inter_mode)
        plt.gca().set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        if names.index(nam)== cut_range[0]:
            plast_mikroni =  str("{:>8.3f}".format((plast-(plast_range[1]-1))* 0.625))+u" \u00B5m"
            plt.text(1,5, plast_mikroni, color='white', fontsize=22)
        elif names.index(nam)== cut_range[1]-1:
            l_edge = 50-int(marker_size/cell)
            plt.hlines(80, l_edge,50, color='white', label='5', linewidth=5)
            plt.text(l_edge, 88,str(marker_size)+' '+u"\u00B5m", color='white', fontsize=22)
        plt.savefig(PATH+mapa+folder+nam+' '+perspective+'='+str(layer)+' '+bes_txt+'.png', bbox_inches='tight', pad_inches = 0)
        plt.clf()
        print()

#X_cutoff_percent =       25       ;  X_cutoff = int(X_cutoff_percent * X / 100)      # Analogy with Hatch Spacing

tracks_RGB = {'Track_0' : {'dir': 'Track_0/cuts_RGB/', 'X_cutoff': 48 , },
                          'Track_1' : {'dir': 'Track_1/cuts_RGB/', 'X_cutoff': 24 , },
                          'Track_2' : {'dir': 'Track_2/cuts_RGB/', 'X_cutoff': 0 , },
              }

tracks_faza = {'Track_0' : {'dir': 'Track_0/cuts_faza/', 'X_cutoff': 48 , },
                          'Track_1' : {'dir': 'Track_1/cuts_faza/', 'X_cutoff': 24 , },
                          'Track_2' : {'dir': 'Track_2/cuts_faza/', 'X_cutoff': 0 , },
              }


zip_folder =           'Tracks (0,1,2) zipped/'         # folder
zip_file_RGB =     'Zip_RGB_'                        # fileroot
zip_file_faza =      'Zip_faza_'                        


zipping_tracks = ['Track_0', 'Track_1', 'Track_2', ]


def Zip_and_Cut(folder, fileroot, perspective, layer):
    names =    [fileroot+'{}'.format(i) for i in range (cut_range[0], cut_range[1])]
    if bessel:
        inter_mode = 'bessel' ; bes_txt = 'bessel'
    else:
        inter_mode = None ;  bes_txt = ''
    for nam in names:
        VC =[]
        for tra in zipping_tracks:
            if RGB_mode:
                npy = np.load(PATH+mapa+tracks_RGB[tra]['dir']+nam+'.npy'); print(nam+'.npy.shape: ', npy.shape)
                VC.append([npy, tracks_RGB[tra]['X_cutoff']])

                if perspective == 'Z':
                    if len(VC)==2:
                        out = np.hstack((VC[1][0][:,VC[1][1]:,:cutoff_limit,:], VC[0][0][:,VC[0][1]:,:cutoff_limit,:] ))
                    
                    elif len(VC)==3:
                        out = np.hstack((VC[2][0][:,VC[2][1]:,:cutoff_limit,:], out ))
                        del VC[0]

            else:
                npy = np.load(PATH+mapa+tracks_faza[tra]['dir']+nam+'.npy'); print(nam+'.npy.shape: ', npy.shape)
                VC.append([npy, tracks_faza[tra]['X_cutoff']])

                if perspective == 'Z':
                    if len(VC)==2:
                        out = np.hstack((VC[1][0][:,VC[1][1]:,:cutoff_limit], VC[0][0][:,VC[0][1]:,:cutoff_limit] ))
                    
                    elif len(VC)==3:
                        out = np.hstack((VC[2][0][:,VC[2][1]:,:cutoff_limit], out ))
                        del VC[0]
                
        #plt.imshow(VC[0][0][:,VC[0][1]:,:,:][0]); plt.title('Track_0')
        #plt.figure(); plt.imshow(VC[1][0][:,VC[1][1]:,:,:][0]); plt.title('Track_1')
        #plt.figure()        
        plt.imshow(out[0], interpolation=inter_mode)
        plt.gca().set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        nam = nam.replace('cut', 'Zip')
        plt.savefig(PATH+mapa+folder+nam+' '+perspective+'='+str(layer)+' '+bes_txt+'.png', bbox_inches='tight', pad_inches = 0)
        

def Combine_Pictures_Horizontally(folder, fileroot, final, perspective, layer):
    if bessel:
        bes_txt = 'bessel'
    else:
        bes_txt = ''
    names =    [fileroot+'{0} {1}={2} {3}'.format(i, perspective,layer,bes_txt) for i in range (cut_range[0], cut_range[1])]
    images = [Image.open(PATH+mapa+folder+nam+'.png') for nam in names]
    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    max_height = max(heights)
    if RGB_mode:
        new_im = Image.new('RGB', (total_width, max_height))
    else:
        new_im = Image.new('RGB', (total_width, max_height))
    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset,0))
        x_offset += im.size[0]
    plt.imshow(new_im)
    new_im.save(PATH+mapa+folder+final)


#*********************************************************************************************************************************************    
''' ------------------------------------------------------------------ execution -------------------------------------------------------------------------'''

for plast in range(plast_range[0], plast_range[1]):
    if RGB_mode:
        if bessel:
            result =     'RGB_final bessel, {0}= {1}.png'.format(perspektiva, plast)
        else:
            result =     'RGB_final, {0}= {1}.png'.format(perspektiva, plast)
        if zipping:
            zip_final =   'zipped_final_RGB.png'
            Zip_and_Cut(zip_folder, cut_file_RGB, perspektiva, plast)
            Combine_Pictures_Horizontally(zip_folder, zip_file_RGB, zip_final, perspektiva, plast)
        else:
            Cut_Pictures(cuts_RGB, cut_file_RGB, perspektiva, plast)
            Combine_Pictures_Horizontally(cuts_RGB, cut_file_RGB, result,perspektiva, plast)
    
    else:
        if bessel:
            result =     'faza_final bessel, {0}= {1}.png'.format(perspektiva, plast)
        else:
            result =     'faza_final, {0}= {1}.png'.format(perspektiva, plast)
        if zipping:
            zip_final =   'zipped_final_faza.png'
            Zip_and_Cut(zip_folder, cut_file_faza, perspektiva, plast)
            #Combine_Pictures_Horizontally(zip_folder, zip_file_faza, zip_final, perspektiva, plast)
        else:
            Cut_Pictures(cuts_faza, cut_file_faza, perspektiva, plast)
            Combine_Pictures_Horizontally(cuts_faza, cut_file_faza, result,perspektiva, plast)

    plt.clf()



    
