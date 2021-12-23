""" Creating a .AVI video of CA grain growth, melt-pool and temperature field -----------------------------"""

import numpy as np
import cv2

case      =   'SLM_2D_Source'     ; subcase =   '0002/'

#mapa       =        'INTER  time=40, space=8  Z[7-8], X[15-27], Y[12-97], 1500°C, N=12/'
mapa         =        'INTER  time=1, space=8  Z[0-9], X[15-27], Y[12-97], 1500°C, N=12/'


track =                 'Track_12_layers/'           # for zipping False only  (combining cut figures of individual track..)


PATH   =    'C:/sm-2018-w64-0-3/WORK/'+case+'_Files/post processing database/'+subcase+mapa+track

mapa_animacija =         'animacija_rgb/'   # 'animacija_rgb_bessel/'
koren_slike =                'RGB_final, Z= '   # 'RGB_final bessel, Z= '

ime_videa =                   '3D 1st Moore .avi'



fps = 1                                                    # slik na sekundo

start_frame = 3                                     # zaporedna številka začetnega okvirja
stop_frame = 0                                  # zaporedna številka zadnjega okvirja
frame_step = 1

image_folder = 'C:/Users/akocjan/Desktop/'+mapa_animacija #PATH+mapa_animacija
video_name = image_folder + ime_videa

frame = cv2.imread(image_folder+koren_slike+str(start_frame)+'.png')
height, width, layers = frame.shape

#fourcc = cv2.VideoWriter_fourcc(*'XVID')
#fourcc = cv2.VideoWriter_fourcc(*'MP4V')
fourcc = cv2.VideoWriter_fourcc('M','J','P','G')

video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))


for snap in range(-start_frame, -stop_frame+1, frame_step):
    im=cv2.imread(image_folder+koren_slike+str(-snap)+'.png')
    #print(i, '     ', im.shape)
    '''.................................................................................................... these lines are used to shape all figures to same size .................................................'''
    #im=im[:,:670,:3]
    #cv2.imwrite(image_folder+'figure_'+str(i)+'.png', im)
    #print(i, '     ', im.shape)
    '''....................................................................................................................................................................................................................................'''
    video.write(im)

cv2.destroyAllWindows()
video.release()


""" ........................................................ Reading the. AVI video ............................................."""
"""
import matplotlib.pyplot as plt
plt.ion()

# read video
video_name = PATH+mapa_video+ime_videa
cap = cv2.VideoCapture(video_name)

while(cap.isOpened()):
    ret, frame = cap.read()
    #gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    #plt.figure()
    cv2.imshow('frame', frame)

    if cv2.waitKey(1) & 0XFF == ord('q'):
        break

cap.release()
cap.destroyAllWindows()
"""
