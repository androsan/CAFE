import numpy as np
import matplotlib.pyplot as plt
plt.ion()

case      =   'SLM_2D_Source'  ;   subcase =   '0002'

#mapa    =      'INTER  time=40, space=8  Z[7-8], X[15-27], Y[12-97], 1500°C, N=12/Track_1/'
#mapa      =     'INTER  time=40, space=8  Z[4-8], X[15-27], Y[16-97], 1500°C, N=80/'

mapa = 'INTER  time=1, space=8  Z[0-9], X[15-27], Y[12-97], 1500°C, N=12/'

PATH   =    'C:/sm-2018-w64-0-3/WORK/'+case+'_Files/post processing database/'+subcase+'/'

flashies = 'flashies_RGB'

time_factor_folder = '/time_factor_10/time_factor_6'

def Load_NPY(yp, tr, x):
    npy = np.load(PATH+mapa+yp+'/'+tr+time_factor_folder+'/salome_'+str(x)+'.npy')
    return npy

def P(t_field, layer, tresh, name, xlab, col, row, column):
    global p, axs
    p=np.zeros(t_field.shape)
    p[t_field>(tresh+273)]=1
    p[:,:,col]=2

    axs[row, column].imshow(p[layer], cmap='jet')
    #axs[row].imshow(p[layer], cmap='jet')
    
    fig.suptitle(name+' :   '+str(t_field.shape[0] )+' x '+str(t_field.shape[1] )+' x '+str(t_field.shape[2] )+'  cells')
    #axs[0].xlabel(xlab, fontsize=16)

def Stack_2 (yp1st, yp2nd, tm, cif, show_figure, layer, col, row, column):
    global tit, xax, axs
    yp_1st =  Load_NPY(yp1st, tm, cif)
    yp_2nd = Load_NPY(yp2nd, tm, cif)
    yp_stacked = np.dstack((yp_1st,yp_2nd ))
    if show_figure:
        tit = 'YP'+str(yp1st[2:4])+'+YP'+str(yp2nd[2:4])+' --- salome_'+str(cif)+'.npy\n'
        xax = tm[:4]+'  real time = '+str(round(cif*dt*1000, 3))+' msec.'
        P(yp_stacked, layer, Tm_C,tit+xax, xax, col, row, column)
    return yp_stacked

def Make_Y_Partitions(Ymin, Ymax, slices):
    dDY = (Ymax-Ymin)/slices
    a=[Ymin+i*dDY for i in range(slices)]
    b=[Ymin+i*dDY+1 for i in range(1,slices+1)]
    YP={}
    for num,(i,j) in enumerate(zip(a,b)):
        YP['YP'+str(num)] = [int(i),int(j)]
    return YP
    

''' ................................. Interpolation data ....................................'''

time_factor = 40
space_factor = 8

dt = 6.25e-06/space_factor              # Real time step (interpolated); unit:  SECONDS [s]
Tm_C = 1507                                 # Melting Point in degrees Celsius
''' ................................................................................................ '''

y_min =     12    #16
y_max =    97
N =            12     # 80

YP = Make_Y_Partitions(y_min, y_max-1, N)
yp_list = [i+'  '+str(YP[i]).replace(' ', '')for i in YP]

tm_list = ['TR{0}  [{0},{1}]'.format(i,i+1) for i in range (0,8)]

#tm0 = 'TR0  [0,1]'
#tm1 = 'TR1  [1,2]'
#tm2 = 'TR2  [2,3]'
#tm3 = 'TR3  [3,4]'
#tm4 = 'TR4  [4,5]'
#tm5 = 'TR5  [5,6]'

#yp0 = 'YP0  [12,27]'
#yp0 = 'YP0  [12,20]'

#yp1 = 'YP1  [26,41]'
#yp1 = 'YP1  [19,27]'

#yp2 = 'YP2  [40,55]'
#yp2 = 'YP2  [26,34]'


#tm_list = [tm0, tm1, tm2, tm3, tm4, tm5]
#yp_list = [yp0, yp1, yp2]


''' 1st PAIR '''
#S0 = Stack_2(yp0, yp1, tm0, 40, True)
#S1 = Stack_2(yp0, yp1, tm0, 79, True)
#S2 = Stack_2(yp0, yp1, tm1, 119, True)
#S3 = Stack_2(yp0, yp1, tm2, 149, True)
#Sa = Stack_2(yp0, yp1, tm3, 183, True)
''' 1st cut --- when melt leaves yp0, i.e. at ace time 183, then in CA yp0 faza is cut off and saved as .npy '''

''' 2nd PAIR '''
#S4 = Stack_2(yp1, yp2, tm3, 183, True)
#S5 = Stack_2(yp1, yp2, tm4, 239, True)


'''**************************'''
START_step =         271                    # Starting time step number
END_step     =        1594                   # Ending time step number
'''**************************'''

plast_T  =                      71                      #  Index of Z-layer to be shown in plt.imshow() ---->  temperature field
plast_CA =                    0                      #  Index of Z-layer to be shown in plt.imshow() ---->  microstructure


Cut_Off_Percent =   50                    #  Percent of Stack_2 domain lenght at which this domain should be cut off  ---> the left part is saved (.npy and/or .png)the right goes on to CA and so on and on and on..

CA_mode = False

animate = False                                 # False for clicking (or pressing Enter) Animate! button on GUI to show individual time step,  True to click it only once and watch the animation
shranjuj_animacijo = False                 # True to save each individual time frame image as .npy in order to make movie later on ..


from tkinter import *
root = Tk(); root.geometry('300x150+1100+600')

korak = IntVar(); korak.set(START_step)

tt = korak.get()
tm_count = int((START_step+1)/(time_factor*space_factor))

yp_count = 0
cut_count = -1

fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=False, figsize=(5,8))

#yps_test=Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], tt, False, plast_T, 0, 0, 0)
#cutoff_limit = int(yps_test.shape[2]*Cut_Off_Percent/100)

cutoff_limit = 64

def ani():
    global tt, tm_count, yp_count, cut_count, fig, axs, ani_run, npy
    try:
        if not CA_mode:
            yps1=Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], tt, True, plast_T, cutoff_limit, 1,0)
            axs[1,0].set_title('salome_'+str(tt)+'.npy')#+' :   '+str(t_field.shape[0] )+' x '+str(t_field.shape[1] )+'  cells')
            axs[1,0].set_xlabel(xax, fontsize=16)

            #yps1=Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], tt, True, plast_T-1, cutoff_limit, 1,1)
            #axs[1,1].set_title('salome_'+str(tt)+'.npy')#+' :   '+str(t_field.shape[0] )+' x '+str(t_field.shape[1] )+'  cells')
            #axs[1,1].set_xlabel(xax, fontsize=16)
        
        '''
        if np.all(p[:,:,:cutoff_limit]==0):
            print(); print(50*'*',' CUT! ',50*'*')
            print('Cut_Off_Percent = ',Cut_Off_Percent,' %    , cutoff limit = ',cutoff_limit)
            cut_text = 'Time step number  '+str(tt)+',  real time: '+str(round(1000*dt*tt, 3))+' msec.'
            print(cut_text)
            print(106*'*'); print()
            with open(PATH+mapa+'/cut_data.txt', 'a')as cuttxt:
                cuttxt.write(cut_text+'\n')

            yp_count+=1
            cut_count+=1
            
            plt.savefig(PATH+mapa+'/cut_'+str(cut_count)+'.png', bbox_inches='tight')
            np.save(PATH+mapa+'/cut_'+str(cut_count)+'.npy', p[:,:,:cutoff_limit])
            
        '''
        if CA_mode:
            npy = np.load(PATH+mapa+flashies+'/flashy_RGB_'+str(tt)+'.npy')
            plt.imshow(npy[plast_CA])

            #axs[0,0].imshow(npy[plast_CA])
            #axs[0,0].set_title('flashy_RGB_'+str(tt)+'.npy')
            #axs[0,0].set_xlabel('Z= '+str(plast_CA))
            
            #axs[0,1].imshow(npy[0])
            #axs[0,1].set_xlabel('Z=46')
             
        
    except FileNotFoundError:
        print(); print('FileNotFound Exception !')
        #print('STEP: ',tt,'  ,  yp_count :  ', yp_count,'     tm_count :  ',tm_count); print()
        
        tm_count+=1
        #yps1=Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], tt, True, plast_T, cutoff_limit)

        yps1=Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], tt, True, plast_T, cutoff_limit, 1,0)
        axs[1,0].set_title('salome_'+str(tt)+'.npy')#+' :   '+str(t_field.shape[0] )+' x '+str(t_field.shape[1] )+'  cells')
        axs[1,0].set_xlabel(xax, fontsize=16)

        yps1=Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], tt, True, plast_T-1, cutoff_limit, 1,1)
        axs[1,1].set_title('salome_'+str(tt)+'.npy')#+' :   '+str(t_field.shape[0] )+' x '+str(t_field.shape[1] )+'  cells')
        axs[1,1].set_xlabel(xax, fontsize=16)

        
    if shranjuj_animacijo:
        try:
            plt.savefig(directory.get()+'/'+animation_figures+'/figure_'+str(tt)+'.png', bbox_inches='tight')
        except FileNotFoundError:
            os.mkdir(directory.get()+'/'+animation_figures)
            plt.savefig(directory.get()+'/'+animation_figures+'/figure_'+str(tt)+'.png', bbox_inches='tight')

    print('STEP: ',tt,'  ,  yp_count :  ', yp_count,'     tm_count :  ',tm_count)            
    tt+=1

    if animate:
        ani_run=root.after(500,ani)


ani_safety_switch = False
def start_animation():
    global ani_safety_switch
    if not ani_safety_switch:
        tt = korak.get()
        ani()
        if animate:
            ani_safety_switch=True
        else:
            ani_safety_switch=False

def stop_animation():
    global ani_safety_switch, tt
    ani_safety_switch=False
    try:
        root.after_cancel(ani_run)
    except NameError:
        pass
    tt= korak.get()
    
    return tt


ent1= Entry(root, textvariable= korak, font=('Arial', '30'), justify='center', width=7); ent1.grid(row=19, column=11, columnspan=2, sticky=W)
but5= Button(root, text= 'Animate!', font=('Arial', '22'), fg='green', width=9, command=start_animation); but5.grid(row=21, column=12, sticky=W)
but6= Button(root, text= 'Stop', font=('Arial', '22'), fg='red', width=9, command=stop_animation); but6.grid(row=22, column=12, sticky=W)

root.bind('<Return>', lambda event=None: but5.invoke())     #  command to invoke but5 (Animate!) with Return (Enter) key on keyboard

def on_enter(e):
    but5['background'] = 'red'

def on_leave(e):
    but5['background'] = 'SystemButtonFace'

but5.bind("<Enter>", on_enter)
but5.bind("<Leave>", on_leave)
