# coding=utf-8
"""  Program za nukleacijo in rast zrn po 'Cellular Automata Mesh Dependency' metodi v 3D,
      avtor:  dr. Andraz Kocjan, Institut za kovinske materiale in tehnologije, Oktober 2020  """
import sys
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import multi_dot
import random, math, time
import psutil
plt.ion()
np.random.seed(90983)

#sys.path.append('C:\\Users\\akocjan\\AppData\\Local\\Programs\\Python\\Python36-32\\Lib\\site-packages\\')

""" =====================================  T  K  I  N  T  E  R  ================================================================"""
from tkinter import *
from tkinter import filedialog, messagebox
root = Tk(); root.geometry('1250x525+200+200'); basic_color = '#666699'
root.title('SLM Microstructure Simulator by Andro'); root.config(bg=basic_color)

'''...textvariables .................................. '''
# OUTPUT
var0=IntVar();  #   Absolute time step index (-) :: i
var1=StringVar();  #   CPU time of current time step (seconds) :: step_cpu
var2=IntVar();  #   Number of active grains (-) :: active_grains
var3=IntVar();  #   Number of ALL grains (-) :: grain_counter
var4=StringVar();  #   Printed length (micrometers) :: real_length
var5=StringVar();  #   Printing time (miliseconds):: real_time
var6=StringVar();  #   CPU % (%) :: CPU 
var7=StringVar();  #   RAM % (%) :: RAM

# INPUT (here, default values of all input variables, i.e. Entries, etc., are defined using var.set() method)
''' FEM inputs '''
var8=DoubleVar();  var8.set(800)#   Laser Scanning Speed (mm/s) :: FEM_scanning_speed
var9=DoubleVar();  var9.set(5e-06)#   FEM Cell Size (m) :: FEM_cell_size
var10=IntVar(); var10.set(8) #   Space Interpolation Factor (-) :: space_factor
var11=IntVar(); var11.set(8) #   Extra Time Interpolation Factor (-) :: extra_time_factor
''' Material inputs '''
var12=IntVar(); var12.set(1507)#   Melting point of metal powder (deg. Celsius):: Tmelt_Celsius
var13=IntVar(); var13.set(50)#   Delta T (liquidus)(deg. Celsius):: dTliquidus

''' CA inputs '''
tracks_database = {
        'iso': ['2D 1st order Moore, iso field/', '2D 2nd order Moore, iso field/', '3D 1st order Moore, iso field/', '3D 2nd order Moore, iso field/'],
        'real': ['2D 1st order Moore, real field/', '2D 2nd order Moore, real field/', '3D 1st order Moore, real field/', '3D 2nd order Moore, real field/'],}
track = 'default_track/'   #track = tracks_database['real'][0]
var20=StringVar(); var20.set(track)#   Name of the Track (-) :: track
var21=IntVar(); var21.set(71)#    Z-index of domain bottom layer :: z_min
var22=IntVar(); var22.set(72)#    Z-index of upper layer :: z_max
var23=IntVar(); var23.set(5)#   Index of Grains and Directions Selection Mechanism (random - ordered):: pick_selection_mechanism
var24=DoubleVar(); var24.set(0.95)#   Frequency of new nuclei formed, BETA ranges from 0 to 1, high to low, respectively :: BETA
var25=IntVar(); var25.set(0)#   Starting time step number, which relates to suffix of salome_{}.npy file :: START_step
var26=IntVar(); var26.set(50)#   End time step number, which relates to suffix of salome_{}.npy file :: END_step
var27=IntVar(); var27.set(1)#   Time Shift = 1, 2 or 3 :: time_shift
var28=IntVar(); var28.set(45)#   Number of possible random alfa, beta, gama choices for grain orientation randomization :: rp
var29=DoubleVar(); var29.set(0)#   Neighbourhood randomness (it's a float between zero and 0.49):: epsilon
var30=IntVar(); var30.set(0)#   Index of YP first half folder :: yp_count
var31=IntVar(); var31.set(1)#   Index of TR folder :: tm_count

zavora = StringVar()
'''...lables ............................................ '''

lab_edge0 = Label(root, text=1*' ', width=2, height=1, bg=basic_color);  lab_edge0.grid(row=99, column=98, columnspan=1)

lab0_col='white'
lab0L = Label(root, anchor=E, text='Absolute time step:', font=('Arial', '14'), bg='black', fg=lab0_col, height=4, width=25); lab0L.grid(row=100, column=99, columnspan=1)
lab0 = Label(root, textvariable=var0, font=('Arial', '29'), bg='black', fg=lab0_col, height=2, width=6); lab0.grid(row=100, column=100, columnspan=1)
lab0R = Label(root, anchor=W, text=' ', font=('Arial', '14'), bg='black', fg=lab0_col, height=4, width=7); lab0R.grid(row=100, column=101, columnspan=1)

lab1_col='yellow'
lab1L = Label(root, anchor=E, text='CPU of current time step:', font=('Arial', '14'), bg='black', fg=lab1_col, height=4, width=25); lab1L.grid(row=100, column=102, columnspan=1)
lab1 = Label(root, textvariable=var1, font=('Arial', '29'), bg='black', fg=lab1_col, height=2, width=6); lab1.grid(row=100, column=103, columnspan=1)
lab1R = Label(root, anchor=W, text='sec.', font=('Arial', '14'), bg='black', fg=lab1_col, height=4, width=7); lab1R.grid(row=100, column=104, columnspan=1)

lab2_col='orange'
lab2L = Label(root, anchor=E, text='Active grains:', font=('Arial', '14'), bg='black', fg=lab2_col, height=4, width=25); lab2L.grid(row=101, column=102, columnspan=1)
lab2 = Label(root, textvariable=var2, font=('Arial', '29'), bg='black', fg=lab2_col, height=2, width=6); lab2.grid(row=101, column=103, columnspan=1)
lab2R = Label(root, anchor=W, text=' ', font=('Arial', '14'), bg='black', fg=lab2_col, height=4, width=7); lab2R.grid(row=101, column=104, columnspan=1)

lab3_col='red'
lab3L = Label(root, anchor=E, text='TOTAL grains:', font=('Arial', '14'), bg='black', fg=lab3_col, height=4, width=25); lab3L.grid(row=101, column=99, columnspan=1)
lab3 = Label(root, textvariable=var3, font=('Arial', '29'), bg='black', fg=lab3_col, height=2, width=6); lab3.grid(row=101, column=100, columnspan=1)
lab3R = Label(root, anchor=W, text=' ', font=('Arial', '14'), bg='black', fg=lab3_col, height=4, width=7); lab3R.grid(row=101, column=101, columnspan=1)

lab4_col='lime'
lab4L = Label(root, anchor=E, text='Printed length:', font=('Arial', '14'), bg='black', fg=lab4_col, height=4, width=25); lab4L.grid(row=102, column=99, columnspan=1)
lab4 = Label(root, textvariable=var4, font=('Arial', '29'), bg='black', fg=lab4_col, height=2, width=6); lab4.grid(row=102, column=100, columnspan=1)
lab4R = Label(root, anchor=W, text='mm', font=('Arial', '14'), bg='black', fg=lab4_col, height=4, width=7); lab4R.grid(row=102, column=101, columnspan=1)

lab5_col='lime'
lab5L = Label(root, anchor=E, text='Printing time:', font=('Arial', '14'), bg='black', fg=lab5_col, height=4, width=25); lab5L.grid(row=102, column=102, columnspan=1)
lab5 = Label(root, textvariable=var5, font=('Arial', '29'), bg='black', fg=lab5_col, height=2, width=6); lab5.grid(row=102, column=103, columnspan=1)
lab5R = Label(root, anchor=W, text='msec.', font=('Arial', '14'), bg='black', fg=lab5_col, height=4, width=7); lab5R.grid(row=102, column=104, columnspan=1)

lab6_col='aqua'
lab6L = Label(root, anchor=E, text='CPU (%):', font=('Arial', '14'), bg='black', fg=lab6_col, height=4, width=25); lab6L.grid(row=103, column=102, columnspan=1)
lab6 = Label(root, textvariable=var6, font=('Arial', '29'), bg='black', fg=lab6_col, height=2, width=6); lab6.grid(row=103, column=103, columnspan=1)
lab6R = Label(root, anchor=W, text='%', font=('Arial', '14'), bg='black', fg=lab6_col, height=4, width=7); lab6R.grid(row=103, column=104, columnspan=1)

lab7_col='aqua'
lab7L = Label(root, anchor=E, text='RAM (%):', font=('Arial', '14'), bg='black', fg=lab7_col, height=4, width=25); lab7L.grid(row=103, column=99, columnspan=1)
lab7 = Label(root, textvariable=var7, font=('Arial', '29'), bg='black', fg=lab7_col, height=2, width=6); lab7.grid(row=103, column=100, columnspan=1)
lab7R = Label(root, anchor=W, text='%', font=('Arial', '14'), bg='black', fg=lab7_col, height=4, width=7); lab7R.grid(row=103, column=101, columnspan=1)

'''...entries ............................................ '''


# ******************************MENUBAR *********************************************************************************************************************************
menubar = Menu(); root.config(menu=menubar)
# Cascades
men0=Menu(); men1=Menu(); men2=Menu(); men3=Menu();  men4=Menu()
menubar.add_cascade(label="FEM data", menu=men0)
menubar.add_cascade(label="Booleans", menu=men1)
menubar.add_cascade(label="Parameters", menu=men2)
menubar.add_cascade(label="Moore", menu=men3)
menubar.add_cascade(label="Help", menu=men4)

# Commands

# 1.1.   FEM data / Find FEM data

case =           'SLM_2D_Source'   ;  subcase = '0002'
PATH =         'C:/sm-2018-w64-0-3/WORK/'+case+'_Files/post processing database/'+subcase+'/'
mapa    =       'INTER  time=1, space=8  Z[0-9], X[15-27], Y[12-97], 1500°C, N=12/'

directory = StringVar(); directory.set(PATH+mapa) # ---> ! when making executable, set directory to something common like C:\User\Desktop\ !
def FEM_data():
    direktorij = filedialog.askdirectory()
    directory.set(direktorij+'/')
men0.add_command(label="Find FEM data", command=FEM_data)

# 2.1.   Booleans
from_beginning=BooleanVar(); from_beginning.set( True)  # True to start from beginning, False to continue from previously saved simulation (KickOff)
save_kickoff=BooleanVar(); save_kickoff.set(True)           # if True it saves kickoff parameters (runs Save_KickOff())
save_last_cut_figure=BooleanVar(); save_last_cut_figure.set(True)
avtomatska_nukleacija=BooleanVar(); avtomatska_nukleacija.set(True)

delete_inactive_grains_ID =         False     

save_flashy_as_RGB=BooleanVar(); save_flashy_as_RGB.set(False)
save_flashy_as_faza=BooleanVar(); save_flashy_as_faza.set(False)

save_cut_as_RGB=BooleanVar();   save_cut_as_RGB.set(True)
save_cut_as_faza=BooleanVar();   save_cut_as_faza.set(True)

men1.add_checkbutton(label="Start From Beginning", onvalue=1, offvalue=0, variable=from_beginning)
men1.add_checkbutton(label="Save KickOff", onvalue=1, offvalue=0, variable=save_kickoff)
men1.add_checkbutton(label="Save Last Cut Figure", onvalue=1, offvalue=0, variable=save_last_cut_figure)
men1.add_checkbutton(label="Auto Nucleation", onvalue=1, offvalue=0, variable=avtomatska_nukleacija)
men1.add_checkbutton(label="Save RGB Images", onvalue=1, offvalue=0, variable=save_flashy_as_RGB)
men1.add_checkbutton(label="Save Grain ID Images", onvalue=1, offvalue=0, variable=save_flashy_as_faza)
men1.add_checkbutton(label="Save RGB Cuts", onvalue=1, offvalue=0, variable=save_cut_as_RGB)
men1.add_checkbutton(label="Save Grain ID Cuts", onvalue=1, offvalue=0, variable=save_cut_as_faza)


# 3.1.   FEM Parameters 

def Safe_Refresh_of_FEM_parameters(subwin, safety_message):
    global FEM_scanning_speed, FEM_cell_size, space_factor, FEM_time_factor, extra_time_factor, Tmelt_Celsius, dTliquidus
    if safety_message:
        smout = messagebox.askquestion("Parameter Change Confirmation", "Do you really want to use these parameters?", icon='warning')
        if smout == 'yes':
            FEM_scanning_speed=var8.get()
            FEM_cell_size=var9.get()
            space_factor=var10.get()
            extra_time_factor=var11.get()
            FEM_time_factor=5
            Tmelt_Celsius=var12.get()
            dTliquidus=var13.get()
        else:
            pass
        subwin.destroy()
    else:
        FEM_scanning_speed=var8.get()
        FEM_cell_size=var9.get()
        space_factor=var10.get()
        extra_time_factor=var11.get()
        FEM_time_factor=5
        Tmelt_Celsius=var12.get()
        dTliquidus=var13.get()

def Enter_FEM_n_Material_parameters():
    sub0 = Toplevel(root)
    sub0.title('FEM Analysis and Material Parameters')
    sub0.geometry('700x700+10+10'); basic_color='lightgreen'
    sub0.config(bg=basic_color)
    label_font=('Arial', '17'); label_bg='black'; label_fg='white'; entry_bg='white'; entry_fg='blue'; entry_font=('Arial', '35')

    lab_edge0 = Label(sub0, text=10*' ', width=5, height=2, bg=basic_color);  lab_edge0.grid(row=5, column=8, columnspan=1)
    lab0 = Label(sub0, text='Laser Scanning Speed (mm/s): ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab0.grid(row=10, column=9, columnspan=1)
    ent0 = Entry(sub0, textvariable=var8, width='8', bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent0.grid(row=10, column=10, columnspan=1)
    lab1 = Label(sub0, text='FEM Cell Size (m): ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab1.grid(row=11, column=9, columnspan=1)
    ent1 = Entry(sub0, textvariable=var9, width='8', bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent1.grid(row=11, column=10, columnspan=1)
    lab2 = Label(sub0, text='Space Interpolation Factor (-): ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab2.grid(row=12, column=9, columnspan=1)
    ent2 = Entry(sub0, textvariable=var10, width='8', bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent2.grid(row=12, column=10, columnspan=1)
    lab3 = Label(sub0, text='Extra Time Interpolation Factor (-): ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab3.grid(row=13, column=9, columnspan=1)
    ent3 = Entry(sub0, textvariable=var11, width='8', bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent3.grid(row=13, column=10, columnspan=1)
    lab4 = Label(sub0, text='Powder Melting Point ('+u"\u00B0C): ", font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab4.grid(row=14, column=9, columnspan=1)
    ent4 = Entry(sub0, textvariable=var12, width='8', bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent4.grid(row=14, column=10, columnspan=1)

    lab5 = Label(sub0, text=u"\u0394T(liquidus)   ("+u"\u00B0C): ", font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab5.grid(row=15, column=9, columnspan=1)
    ent5 = Entry(sub0, textvariable=var13, width='8', bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent5.grid(row=15, column=10, columnspan=1)

    lab_edge1 = Label(sub0, text=10*' ', width=5, height=3, bg=basic_color);  lab_edge1.grid(row=49, column=9, columnspan=1)
    but0 = Button(sub0, text='Change Parameters', command=lambda:Safe_Refresh_of_FEM_parameters(sub0, True), font=('Arial', '20'), fg='black', width=22); but0.grid(row=50, column=9, columnspan=1)
    but1 = Button(sub0, text='Close', command=lambda:sub0.destroy(), font=('Arial', '20'), fg='black', width=11); but1.grid(row=50, column=10, columnspan=1)



# 3.2.   CA Parameters 

def Safe_Refresh_of_CA_parameters(subwin, safety_message):
    global track, z_min, z_max, pick_selection_mechanism, BETA, START_step, END_step, time_shift, rp, epsilon, yp_count, tm_count
    if safety_message:
        smout = messagebox.askquestion("Parameter Change Confirmation", "Do you really want to use these parameters?", icon='warning')
        if smout == 'yes':
            track=var20.get()
            z_min=var21.get()
            z_max=var22.get()
            pick_selection_mechanism=var23.get()
            BETA=var24.get()
            START_step=var25.get()
            END_step=var26.get()
            time_shift=var27.get()
            rp=var28.get()
            epsilon=var29.get()
            yp_count=var30.get()
            tm_count=var31.get()
            
        else:
            pass
        subwin.destroy()
    else:
        track=var20.get()
        z_min=var21.get()
        z_max=var22.get()
        pick_selection_mechanism=var23.get()
        BETA=var24.get()
        START_step=var25.get()
        END_step=var26.get()
        time_shift=var27.get()
        rp=var28.get()
        epsilon=var29.get()
        yp_count=var30.get()
        tm_count=var31.get()
        

def Enter_CA_parameters():
    sub1 = Toplevel(root)
    sub1.title('CA Parameters')
    sub1.geometry('1100x800+50+100'); basic_color = 'lightblue'
    sub1.config(bg=basic_color)
    label_font=('Arial', '14'); label_bg='black'; label_fg='white'; entry_bg='white'; entry_fg='blue'; entry_font=('Arial', '30')

    lab_edge0 = Label(sub1, text=10*' ', width=5, height=1, bg=basic_color);  lab_edge0.grid(row=9, column=8, columnspan=1)

    lab0 = Label(sub1, text='Track Name: ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab0.grid(row=10, column=9, columnspan=1)
    ent0 = Entry(sub1, textvariable=var20, width=40, bg=entry_bg, fg=entry_fg, font=('Arial', '20'), justify='center'); ent0.grid(row=10, column=10, columnspan=1, sticky='SW')
    lab_edge1 = Label(sub1, text=10*' ', width=5, height=1, bg=basic_color);  lab_edge1.grid(row=11, column=9, columnspan=1)
    lab1 = Label(sub1, text='Z-bottom layer index: ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab1.grid(row=12, column=9, columnspan=1)
    ent1 = Entry(sub1, textvariable=var21, width=8, bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent1.grid(row=12, column=10, columnspan=1, sticky='W')
    lab2 = Label(sub1, text='Z-upper layer index: ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab2.grid(row=13, column=9, columnspan=1)
    ent2 = Entry(sub1, textvariable=var22, width=8, bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent2.grid(row=13, column=10, columnspan=1, sticky='W')
    lab3 = Label(sub1, text='Selection Mechanism: ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab3.grid(row=14, column=9, columnspan=1)
    ent3 = Entry(sub1, textvariable=var23, width=8, bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent3.grid(row=14, column=10, columnspan=1, sticky='W')
    lab_edge2 = Label(sub1, text=10*' ', width=5, height=1, bg=basic_color);  lab_edge2.grid(row=15, column=8, columnspan=1)
    lab4 = Label(sub1, text='Nucleation Frequency (BETA): ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab4.grid(row=16, column=9, columnspan=1)
    ent4 = Entry(sub1, textvariable=var24, width=8, bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent4.grid(row=16, column=10, columnspan=1, sticky='W')
    lab5 = Label(sub1, text='START_step: ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab5.grid(row=17, column=9, columnspan=1)
    ent5 = Entry(sub1, textvariable=var25, width=8, bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent5.grid(row=17, column=10, columnspan=1, sticky='W')
    lab6 = Label(sub1, text='END_step: ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab6.grid(row=18, column=9, columnspan=1)
    ent6 = Entry(sub1, textvariable=var26, width=8, bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent6.grid(row=18, column=10, columnspan=1, sticky='W')
    lab10 = Label(sub1, text='Index of YP folder: ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab10.grid(row=19, column=9, columnspan=1)
    ent10 = Entry(sub1, textvariable=var30, width=8, bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent10.grid(row=19, column=10, columnspan=1, sticky='W')
    lab11 = Label(sub1, text='Index of TR folder: ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab11.grid(row=20, column=9, columnspan=1)
    ent11 = Entry(sub1, textvariable=var31, width=8, bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent11.grid(row=20, column=10, columnspan=1, sticky='W')
    lab_edge3 = Label(sub1, text=10*' ', width=5, height=1, bg=basic_color);  lab_edge3.grid(row=21, column=8, columnspan=1) 
    lab7 = Label(sub1, text='Time Shift: ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab7.grid(row=22, column=9, columnspan=1)
    ent7 = Entry(sub1, textvariable=var27, width=8, bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent7.grid(row=22, column=10, columnspan=1, sticky='W')
    lab8 = Label(sub1, text='Possible Grain Orientations: ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab8.grid(row=23, column=9, columnspan=1)
    ent8 = Entry(sub1, textvariable=var28, width=8, bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent8.grid(row=23, column=10, columnspan=1, sticky='W')
    lab9 = Label(sub1, text='Neighbourhood Randomness: ', font=label_font, bg=label_bg, fg=label_fg, height=2, width=30); lab9.grid(row=24, column=9, columnspan=1)
    ent9 = Entry(sub1, textvariable=var29, width=8, bg=entry_bg, fg=entry_fg, font=entry_font, justify='center'); ent9.grid(row=24, column=10, columnspan=1, sticky='W')
    
    lab_edge9 = Label(sub1, text=10*' ', width=5, height=2, bg=basic_color);  lab_edge9.grid(row=49, column=9, columnspan=1)
    but0 = Button(sub1, text='Change Parameters', command=lambda:Safe_Refresh_of_CA_parameters(sub1, True), font=('Arial', '20'), fg='black', width=22); but0.grid(row=50, column=9, columnspan=1)
    but1 = Button(sub1, text='Close', command=lambda:sub1.destroy(), font=('Arial', '20'), fg='black', width=11); but1.grid(row=50, column=10, columnspan=1, sticky='W')

    
men2.add_command(label="FEM & Material", command=Enter_FEM_n_Material_parameters)
men2.add_command(label="CA", command=Enter_CA_parameters)


# 4.1.   Moore

moore_var = StringVar(); moore_var.set('1st')

def Select_Groups(kriterij):
    global smeri_database, smeri
    if kriterij=='1st':
        smeri_database = np.array([np.concatenate((G1,G4))])# 1 segment --- 2D , I. order Moore
        smeri = smeri_database.flatten()
    elif kriterij=='2nd':
        smeri_database = np.array([np.concatenate((G1,G4, G9, G14, G17))])# 1 segment --- 2D , II. order Moore
        smeri = smeri_database.flatten()
    elif kriterij=='3rd':
        smeri_database = np.array([np.concatenate([G1, G4]), np.concatenate([G9, G14, G17]),], dtype=object)# 2 segments --- 2D , II. order Moore
        smeri = np.array([np.concatenate((G1,G4, G9, G14, G17))]).flatten()
    elif kriterij=='4th':
        smeri_database = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,9)]))])  # 1 segment --- 3D , I. order Moore
        smeri = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,9)]))]).flatten()
    elif kriterij=='5th':
        smeri_database = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,30)]))]) # 1 segment --- 3D , II. order Moore
        smeri = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,30)]))]).flatten()
    elif kriterij=='6th':
        smeri_database = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,9)])),
                                                    np.concatenate(tuple([eval('G{}'.format(i)) for i in range(9,30)])),], dtype=object) # 2 segments --- 3D , II. order Moore
        smeri = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,30)]))]).flatten()

men3.add_radiobutton(label="2D, I. order Moore (1 segment)", variable=moore_var, value='1st')
men3.add_radiobutton(label="2D, II. order Moore (1 segment)", variable=moore_var, value='2nd')
men3.add_radiobutton(label="2D, II. order Moore (2 segments)", variable=moore_var, value='3rd')
men3.add_radiobutton(label="3D, I. order Moore (1 segment)", variable=moore_var, value='4th')
men3.add_radiobutton(label="3D, II. order Moore (1 segment)", variable=moore_var, value='5th')
men3.add_radiobutton(label="3D, II. order Moore (2 segments)", variable=moore_var, value='6th')

# 5.1. Multiprocessing
########################################################################################## MPI
run_1st = True
run_2nd = True
run_3rd = True

MPI  =       True
########################################################################################## MPI


''' . . . . . . . . . . . . . . . . . . . . . STRUCTURING (subarrays)groups (G1,..)within smeri_database to form segments . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .'''
#smeri_database = np.array([np.concatenate(tuple([eval('G{}'.format(i)) for i in range(1,30)]))])
#                            .reshape((4,31))                                                                                                              # 4 segments (equal size)--- 3D , II. order Moore


#smeri_database = np.array([np.concatenate((G1,G2,G3)),np.concatenate((G4,G5,G6)),np.concatenate((G7,G8)),  # 9 segments --- 3D , II. order Moore
#        np.concatenate((G9,G10,G11,G12,G13)), np.concatenate((G14,G15,G16)), np.concatenate((G17,G18,G19)),
#        np.concatenate((G20,G21,G22,G23)), np.concatenate((G24,G25,G26,G27)), np.concatenate((G28,G29)), ])

''' . . . . . . . . . . . . . . . . . . . . >>> smeri_database.flatten() to get ::: 'smeri' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .'''
#smeri_database = np.array([['001', '00_1', '010', '0_10',]])

#smeri_database = np.array([np.array(['001']), np.array(['00_1'])])
#smeri = smeri_database.flatten()

# *********************************************************************************************************************************************************************************


""" =====================================  F  U  N  C  T  I  O  N  S ================================================================"""
def CPU_RAM():
    CPU=psutil.cpu_percent(); RAM= psutil.virtual_memory().percent
    print('CPU (%): ', CPU,' %       RAM (%): ', RAM,' %')
    return CPU, RAM

def f(n):
    f=np.load('C:/sm-2018-w64-0-3/WORK/SLM_2D_Source_Files/post processing database/0002/'+
              'INTER  time=1, space=8  Z[0-9], X[15-27], Y[12-97], 1500°C, N=12/2D 1st order Moore, real field/flashies_RGB/flashy_RGB_{}.npy'.format(n))
    return f

def time_counter_original(neg_index):
    global Negatives
    Negatives[neg_index] += dt
    return Negatives

        
# Loading temperature field (unit: KELVIN)in corresponding space (YP) and time (TR) related folder
def Load_NPY(yp, tr, x):
    npy = np.load(directory.get()+yp+'/'+tr+time_factor_folder+'/salome_'+str(x)+'.npy')
    return npy

def Stack_2 (yp1st, yp2nd, tm, cif):
    yp_1st =  Load_NPY(yp1st, tm, cif)
    yp_2nd = Load_NPY(yp2nd, tm, cif)
    yp_stacked = np.dstack((yp_1st,yp_2nd ))
    return yp_stacked

def Make_Y_Partitions(Ymin, Ymax, slices):
    dDY = (Ymax-Ymin)/slices
    a=[Ymin+i*dDY for i in range(slices)]
    b=[Ymin+i*dDY+1 for i in range(1,slices+1)]
    YP={}
    for num,(i,j) in enumerate(zip(a,b)):
        YP['YP'+str(num)] = [int(i),int(j)]
    return YP

# **** Domain Size Reduction ****
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

# **** Nucleation functions ****
def nukleacija_povrsina(pp):                                                   # heterogeneous nucleation (at the liquid-solid interface ---> melt-pool border)
    ''' parameters of heterogeneous nucleation'''
    ns_max =                      5e10                                                 #  Number of new nuclei at optimal dT; unit: per SQUARE METER  [m^-2]
    dTs_max0 =                       2                                                 #   Optimal undercooling dT, where the formation of nuclei is the highest; unit:  KELVIN  [K]
    dTs_sigma =                    0.5                                                 #   Standard deviation of Gaussian nuclei distribution [K]
    Ns_Nas =  ns_max/(math.sqrt(2*math.pi))*math.e**(-(pp-dTs_max0)**2/(2*dTs_sigma**2))
    #Ns_Nas = 0
    return Ns_Nas

def nukleacija_volumen(vv):                                                   # homogeneous nucleation
    ''' parameters of homogeneous nucleation '''
    nv_max =                      5e14                                             
    dTv_max0 =                       2
    dTv_sigma =                    0.5 
    Nv_Nav =  nv_max/(math.sqrt(2*math.pi))*math.e**(-(vv-dTv_max0)**2/(2*dTv_sigma**2))
    Nv_Nav=0
    return Nv_Nav

def nakljucje(msm):
   rand = np.random.random_sample(msm.shape)
   rand[taula==-1]=-1  
   return rand

def taljenje(u, Tm):
   global taula
   plt.clf()
   taula=np.zeros((Z,X,Y))
   taula[u<(Tm)] = -1
   return taula

# **** Growth functions ****
def liquidus(u, liquidus_temp):
   global likvid
   likvid=np.zeros((Z,X,Y))
   likvid[u>liquidus_temp]=-1
   return likvid

def growth_speed(dT_field):
    vg = 2.03e-4 * dT_field**2 - 0.544e-4 * dT_field
    return vg


# **** Grain orientation and Cube Rotation Matrices ****

def Rotate_the_Cube_XY(alpha, beta):              # rotation of cube around Y- axis (alpha angle) and  X-axis (beta angle)
    global RmX, RmY
    cube = np.array([
                              [1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1]
                            ])
    RmY = np.array([ # rotation matrix around Y-axis
        [math.cos(math.radians(alpha)), 0, math.sin(math.radians(alpha))],
        [0, 1, 0],
        [-math.sin(math.radians(alpha)), 0, math.cos(math.radians(alpha))],
        ])
    RmX = np.array([ # rotation matrix around X-axis
        [1, 0, 0],
        [0, math.cos(math.radians(beta)), math.sin(math.radians(beta))],
        [0, - math.sin(math.radians(beta)), math.cos(math.radians(beta))],
        ])
    oix = multi_dot([cube, RmY, RmX])
    return oix

def Rotate_the_Cube_Z(xy_cube, gamma):                   # rotation of cube around Z- axis (gamma angle)
    RmZ = np.array([# rotation matrix around Z-axis
        [math.cos(math.radians(gamma)), math.sin(math.radians(gamma)), 0],
        [- math.sin(math.radians(gamma)), math.cos(math.radians(gamma)), 0],
        [0, 0, 1],
        ])
    oiz = np.dot(xy_cube, RmZ)
    oiz_switched_rows_and_columns = oiz[np.ix_([2,0,1], [2,0,1])]
    return oiz_switched_rows_and_columns

def get_color(x,y):
    u = 1-x
    v = x-y
    w = y
    rgb = np.array([u,v,w])
    RGB = 255*rgb/np.max(rgb)
    return RGB.astype('int')

def random_xy_tilt(posibilities):
    xran1 = np.random.randint(0, posibilities+1)
    xran2 = np.random.randint(0, posibilities+1)
    xran3 = np.random.randint(0, posibilities+1)
    xran = (xran1+xran2+xran3)/3  
    x = xran / posibilities
    y = np.random.randint(0, xran+1)/posibilities
    return x,y

def merging(prva,druga):
    tot = np.where(druga==0, prva, druga)
    return tot

def time_counter(attr0, neg_index):      # attr0 = manager.(Negatives)
        attr0[neg_index] += dt
        return attr0
        
def random_selection(x, erase):
    global item
    if len(x)!=0:
        item = random.choice(x)
        if erase:
            x=np.delete(x,np.where(x==item))
    else:
        pass
    return x

def Ordered_Grains_n_Directions(g_ID, smeri_arrays):
    selekcija=[]
    for gid in g_ID:
        for sa in smeri_arrays:
            for sm in sa:
                selekcija.append((gid, sm))
    return selekcija

def Random_Grains_Ordered_Directions(g_ID, smer):
    selekcija=[]
    while True:
        if any(g_ID):
            g_ID = random_selection(g_ID, True)
            grain = item
            for j in smer:
                selekcija.append((grain,j))
        else:
            break     
    return selekcija

def Random_Grains_n_Directions(g_ID, smer):
    GD = {}
    for gid in g_ID:
       GD[gid] = smer.copy()
    selekcija=[]
    while True:
        if any(g_ID):
            g_ID = random_selection(g_ID, False)
            grain = item
            counter = 0
        else:
            break
        while counter==0:
            if any(GD[grain]):
                GD[grain] = random_selection(GD[grain], True)
                s = item
                selekcija.append((grain,s))
                counter+=1
                break
            else:
                g_ID=g_ID[g_ID!=grain]
                break
        continue
    return selekcija

def Random_Grains_Directions_Segmented(g_ID, smeri_arrays):
    SM=[]
    for sm in smeri_arrays:
        SM+=Random_Grains_n_Directions(g_ID, sm)
    return SM

def TT(x):
   temp = np.load(directory.get()+'salome_'+str(x)+'.npy')[z_min:z_max, x_min:x_max, y_min:y_max]
   return temp

def W(O, L):
    O = np.array(O)
    cos_Z = np.dot(O[0], L) / (np.linalg.norm(O[0]* np.linalg.norm(L)))
    cos_X = np.dot(O[1], L) / (np.linalg.norm(O[1] * np.linalg.norm(L)))
    cos_Y = np.dot(O[2], L) / (np.linalg.norm(O[2] * np.linalg.norm(L)))
    return np.max(np.absolute(np.array([cos_Z, cos_X, cos_Y])))

def Dissipate_Distances(eps):
    R_min  = 0.5 - epsilon
    R_max = 1-R_min
    RAN = np.random.uniform(R_min, R_max, (Z,X,Y))
    return RAN

def shape_double(nuum, zao):
    celi = int(nuum)
    if celi < 10:
        jh='%.'+str(zao)+'f'
        ven = jh % round(nuum, zao)
    elif 10 <= celi < 100 :
        jh='%.'+str(zao-1)+'f'
        ven = jh % round(nuum, zao-1)
    elif 100 <= celi < 1000:
        jh='%.'+str(zao-2)+'f'
        ven = jh % round(nuum, zao-2)
    elif 1000 <= celi:
        if zao>3:
            jh='%.'+str(zao-3)+'f'
            ven=jh % round(nuum, zao-3)
        else:
            ven=celi
    return ven

def plot_CPU_i():
    plt.figure('CPU plot vs. Time Step')
    plt.plot(PD['h'], PD['current step CPU'], marker='o')
    plt.xlabel('Time step'); plt.ylabel('CPU time (sec.)')

def plot_ALLgrains_i():
    plt.figure('ALL grains plot vs. Time Step')
    plt.plot(PD['h'], PD['ALL grains'])

# Writing the log..
def Writing_the_log():
    with open(directory.get()+track+'Logfile.txt', 'w')as cuttxt:
        cuttxt.write(100*'*'+'\n'+'z_min = '+str(z_min)+' ,  z_max = '+str(z_max)+' ,  Z = '+str(Z)+' ,\n'+
         'X = '+str(X)+' ,  y_min = '+str(y_min)+' ,  y_max = '+str(y_max)+' ,  Y = '+str(Y)+' ,\n'+
         'N = '+str(N)+' ,\n\n'+

         'FEM_time_step = '+str(FEM_time_step)+' sec.'+' ,\n'+
         'dt = '+str(dt)+' sec. ,  FEM_scanning_speed = '+str(FEM_scanning_speed)+' mm/s'+' ,\n'+
         'FEM_cell_size = '+str(1e6*FEM_cell_size)+u'\u00B5m ,  cell = '+str(1e6*cell)+u'\u00B5m'+' ,\n'+
         'space_factor = '+str(space_factor)+' ,  FEM_time_factor = '+str(FEM_time_factor)+' ,  extra_time_factor = '+str(extra_time_factor)+' ,\n\n'+

         'START_step = '+str(START_step)+' ,\n\n'+

         'smeri = '+str(list(smeri))+' ,\n\n'+

         'BETA = '+str(BETA)+' ,\n\n'+

         'Tmelt_Celsius = '+str(Tmelt_Celsius)+u'\N{DEGREE SIGN}C ,  dTliquidus = '+str(dTliquidus)+u'\N{DEGREE SIGN}C'+' ,\n'+
         'delete_inactive_grains_ID = '+str(delete_inactive_grains_ID)+'\n\n'+
         100*'*'+'\n\n')

def Save_KickOff():
    kickoff_folder = 'kickoff_data/'
    if not os.path.isdir(directory.get()+track+kickoff_folder):
       os.mkdir(directory.get()+track+kickoff_folder)

    np.save(directory.get()+track+kickoff_folder+'faza_kickoff.npy', faza)
    np.save(directory.get()+track+kickoff_folder+'cas_kickoff.npy', cas)
    np.save(directory.get()+track+kickoff_folder+'rgb_snap_kickoff.npy', rgb_snap)
    try:
        np.save(directory.get()+cuts_RGB+'cut_RGB_'+str(cut_count)+'.npy', rgb_snap[:,:,:cutoff_limit,:])
        np.save(directory.get()+cuts_RGB+'cut_RGB_'+str(cut_count+1)+'.npy', rgb_snap[:,:,cutoff_limit:,:])
        np.save(directory.get()+cuts_faza+'cut_faza_'+str(cut_count)+'.npy', faza[:,:,:cutoff_limit])
        np.save(directory.get()+cuts_faza+'cut_faza_'+str(cut_count+1)+'.npy', faza[:,:,cutoff_limit:])
    except FileNotFoundError:
        os.mkdir(directory.get()+cuts_RGB)
        os.mkdir(directory.get()+cuts_faza)
        np.save(directory.get()+cuts_RGB+'cut_RGB_'+str(cut_count)+'.npy', rgb_snap[:,:,:cutoff_limit,:])
        np.save(directory.get()+cuts_RGB+'cut_RGB_'+str(cut_count+1)+'.npy', rgb_snap[:,:,cutoff_limit:,:])
        np.save(directory.get()+cuts_faza+'cut_faza_'+str(cut_count)+'.npy', faza[:,:,:cutoff_limit])
        np.save(directory.get()+cuts_faza+'cut_faza_'+str(cut_count+1)+'.npy', faza[:,:,cutoff_limit:])

    with open(directory.get()+track+kickoff_folder+'nuclei_kickoff.json', 'w') as nuks:              # Writing data of ALL nuclei as .json file, but values must be list NOT np.array !!!
        asc_list = asc.copy()
        for nuk in asc:
            try:
                asc_list[nuk]['oi']=asc[nuk]['oi'].tolist()
                asc_list[nuk]['rgb']=asc[nuk]['rgb'].tolist()
            except AttributeError:
                asc_list[nuk]['oi']=asc[nuk]['oi']
                asc_list[nuk]['rgb']=asc[nuk]['rgb']
        json.dump(asc_list, nuks)
    
    with open(directory.get()+track+kickoff_folder+'negatives_kickoff.json', 'w') as negs:
        json.dump(Negatives, negs)
    with open(directory.get()+track+kickoff_folder+'S_kickoff.json', 'w') as S_kick:
        S_list = S.copy()
        for _s_ in S:
            try:
                S_list[_s_]=S[_s_].tolist()
            except AttributeError:
                S_list[_s_]=S[_s_]
        json.dump(S_list, S_kick)
    
    np.save(directory.get()+track+kickoff_folder+'grain_ID_kickoff.npy', grain_ID_)
    np.save(directory.get()+track+kickoff_folder+'FF_kickoff.npy', FF)
    np.save(directory.get()+track+kickoff_folder+'inactive_grains_kickoff.npy', inactive_grains)

    with open(directory.get()+track+kickoff_folder+'AG_kickoff.json', 'w') as ag:
        json.dump(AG, ag)
    with open(directory.get()+track+kickoff_folder+'IG_kickoff.json', 'w') as ig:
        json.dump(IG, ig)

    counts={}; counts['tm_count']=tm_count ; counts['yp_count']=yp_count ; counts['cut_count']=cut_count; counts['start_step_intermediate'] = i
    counts['h_intermediate']=h-1
    with open(directory.get()+track+kickoff_folder+'counts_kickoff.json', 'w') as cnt:
        json.dump(counts, cnt)
    with open(directory.get()+track+kickoff_folder+'PD.json', 'w') as pd:
        json.dump(PD, pd)
        
def Load_KickOff():
    kickoff_folder = 'kickoff_data/'
    faza=np.load(directory.get()+track+kickoff_folder+'faza_kickoff.npy')
    rgb_snap=np.load(directory.get()+track+kickoff_folder+'rgb_snap_kickoff.npy')
    cas= np.load(directory.get()+track+kickoff_folder+'cas_kickoff.npy')
    with open(directory.get()+track+kickoff_folder+'nuclei_kickoff.json', 'r') as nuks:              # Writing data of ALL nuclei as .json file, but values must be list NOT np.array !!!
        asc=json.load(nuks)
        asc ={int(k):v for k,v in asc.items()}
        for nuk in asc:
            asc[nuk]['oi']=np.array(asc[nuk]['oi'])
            asc[nuk]['rgb']=np.array(asc[nuk]['rgb'])

    grain_counter = len(list(asc.keys()))

    grain_ID_ = np.load(directory.get()+track+kickoff_folder+'grain_ID_kickoff.npy')
    FF = list(np.load(directory.get()+track+kickoff_folder+'FF_kickoff.npy'))
    inactive_grains = np.load(directory.get()+track+kickoff_folder+'inactive_grains_kickoff.npy').tolist()
           
    with open(directory.get()+track+kickoff_folder+'negatives_kickoff.json', 'r') as negs:
        Negatives=json.load(negs)
        Negatives ={int(k):v for k,v in Negatives.items()}
    with open(directory.get()+track+kickoff_folder+'S_kickoff.json', 'r') as S_kick:
        S=json.load(S_kick)
        for _s_ in S:
            S[_s_]=np.array(S[_s_])
    with open(directory.get()+track+kickoff_folder+'AG_kickoff.json', 'r') as ag:
        AG=json.load(ag)
        AG ={int(k):v for k,v in AG.items()}
    with open(directory.get()+track+kickoff_folder+'IG_kickoff.json', 'r') as ig:
        IG=json.load(ig)
        IG ={int(k):v for k,v in IG.items()}
    with open(directory.get()+track+kickoff_folder+'counts_kickoff.json', 'r') as cnt:
        counts = json.load(cnt)
    with open(directory.get()+track+kickoff_folder+'PD.json', 'r') as pd:
        PD = json.load(pd)
    
    return faza, rgb_snap, cas, asc, Negatives, grain_counter, S, AG, IG, FF, inactive_grains, counts['tm_count'], counts['yp_count'], counts['cut_count'], counts['start_step_intermediate'], grain_ID_, counts['h_intermediate'], PD

def Selection_Mechanism(zrno_ID, smeri_podatkovna, pick_one):
    selection_mechanisms ={
            1: 'Ordered_Grains_n_Directions(zrno_ID, smeri_podatkovna)',                                           # fully ORDERED ........OK
            2: 'Ordered_Grains_n_Directions(zrno_ID, np.flip(smeri_podatkovna))',                            # fully ORDERED, flipped

            3: 'Random_Grains_Ordered_Directions(zrno_ID,  smeri)',                                                    # random selection of grains, ordered selection of directions, from low to high order ......... not ok :( 
            4: 'Random_Grains_Ordered_Directions(zrno_ID, np.flip( smeri_podatkovna.flatten()))', # random selection of grains, ordered selection of directions, from high to low order ......... OK
            
            5: 'Random_Grains_Directions_Segmented(zrno_ID, smeri_podatkovna)',                             # random selection of grains, ordered selection of segments of random directions within, from low to high order ......... not ok :(
            6: 'Random_Grains_Directions_Segmented(zrno_ID, np.flip(smeri_podatkovna))',              # random selection of grains, ordered selection of segments {[001,010,011],  [012,021],  [002,020],  [022]}of random directions within, from high to low order ......... OK

            7: 'Random_Grains_n_Directions(zrno_ID, smeri_podatkovna.flatten())', }                       # fully RANDOM ....... not OK

    selekcija = eval(selection_mechanisms[pick_one])
    return selekcija

def make_mpi_PARTITIONS_Y(faza, cas, vg, S, taula, T_next, T, likvid, begin, end):
    #global faza, cas, vg, S, taula, T_next, T, likvid

    faza_py = faza[:,:,begin:end]                         # faza ::  partition (in Y-direction)
    cas_py=cas[:,:,begin:end]                             # cas ::  partition (in Y-direction)
    vg_py=vg[:,:,begin:end]
    S_py={}
    for sas in S:
        S_py[sas]=S[sas][:,:,begin:end]
    taula_py=taula[:,:,begin:end]
    T_next_py=T_next[:,:,begin:end]
    T_py=T[:,:,begin:end]
    likvid_py=likvid[:,:,begin:end]
    
    return faza_py, cas_py, vg_py, S_py, taula_py, T_next_py, T_py, likvid_py

def make_shared_memory_Array(tip, shz, shx, shy):
    if tip == 'integer':
        sharr = Array('i', shz*shx*shy)
        sharr = np.frombuffer(sharr.get_obj(), ctypes.c_int)
    elif tip == 'double':
        sharr = Array('d', shz*shx*shy)
        sharr = np.frombuffer(sharr.get_obj(), ctypes.c_double)
    sharr = sharr.reshape((shz, shx, shy))
    return sharr


'''======================================================================================================================================================================'''
'''                                                     Moore Neighbourhood - Crystallographic Orientations, Directions & GROUPS                                                                                                                                                    '''
'''======================================================================================================================================================================'''
'''........................................................................................................................................... I. order Moore neighbourhood '''
G1  =  np.array(['001', '00_1', '010', '0_10']) 
G2  =  np.array(['100']) 
G3  =  np.array(['_100']) 
G4  =  np.array(['011', '01_1', '0_11', '0_1_1'])
G5  =  np.array(['101', '10_1', '110', '1_10'])
G6  =  np.array(['_101', '_10_1', '_110', '_1_10'])
G7  =  np.array(['111', '11_1', '1_11', '1_1_1'])
G8  =  np.array(['_111', '_11_1', '_1_11', '_1_1_1'])

'''.......................................................................................................................................... II. order Moore neighbourhood '''
G9  =  np.array(['012', '01_2', '0_12', '0_1_2', '021', '02_1', '0_21', '0_2_1'])
G10 =  np.array(['102', '10_2', '120', '1_20'])
G11 =  np.array(['_102', '_10_2', '_120', '_1_20'])
G12 =  np.array(['201', '20_1', '210', '2_10'])
G13 =  np.array(['_201', '_20_1', '_210', '_2_10'])
G14 = np.array(['002', '00_2', '020', '0_20'])
G15 = np.array(['200'])
G16 = np.array(['_200'])
G17 = np.array(['022', '02_2', '0_22', '0_2_2'])
G18 = np.array(['202', '20_2', '220', '2_20'])
G19 = np.array(['_202', '_20_2', '_220', '_2_20'])
G20 = np.array(['112', '11_2', '1_12', '1_1_2',    '121', '12_1', '1_21', '1_2_1'])
G21 = np.array(['_112', '_11_2', '_1_12', '_1_1_2',    '_121', '_12_1', '_1_21', '_1_2_1'])
G22 = np.array(['211', '21_1', '2_11', '2_1_1'])
G23 = np.array(['_211', '_21_1', '_2_11', '_2_1_1'])
G24 = np.array(['122', '12_2', '1_22', '1_2_2'])
G25 = np.array(['_122', '_12_2', '_1_22', '_1_2_2'])
G26 = np.array(['212', '21_2', '2_12', '2_1_2',    '221', '22_1', '2_21', '2_2_1'])
G27 = np.array(['_212', '_21_2', '_2_12', '_2_1_2',    '_221', '_22_1', '_2_21', '_2_2_1'])
G28 = np.array(['222', '22_2', '2_22', '2_2_2'])
G29 = np.array(['_222', '_22_2', '_2_22', '_2_2_2'])


'''  M P I  $$$$$$$$$$$$$$$$$$$  MPI  tryin' to make the code run on several cores using multiprocessing module  MPI  $$$$$$$$$$$$  M P I  '''

from multiprocessing import Pool, Manager, Process, Array
import ctypes
from functools import partial
from inspect import signature
import itertools
import MPI_pool_functions as mpf       # MPI pool functions

''' Instruction on HOW TO USE PARTIAL FUNCTION  ---> check the MPI_pool_toycode.py '''

""" ====================================  E X E C U T I O N  ========================================================"""

def ALL_variables_refresh():
    global cell, FEM_time_step, time_factor, dt, time_factor_folder, flashies_RGB, flashies_faza, cuts_RGB, cuts_faza, taula, likvid
    global Z,X,Y, y_min,y_max, yp_list, tm_list, Tmelt, yp_count, tm_count, random_distances, R, i, start_time, PD, START_step, T, T_next
    global vg, NP, NV, faza, cas, Negatives, IG, AG, FF, cut_count, h, S, asc, rgb_snap, grain_counter, inactive_grains, grain_ID_, grain_ID
    global Ycut_indices, pool, Q1, Q2, manager

    Safe_Refresh_of_FEM_parameters(None, False)
    Safe_Refresh_of_CA_parameters(None, False)
    
    cell = FEM_cell_size/space_factor          
    FEM_time_step = FEM_cell_size * 1000 / FEM_scanning_speed           #  time_step to define loads in FEM; unit: SECOND [s]
    time_factor =     extra_time_factor * FEM_time_factor
    dt = FEM_time_step / extra_time_factor               #  CA time differential (after space & time interpolation); unit: SECOND [s]
    Tmelt= Tmelt_Celsius + 273                                    #   Melting point; unit:  KELVIN [K]

    if not os.path.isdir(directory.get()+track):
        os.mkdir(directory.get()+track)

    time_factor_folder = '/time_factor_24/time_factor_3'

    flashies_RGB =       track+'flashies_RGB/'                          #  Subfolder with time-snap 3D matrices, i.e. flashies
    flashies_faza =        track+'flashies_faza/'

    cuts_RGB =            track+'cuts_RGB/'                               #  Subfolder with cut 3D matrices, i.e. cuts
    cuts_faza =             track+'cuts_faza/'

    ''' ~~~~~~~~~~~~~~~~ Long_track_CA  ::: Domain Constraints ::: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~''' 

    Z =                  z_max - z_min

    X=                  96                     #  [27-15] =  12 cells before space interpolation, and 12*8 = 96 cells after space interpolation
    Y=                 128                   #  N=12 ; length of two YPs (pair), each 120 cells long: YP0 [12,27] --->>> 27-12 = 15*8 = 120 + 8
    #Y =                32                   #   N=80 ;  

    ''' ~~~~~~~~~~~~~~~~ Making of  ::: Time (tm_list):::  and  ::: Space (yp_list):::  Partitions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
    #y_min =           16       # N = 80 
    y_min =             12       # N = 12
    y_max =            97

    N =                    12         # Number of equally sized segments along Y-axis of 4D matrix

    YP = Make_Y_Partitions(y_min, y_max-1, N)
    yp_list = [i+'  '+str(YP[i]).replace(' ', '')for i in YP]                                  #  ::: Creation of Y partitions names ::: list of strings  (yp_list)
    tm_list = ['TR{0}  [{0},{1}]'.format(i,i+1) for i in range (0,16)]         #  ::: Creation of time ranges names ::: list of strings  (tm_list)
    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

    Select_Groups(moore_var.get())

    if epsilon == 0:
       random_distances = False
    elif epsilon > 0 and epsilon <= 0.49:
       random_distances = True
    R=Dissipate_Distances(epsilon)if random_distances else 1
    
    #  ~~~~ AUTO Nucleation (AN)~~~~
    if avtomatska_nukleacija.get():
       vg = np.zeros((Z,X,Y))        # matrika hitrosti rasti
       NP = np.vectorize(nukleacija_povrsina)
       NV = np.vectorize(nukleacija_volumen)
       if from_beginning.get():
          faza = np.zeros((Z,X,Y))  # fazna matrika
          cas = np.zeros((Z,X,Y))   # časovna matrika
          Negatives = {-1:0}; asc = {}; grain_counter = 0; inactive_grains=[]; S={}; grain_ID=[] ; grain_ID_=[]
          IG = {}     #  IG (saves time-steps when indivudual grain didn't grow IN ANY DIRECTION, IG stands for Inactive Grains)
          AG = {}   #  AG (saves time-steps when indivudual grain did grow IN ANY DIRECTION, AG stands for Active Grains)
          FF = []
          PD = {'h': [], 'current step CPU': [], 'ALL grains': []}  #   Plot Database
          i = START_step-1
          Writing_the_log()
          
       elif not from_beginning.get():
          faza, rgb_snap, cas, asc, Negatives, grain_counter, S, AG, IG, FF, inactive_grains, tm_count, yp_count, cut_count, START_step, grain_ID_, h, PD  = Load_KickOff()
          i = START_step
          #if negind<=critical_negind:
                #Negatives ={key:val for key, val in Negatives.items() if key < (negind+negatives_thresh)}
          #Negatives ={key:val for key, val in Negatives.items() if val <  dt_thresh}

    #  ~~~~ Manual Nucleation (MN)~~~~
    elif not avtomatska_nukleacija.get():
        print('NOT auto nucleation')
        '''||||||||||||||||||||||||||| nucleation - manual |||||||||||||||||||||||||| '''
        Z,X,Y = 1, 101, 101                                              # Size of domain in terms of cells in Z,X,Y directions, respectively, for testing and code development
        faza = np.zeros((Z,X,Y))                               # fazna matrika
        cas = np.zeros((Z,X,Y))                                # časovna matrika
        vg =   1                                                                # Value of homogeneous 'growing velocity' field, a.u., for testing and code development
        cell =  1                                                               # for mesh dependency development (MDD)
        dt = cell/8
        T=1; T_next=0
        Negatives = {-1:0};  S = {}; asc = {}; grain_counter = 0; inactive_grains=[]; FF = []
        #critical_negind = - 2000
        #negatives_thresh = 15

        M = { 1: {'BETA':(0, 50, 50), 'Ł': (0,0,0)},
                   2: {'BETA':(0, 40, 40), 'Ł': (0,0,45)},
                   3: {'BETA':(0, 60, 60), 'Ł': (1,1,0)},                    #  data of manually created nuclei, 'BETA' are the (Z,X,Y) coordinates, 'Ł' are tilting parameters (x,y,gama)
                   4: {'BETA':(0, 50, 30), 'Ł': (0,0,45*0.5)},
                   5: {'BETA':(0, 32, 42), 'Ł': (0,0,45*0.75)},
                   #6: {'BETA':(0, 150, 150), 'Ł': (0,0,45*0.875)},
                   #7: {'BETA':(0, 175, 175), 'Ł': (0,0,45)},
                 }

        for i in M:
            faza[M[i]['BETA'][0],M[i]['BETA'][1],M[i]['BETA'][2]]=i                               # define nucleus ID in faza matrix
            x,y = M[i]['Ł'][0], M[i]['Ł'][1]
            rgb = get_color(x,y)
            alfa  = math.degrees(math.atan(x))- 9.7356103173*x*y
            beta = math.degrees(math.atan(y))- 9.7356103173*x*y
            cub_xy = Rotate_the_Cube_XY(alfa, beta)
            gama = M[i]['Ł'][2]
            oi = Rotate_the_Cube_Z(cub_xy, gama)
            asc[i] ={'oi': oi, 'alfa':alfa, 'beta':beta, 'gama':gama, 'rgb':rgb,}

        grain_ID = np.array(list(asc.keys()))
        Selection = Selection_Mechanism(grain_ID, smeri_database, pick_selection_mechanism)
        cas[np.isin(faza, grain_ID, invert=False)] = -1
        taula=0;likvid=0; grain_counter=len(grain_ID)
        PD = {'h': [], 'current step CPU': [], 'ALL grains': []}  #   Plot Database

    start_time = time.time()


    if MPI:
        Ycut_indices = [0, int(Y/2), Y]            # splitting basic matrix into 2 partitions of equal size in Y-direction
        if __name__ == '__main__':
    
            #set_start_method('spawn')
            #manager=Manager()
            #Q1=manager.Queue()
            #Q2=manager.Queue()
            pool = Pool()
   
            faza = make_shared_memory_Array('integer', Z,X,Y)
            cas = make_shared_memory_Array('double', Z,X,Y)
            vg = make_shared_memory_Array('double', Z,X,Y)
            taula = make_shared_memory_Array('double', Z,X,Y)
            T_next = make_shared_memory_Array('double', Z,X,Y)
            T = make_shared_memory_Array('double', Z,X,Y)
            likvid = make_shared_memory_Array('double', Z,X,Y)
            #global faza, cas, vg, S, taula, T_next, T, likvid



            #manager = Manager()
            #faza = manager.list(faza)                               
            #cas = manager.list(cas)                                 
            #Negatives = manager.dict({-1:0})
            #S = manager.dict()
            #asc = manager.dict(asc)
            #if avtomatska_nukleacija.get():
                #pass
                #vg = manager.list(vg)
                #FF = manager.list(FF)

    elif not MPI:
        pool=None

 
ALL_variables_refresh()


""" =====================================  I  N  P  U  T    D  A  T  A  ======================================================================="""
# ````` DOMAIN Input variables `````

GDSM = 'new'        #  choose  'old'  or   'new'   , Grain and Directions Selecting Mechanism (GDSM)- see explanation below:
'''
GDSM explanation:

GDSM OLD =    Picks random grain and performs growth in ALL directions by merging function. Here a list
                        seznam premikov has the length of number of given directions.

GDSM NEW =   Picks random grain and performs growth in ONE direction, which can be choosed either in
                        order or randomly. Then, matrix faza is refreshed and the process is repeated for the rest of
                        the directions. '''


'''*************************  FIXED VARIABLES  *********************************************'''
h =                 0                                 #  Relative time step counter
cut_count =   0
Cut_Off_Percent =   50                    #  Percent of Stack_2 domain lenght at which this domain should be cut off  ---> the left part is saved (.npy and/or .png)the right goes on to CA and so on and on and on..
cutoff_limit =  64

FF_length = 30     #  Inactive grains are deleted every FF_length time step
# Filter Negatives by keys
critical_negind = - 100
#negatives_thresh = 500
# Filter Negatives by values
dt_thresh = 500*dt
'''************************************************************************************************'''

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ... STARTING the CA ..
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Growth (time FOR-LOOP, grains FOR-LOOP, directions FOR-LOOP)~~~~~~~~~~~~~~~~~~~~~~~~~
def Run_CA(p):
    global rgb_snap, tm_count, yp_count, cut_count, run_ca, F, TC_list, nega, GGL, i, MM_list, q, h, grain_ID, grain_ID_
    global Negatives, faza, cas, S, FF, vg, T_next, T, grain_counter, GDSM,  asc,  taula, likvid, R, cell, W,  random_distances, Selection
    global partition_list, p1
    i+=1
    step_time_start = time.time(); h+=1
    # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   1.   N  U  C  L  E  A  T  I  O  N   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    
    if avtomatska_nukleacija.get():
        ''' avtomatska nukleacija '''

        T = Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], i)[z_min:z_max]
        taljenje(T, Tmelt)                              # condition for melting  ----> taula matrix of cell values;  if value -1 then solid (or powder)elif value 0 then liquid
        
        if np.all(taula[:,:,:cutoff_limit]==-1):   # CUT condition and consequences
            print(); print(50*'*',' CUT! ',50*'*')
            print('Cut_Off_Percent = ',Cut_Off_Percent,' %    , cutoff limit = ',cutoff_limit)
            cut_text = 'Time step number  '+str(i)+',  real time: '+str(round(1000*dt*h, 3))+' msec.'
            print(cut_text)
            print(106*'*'); print()

            with open(directory.get()+track+'cut_data.txt', 'a')as cuttxt:
                cuttxt.write(cut_text+'\n')

            if save_cut_as_RGB.get():
                try:
                    np.save(directory.get()+cuts_RGB+'cut_RGB_'+str(cut_count)+'.npy', rgb_snap[:,:,:cutoff_limit])       # Saves the first half of RGB snap as .NPY
                except FileNotFoundError:
                    os.mkdir(directory.get()+cuts_RGB)
                    np.save(directory.get()+cuts_RGB+'cut_RGB_'+str(cut_count)+'.npy', rgb_snap[:,:,:cutoff_limit])
            if save_cut_as_faza.get():
                try:
                    np.save(directory.get()+cuts_faza+'cut_faza_'+str(cut_count)+'.npy', faza[:,:,:cutoff_limit])                 # Saves the first half of faza snap as .NPY
                except FileNotFoundError:
                    os.mkdir(directory.get()+cuts_faza)
                    np.save(directory.get()+cuts_faza+'cut_faza_'+str(cut_count)+'.npy', faza[:,:,:cutoff_limit])
                    
            yp_count+=1
            cut_count+=1
            
            faza = np.dstack((faza[:,:,cutoff_limit:],np.zeros((Z,X,cutoff_limit) )))
            for s in smeri:
                S[s] = np.dstack((S[s][:,:,cutoff_limit:],np.zeros((Z,X,cutoff_limit) )))
            FF=[]
            cas = np.dstack((cas[:,:,cutoff_limit:],np.zeros((Z,X,cutoff_limit) )))
                
            T = Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], i)[z_min:z_max]
            taljenje(T, Tmelt)

            if save_kickoff.get():
                Save_KickOff()
                    
        liquidus(T,Tmelt+dTliquidus)              # absolute liquidus line
        
        dTt =  T  -  Tmelt                                  # undercooling [K]
        
        interface = NP(dTt)                                              
        bulk = NV(dTt)
        live= nakljucje(taula)

        try:
            T_next = Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], i+time_shift)[z_min:z_max]
        except FileNotFoundError:
            print(); print('FileNotFound Exception !'); print()
            tm_count+=1
            T_next = Stack_2(yp_list[yp_count], yp_list[yp_count+1], tm_list[tm_count], i+time_shift)[z_min:z_max]
        except IndexError:
            raise IndexError('Please, correct the time range!')

        ''' ======================================= N  U  C  L  E  A  T  I  O  N ============================================================== '''
        new_grains_ID = []
        
        for k in range(Z):
         for ii in range(X):
            for j in range(Y):
               if faza[k][ii][j]==0 and (BETA<live[k][ii][j]<interface[k][ii][j] or (live[k][ii][j]<bulk[k][ii][j] and bulk[k][ii][j]>BETA)) and T_next[k][ii][j] < T[k][ii][j]:              
               
                  grain_counter +=1
                  new_grains_ID.append(grain_counter)
                  #IG[grain_counter]=[]
                  #AG[grain_counter]=[]
                  faza[k][ii][j]=grain_counter
                  """ generation of random grain orientation """
                  x,y = random_xy_tilt(rp)
                  rgb = get_color(x,y)
                  alfa  = math.degrees(math.atan(x))- 9.7356103173*x*y
                  beta = math.degrees(math.atan(y))- 9.7356103173*x*y
                  cub_xy = Rotate_the_Cube_XY(alfa, beta)
                  gama = math.degrees(math.atan(np.random.randint(0, rp+1)/rp))
                  oi = Rotate_the_Cube_Z(cub_xy, gama)
                  asc[grain_counter] ={'oi': oi, 'alfa':alfa, 'beta':beta, 'gama':gama, 'rgb':rgb, 'coords':(i,k,ii,j), 'temp': T[k,ii,j]-273, }    # ALL data about nuclei
        ''' ========================================================================================================================== '''
        try:
            Selection = Selection_Mechanism(grain_ID_, smeri_database, pick_selection_mechanism)   
            grain_ID = grain_ID_.copy()#; print('grain_ID after nucleation (in try): ',grain_ID,'   (copy of grain_ID_, which is ',grain_ID_,')')
        except NameError:
            grain_ID = np.array(list(asc.keys()))#; print('grain_ID after nucleation (in NameError exception): ',grain_ID)
            Selection = Selection_Mechanism(grain_ID, smeri_database, pick_selection_mechanism) 
        
        cas[np.isin(faza, new_grains_ID, invert=False)] = -1          # vrednost časovne matrike vseh novih nukleusov je -1
        vg = growth_speed(dTt)

    elif not avtomatska_nukleacija.get():
        Selection = Selection_Mechanism(grain_ID, smeri_database, pick_selection_mechanism)
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 

    # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   2.   T  I  M  E  -  MPI   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    
    ''' ......................................... 1st RUN --- time counter PART . . . '''
    if run_1st:
        run1_time_start = time.time()
        #if not MPI:
        if True:
            for negind in Negatives:
                time_counter_original(negind)
        '''
        elif MPI:
            nega = sorted(list(Negatives.keys()))
            TC_list = p.map(partial(time_counter, Negatives), nega)   # --------> returns :: Negatives
            Negatives = TC_list[-1]
        '''
        run1_time_end = time.time()
        print('1st run time:  ',round(run1_time_end-run1_time_start, 4),' sec.')
           
    ''' ......................................... 2nd RUN --- MM_ multiprocessing pool function PART . . . '''
    if run_2nd:
        run2_time_start = time.time()
        #if not MPI:
        if True:
            q=list(itertools.product(smeri, Negatives))
            for qq in q:
                S=mpf.MM(Negatives, faza, grain_ID, cas, S, qq)
        ''' 
        elif MPI:
            nega = sorted(list(Negatives.keys()))
            q=list(itertools.product(smeri, nega))
            MM_list = p.map(partial( mpf.MM, Negatives, faza, grain_ID, cas, S), q)   # --------> returns :: S
            for j in MM_list:
                S.update(j)
        '''
        run2_time_end = time.time()
        print('2nd run time:  ',round(run2_time_end-run2_time_start, 4),' sec.')

        
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 

    # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   3.   P  H  A  S   E   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    ''' ......................................... 3rd RUN --- GG_ multiprocessing pool function PART :: execution . . . '''
    if run_3rd:
        if not MPI:
            F = faza.copy()
            #FF.append(faza)
            for selsel in Selection:
                    faza, cas = mpf.GG(GDSM, faza, asc, vg, S, taula, T_next, T, likvid, R, cell, W, cas, random_distances, selsel)
            
            if not np.all(faza==F):
                #negind -=1
                negind = min(list(Negatives.keys()))- 1             #  Here, negind should be defined as absolute counter of negative numbers to avoid erasing existing Negatives items
                for s in S:
                    S[s][np.isin(faza-F, grain_ID)]= negind
                Negatives[negind]=0
                
                ''' Deleting Inactive Items from  ::  Negatives '''
                #negatives_thresh=100
                #Negatives ={key:val for key, val in Negatives.items() if key < (negind+negatives_thresh)}

    
        elif MPI:
            F=faza.copy()
            
            if avtomatska_nukleacija.get()and bool(Selection):

                partlist = []
                #if h==1:
                    #faza = manager.list(faza)                               
                    #cas = manager.list(cas) 

                ''' NO partitions '''
                partlist.append(make_mpi_PARTITIONS_Y(faza, cas, vg, S, taula, T_next, T, likvid, None, None))

                
                
                ''' PARTITION # 1 '''
                #faza1, cas1 = make_mpi_PARTITIONS_Y(Ycut_indices[0], Ycut_indices[1]+1)
                #partlist.append(make_mpi_PARTITIONS_Y(faza, cas, vg, S, taula, T_next, T, likvid, Ycut_indices[0], Ycut_indices[1]+1))
                ''' PARTITION # 2 '''
                #faza2, cas2 = make_mpi_PARTITIONS_Y(Ycut_indices[1]-1, Ycut_indices[2])
                #partlist.append(make_mpi_PARTITIONS_Y(faza, cas, vg, S, taula, T_next, T, likvid, Ycut_indices[1]-1, Ycut_indices[2]))
                
                GGL = pool.map(partial( mpf.GG_pool, GDSM, asc, R, cell, W, random_distances, Selection), [partlist[0],])
                #GGL = pool.map(partial( mpf.GG_pool, GDSM, asc, R, cell, W, random_distances, Selection), [partlist[0], partlist[1]])

                #GGL = pool.map(partial( mpf.GG_pool_no_partitions, GDSM, asc, vg, S, taula, T_next, T, likvid, R, cell, W, random_distances, Selection), [(faza, cas)])
                faza=GGL[0][0]
                cas=GGL[0][1]

                #for attr2 in Selection:
                 #   faza, cas = mpf.GG(GDSM, faza, asc, vg, S, taula, T_next, T, likvid, R, cell, W, cas, random_distances, attr2)

                #faza, cas = mpf.for_loop_GG(GDSM, asc, vg, S, taula, T_next, T, likvid, R, cell, W, random_distances, Selection, [faza, cas])
                #faza, cas = mpf.for_loop_GG_attrs(GDSM, asc, vg, S, taula, T_next, T, likvid, R, cell, W, random_distances, Selection, [faza, cas])

                #faza1=GGL[0][0]
                #cas1=GGL[0][1]

                #faza2=GGL[1][0]
                #cas2=GGL[1][1]


                #faza=GGL[0][1]
                #cas=GGL[0][12]

                #faza=np.dstack((GGL[0][0][:,:,Ycut_indices[0]:Ycut_indices[1]], GGL[1][0][:,:,1:]))   #  Stacking of arrays horizontally (in Y-direction)
                #cas=np.dstack((GGL[0][1][:,:,Ycut_indices[0]:Ycut_indices[1]], GGL[1][1][:,:,1:]))       #  Stacking of arrays horizontally (in Y-direction)


                #faza, cas = mpf.GG_process(GDSM, asc, R, cell, W, random_distances, Selection, partlist[0], None)

                #p1 = Process(target=mpf.GG_process, args= (GDSM, asc, R, cell, W, cas, random_distances, Selection, partlist[0], None)); print('p1 = Process....DONE')
                #p1.start(); print('p1.start()....DONE')
                #p1.join(); print('p1.join()....DONE')
                

                
                #p1 = Process(target=mpf.arb_func, args=(2, Q1)); print('p1 = Process....DONE')
                #p1 = Process(target=mpf.GG_process, args= (GDSM, asc, R, cell, W, cas, random_distances, Selection, partlist[0], Q1)); print('p1 = Process....DONE')
                #p1.start(); print('p1.start()....DONE')
                #faza=Q1.get()[0]  ; print('faza=Q1.get()[0] ....DONE')
                #cas=Q1.get()[1]   ; print('cas=Q1.get()[1] ....DONE')
                #p1.join()
                #out = Q1.get(); print(out)

                '''
                faza1=faza[:,:,Ycut_indices[0]:Ycut_indices[1]+1]                        # 1st partition (in Y-direction)
                faza2=faza[:,:,Ycut_indices[1]-1:Ycut_indices[2]]                         # 2nd partition (in Y-direction)

                cas1=cas[:,:,Ycut_indices[0]:Ycut_indices[1]+1]                        # 1st partition (in Y-direction)
                cas2=cas[:,:,Ycut_indices[1]-1:Ycut_indices[2]]                         # 2nd partition (in Y-direction)

                p1 = Process(target=mpf.GG_process, args= (GDSM, faza1, asc, vg, S, taula, T_next, T, likvid, R, cell, W, cas1, random_distances, Selection, Q1))
                p2 = Process(target=mpf.GG_process, args= (GDSM, faza1, asc, vg, S, taula, T_next, T, likvid, R, cell, W, cas2, random_distances, Selection, Q2))
               
                p1.start()
                p2.start()
            
                faza1=Q1.get()[0]
                faza2=Q2.get()[0]

                cas1=Q1.get()[1]
                cas2=Q2.get()[1]
            
                p1.join()
                p2.join()
            
                faza=np.dstack((faza1[:,:,Ycut_indices[0]:Ycut_indices[1]], faza2[:,:,1:]))   #  Stacking of arrays horizontally (in Y-direction)
                cas=np.dstack((cas1[:,:,Ycut_indices[0]:Ycut_indices[1]], cas2[:,:,1:]))       #  Stacking of arrays horizontally (in Y-direction)
                '''

                #GG_list=p.map(partial( mpf.GG, GDSM, faza, asc, vg, S, taula, T_next, T, likvid, R, cell, W, cas, random_distances), Selection)
                #CPU_RAM()
                # ... GG_list = [(attr1, attr3, attr4, attr5, attr6, attr7, attr8, attr9, attr10, attr11, attr12, attr13, attr14, attr15, attr2)]
                # ... attr1= GDSM, attr2= Selection, attr3= faza, attr4= asc, attr5= vg, attr6= S, attr7= taula, attr8= T_next, attr9= T,
                # ... attr10= likvid, attr11= R, attr12= cell, attr13= W, attr14= cas, attr15= random_distances 

            '''
            elif not avtomatska_nukleacija.get():
                
                GG_list = p.map(partial( mpf.GG_, faza), Selection)   # --------> returns :: (faza, cas)
                for selsel in Selection:
                    faza, cas = mpf.GG(GDSM, faza, asc, vg, S, taula, T_next, T, likvid, R, cell, W, cas, random_distances, selsel)
            
            for gp in range(len(GG_list)-1):
                if gp==0:
                    faza=np.where(GG_list[gp+1]==0, GG_list[gp], GG_list[gp+1])
                else:
                    faza=np.where(GG_list[gp+1]==0, faza, GG_list[gp+1])

            #cas=np.where(S['00_1']==0, S['001'], S['00_1'])
            cas=S['001']
            '''

            if not np.all(faza==F):
                #negind -=1
                negind = min(list(Negatives.keys()))- 1             #  Here, negind should be defined as absolute counter of negative numbers to avoid erasing existing Negatives items
                #for s in S._getvalue():
                for s in S:
                    S[s][np.isin(faza-F, grain_ID)]= negind
                Negatives[negind]=0


    '''
    if GDSM == 'new':
        #*********************************************************************** NEW ***************************************************************************************************
        grain_size = np.count_nonzero(faza == grain)
        asc[grain]['grain size'] = grain_size 
        # ----------------------------------------- this function selects inactive grains (i.e. grains, which no longer grow -
        # ----------------------------------------- either domain limits were reached and/or has reached another grain)
       
        if delete_inactive_grains_ID and len(FF)== FF_length:
            active=np.isin(faza-FF[0], grain)
            if np.all(active==False):
               inactive_grains.append(grain)
       
               #IG[grain].append(i)
            #else:
               #AG[grain].append(i)
        # -----------------------------------------
        #************************************************************************* NEW ***************************************************************************************************
    '''
    if GDSM == 'old':
        #------------------------------------------------------------------------------------  old  ---------------------------------------------------------------------------------------------------------------
        total=merging(seznam_premikov[0], seznam_premikov[1])
        if len(smeri)> 2:
            for move in range(2, len(smeri)):
                total=merging(total, seznam_premikov[move])
        faza = total.copy()
        grain_size = np.count_nonzero(faza == grain)
        asc[grain]['grain size'] = grain_size
        '''
        with open(directory.get()+track+'nuclei_data.json', 'w') as nuks:              # Writing data of ALL nuclei as .json file, but values must be list NOT np.array !!!
            asc_list =asc.copy()
            for nuk in asc:
                asc_list[nuk]['oi']=asc[nuk]['oi'].tolist()                                
                asc_list[nuk]['rgb']=asc[nuk]['rgb'].tolist()
            json.dump(asc_list, nuks)
        ''' 
        # ----------------------------------------- this function selects inactive grains (i.e. grains, which no longer grow - either domain limits were reached and/or has reached another grain)
        if delete_inactive_grains_ID and len(FF)==FF_length:
            active=np.isin(faza-FF[0], grain)
            if np.all(active==False):
               inactive_grains.append(grain)
               #IG[grain].append(i)
            #else:
               #AG[grain].append(i)
        # -----------------------------------------
        #------------------------------------------------------------------------------------  old  ---------------------------------------------------------------------------------------------------------------

    grain_ID = np.array(list(asc.keys()))
    grain_ID_ = grain_ID
    
    if delete_inactive_grains_ID and len(FF)== FF_length:

        '''............................................ erasing inactive grains from grain_ID register ........................................'''
        #for ggg in grain_ID:
         #   active=np.isin(faza-FF[0], ggg)  
          #  if np.all(active==False):
           #     inactive_grains.append(ggg)
                
        grain_ID_ = np.array(list(set(grain_ID)-set(inactive_grains)))    #  array of active grains (grain_ID_)
        
        '''.......................................................................................................................................................................'''
        del FF[0]

        '''............................................ erasing Negatives items with values less than dt_thresh ........................................'''

        if avtomatska_nukleacija.get():
            Negatives ={key:val for key, val in Negatives.items() if val < dt_thresh}
            #Negatives ={key:val for key, val in Negatives.items() if val < 1000*dt}
            
        else:
            if negind<=critical_negind:
                Negatives ={key:val for key, val in Negatives.items() if key < (negind+negatives_thresh)}
    
        '''.......................................................................................................................................................................'''

    else:
        grain_ID_ = grain_ID

    if GDSM == 'old':
        #--------------------------------------------------------------------------------------------- TIME MATRIX part - OLD --------------------<<<
        times = np.array(list(S.values()))
        cas_total = merging(times[0], times[1])
        if len(smeri)> 2:
            for move in range(2, len(smeri)):
                cas_total = merging(cas_total, times[move])
        try:
            cas = cas_total.copy()
        except NameError:
            pass
        #--------------------------------------------------------------------------------------------- TIME MATRIX part - OLD --------------------<<<
    
    rgb_snap = np.zeros((faza.shape[0], faza.shape[1], faza.shape[2], 3)).astype('int')
    for zrno in asc:
        rgb_snap[faza==zrno]=asc[zrno]['rgb']
    
    if save_flashy_as_RGB.get():
        try:
            np.save(directory.get()+flashies_RGB+'flashy_RGB_'+str(i)+'.npy', rgb_snap)
        except FileNotFoundError:
            os.mkdir(directory.get()+flashies_RGB)
            np.save(directory.get()+flashies_RGB+'flashy_RGB_'+str(i)+'.npy', rgb_snap)
  
    if save_flashy_as_faza.get():
        try:
            np.save(directory.get()+flashies_faza+'flashy_faza_'+str(i)+'.npy', faza)
        except FileNotFoundError:
            os.mkdir(directory.get()+flashies_faza)
            np.save(directory.get()+flashies_faza+'flashy_faza_'+str(i)+'.npy', faza)

    step_time_end = time.time()
    step_cpu = shape_double(step_time_end- step_time_start, 3)
    try:
        rate = round((step_time_end- start_time)/h, 4)
    except ZeroDivisionError:
        rate = 'inf.'
        
    print('3rd run time:  ',round(step_time_end-run2_time_end, 4),' sec.')#; print('Selection: ',Selection,'    grain_ID: ',grain_ID); print()
    active_grains= grain_ID_.shape[0]; real_length= shape_double(FEM_scanning_speed*dt*h, 4); real_time=shape_double(1000*dt*h, 4)
    print('TIME STEP # ',i,'    |    current step CPU: ',step_cpu,' s    |    RATE = ',rate,' s / step    |    active grains / ALL grains :  ',str(active_grains),' / ',grain_counter)
    print(34*' ','length = ',real_length,' mm              |      time = ',real_time,' msec.','      |    inactive grains: ',len(list(set(inactive_grains))))
    CPU, RAM = CPU_RAM(); CPU='%.1f' % CPU ;  RAM='%.1f' % RAM
    print(30*' .  . ')


    # Tkinter display output ---------------------------------------------------------------------------
    var0.set(i); var1.set(step_cpu); var2.set(active_grains); var3.set(grain_counter); var4.set(real_length); var5.set(real_time)
    var6.set(CPU); var7.set(RAM)
    #  -------------------------------------------------------------------------------------------------------
    
    PD['h'].append(h)
    PD['current step CPU'].append(step_cpu)
    PD['ALL grains'].append(grain_counter)

    #if not bool(grain_ID.shape[0]) and grain_counter:
     #   root.after_cancel(run_ca)

    if i>=END_step:
        #root.after_cancel(run_ca)
        # ......................................... PRIKAZ in IZPIS REZULTATOV ..................................................#
        Save_Last()
        
        process_time = time.time() - start_time
        print();print(); print()
        print(20*' ',100*'*')
        print(25*' ','DONE!!  Total computing time =  ',round(process_time, 3),'  seconds.   (RATE = ', round(process_time/h, 4),' s / step)'); print(25*' ',100*'-')
        print(25*' ','Printed length real dimension = ',round(FEM_scanning_speed*dt*h, 5),u'\u00B5m      |    Printing time = ',round(1000*dt*h, 5),' msec.')
        print(25*' ',100*'-')
        print(65*' ','Number of grains =  ',grain_counter,' !'); print(20*' ',100*'*'); print(); print(); print()
        
    #elif i<END_step:
     #   run_ca=root.after(1000, lambda:Run_CA(p))

    root.update()

    
run_safety_switch = False
def RUN():
    global run_safety_switch
    if not run_safety_switch:
        run_safety_switch=True
        zavora.set('pass')
        ALL_variables_refresh()

        for i in range(START_step, END_step):
            if zavora.get()=='pass':
                Run_CA(pool)
            else:
                run_safety_switch=False
                break
        
def Stop_CA():
    global run_safety_switch
    try:
        zavora.set('break')
        run_safety_switch=False
        #root.after_cancel(run_ca)
        Save_Last()
    except NameError:
        pass

def Save_Last():
    global cuts_RGB, cut_count, rgb_snap, cutoff_limit
    #plt.imshow(rgb_snap[0]); plt.figure(); plt.imshow(faza[0])
    if save_kickoff.get():
        Save_KickOff()
    if save_last_cut_figure.get():
        np.save(directory.get()+cuts_RGB+'cut_RGB_'+str(cut_count)+'.npy', rgb_snap[:,:,:cutoff_limit,:])
        np.save(directory.get()+cuts_RGB+'cut_RGB_'+str(cut_count+1)+'.npy', rgb_snap[:,:,cutoff_limit:,:])
        np.save(directory.get()+cuts_faza+'cut_faza_'+str(cut_count)+'.npy', faza[:,:,:cutoff_limit])
        np.save(directory.get()+cuts_faza+'cut_faza_'+str(cut_count+1)+'.npy', faza[:,:,cutoff_limit:])

Z_layer = 0
def plot_RGB(z_layer):
    plt.figure('Grain Oriented Microstructure (EBSD)at  Z= '+str(z_layer))
    plt.imshow(rgb_snap[z_layer])
    plt.show()

lab_but0 = Label(root, text=10*' ', width=5, height=2, bg=basic_color);  lab_but0.grid(row=199, column=100, columnspan=1)
but0 = Button(root, text='RUN', command=lambda: RUN(), font=('Arial', 20), fg='green' , width=15); but0.grid(row=200, column=100, columnspan=3, sticky='W')
but1 = Button(root, text='Stop', command=Stop_CA, font=('Arial', 20), fg='red', width=15); but1.grid(row=200, column=101, columnspan=3)
lab_but2 = Label(root, text=1*' ', width=1, height=1, bg=basic_color);  lab_but2.grid(row=100, column=105, columnspan=1)

frame = Frame(root, bg=basic_color); frame.grid(row=100, column=106, sticky='NW')

but2 = Button(frame, text='Plot CPU vs. Time Step', command=plot_CPU_i, font=('Arial', '11'), width=20); but2.grid(row=10, column=11)
#but3 = Button(frame, text='Plot CPU vs. Time Step', command=plot_CPU_i, font=('Arial', '11'), width=20); but3.grid(row=10, column=11)

but4 = Button(frame, text='Show EBSD', command=lambda: plot_RGB(Z_layer), font=('Arial', '11'), width=20); but4.grid(row=12, column=11)


# IMPORTANT ! -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''... uncomment the following line when making executable file .........................'''
#root.mainloop()
''' ---------------------------------------------------------------------------------- the end -------------------------------------------------------------------------------------------------------------- '''






