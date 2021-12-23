import numpy as np
import random

def random_selection(x, erase):
    global item
    if len(x)!=0:
        item = random.choice(x)
        if erase:
            x=np.delete(x,np.where(x==item))
    else:
        pass
            
    return x


smeri = np.array(['001', '00-1', '010', '0-10', '100', '-100',])
grain_ID =      np.array([1,2,])
GD = {}
for gid in grain_ID:
    GD[gid] = smeri.copy()


selekcija=[]


while True:
    if any(grain_ID):
        grain_ID = random_selection(grain_ID, False)
        grain = item
        counter = 0
    else:
        break
    while counter==0:
        if any(GD[grain]):
            GD[grain] = random_selection(GD[grain], True)
            s = item
            print(20*'-')
            print('zrno: ',grain)
            print('smer: ',s)
            selekcija.append((grain,s))
            counter+=1
            break
        else:
            grain_ID=grain_ID[grain_ID!=grain]
            break
    
    continue
    


