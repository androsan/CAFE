import numpy as np

seznam_premikov = ['left', 'right', 'down', 'up']

init = np.zeros((3,3)); init[1,1]=1

left = init.copy()   ; left[1,0]=1
right = init.copy() ; right[1,2]=1
down = init.copy(); down[2,1]=1
up = init.copy()    ; up[0,1]=1


def merging_meshes(prva,druga):
    tot = np.where(druga==0, prva, druga)
    return tot

right=merging_meshes(left,right)
down=merging_meshes(right,down)
up=merging_meshes(down,up)

print(up)


