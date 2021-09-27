# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 22:31:07 2021

@author: joyce
"""


from PIL import Image as im
import numpy as np

# Open the original image, the grafiti image and the mask image. 
Im1=im.open('Original.jpg')
original=Im1.convert('L')
Im2=im.open('Graffiti.jpg')
graffiti=Im2.convert('L')
# The mask is bigger than the actual grafiti to avoid problems at the boundary (Isabel's findings)
Im3=im.open('Mask.jpg')
mask=Im3.convert('L')

# Converting images to arrays
o_array = np.array(original)
g_array = np.array(graffiti)
m_array = np.array(mask)

shape = g_array.shape
n = int(shape[0])
m = int(shape[1])

# Initializing empty arrays with the shape of the picture
boundary_array = np.zeros((n,m))
inside_array = np.zeros((n,m))
u0 = np.zeros((n,m))

c = 230
shape = g_array.shape
n = shape[0]
m = shape[1]
for i in range (1, (n-1)):
    for j in range (1, (m-1)):
        # make matrix where the elements corresponding to inside have value 1
        if m_array[i,j]>c:
            m_array[i,j]=255
            inside_array[i,j]=1
        else:
            # make matrix where the elements corresponding to boundaries have value 2
            #m_array[i,j] = 0
            if m_array[i-1][j]>c or m_array[i+1][j]>c or m_array[i][j-1]>c or m_array[i][j+1]>c:
                inside_array[i,j]= 2

for i in range (0,n):
    for j in range(0,m):
        # Fill u only with the values of coordinates corresponding to the area to be restored 
        if inside_array[i,j] == 1:
            u0[i,j] = float(g_array[i,j])
        if inside_array[i,j] == 2:
            u0[i,j] = float(g_array[i,j])
    
count = 50
u = u0                      
while count > 0:
    for i in range(0,n):
        for j in range(0,m):
            if inside_array[i,j] == 1 :
                u[i,j] = (1-1.5)*u[i,j] - (1.5/4)*(u[(i+1), j] + u[(i-1), j] + u[i, (j+1)] + u[(i, (j-1))])
    count-=1

rs = g_array
for i in range(0,n):
    for j in range(0,m):
        if inside_array[i,j] ==1:
            rs[i,j]=u[i,j]
restored=im.fromarray(rs)
restored.show()

'''
# Check quality of restoration
original = o_array
restored = rs

# Initiate variables
total = 0
nb_of_elements = 0
index_array = []
diff_array = []
sol_array = []

# Compute variables
for i in range (0,n):
    for j in range (0,m):
        if inside_array[i,j] == 1: 
            total += original[i,j]
            nb_of_elements += 1
            index_array.append(original[i,j])
            diff_squared = (restored[i,j] - original[i,j])**2
            diff_array.append(diff_squared)
            
mean = total/nb_of_elements

for i in range (0, len(index_array)):
    sol_squared = (index_array[i]-mean)**2
    sol_array.append(sol_squared)
    
sigma_squared = (1/(nb_of_elements - 1)) * (sum(sol_array))
chi_squared = ((1/nb_of_elements)*(np.sum(diff_array)))/(sigma_squared)

print("chi squared: " + str(chi_squared))

# Compute Chi^2
'''
            
                