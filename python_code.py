# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 16:27:16 2021

@author: joyce
"""


from PIL import Image
from numpy import *
import sys
#set_printoptions(threshold=sys.maxsize)

# Change location to image location (if the image is not in the same folder as the python file)!
img = Image.open('rainbow_circle.jpg')
imgGray = img.convert('L')
imgGray.save('grey_circle.jpg')
img_grey = Image.open('grey_circle.jpg')
#img_grey.show()

# Change location to image location
img = Image.open('grey_circle.jpg') 
img_array = array(img)
#print(img_array)
print(img_array.shape)

'''
n = 2
b = img_array.shape[0]//n
img_reshape = img_array.reshape(-1, n, b, n).sum((-1, -3)) / n 
#print(img_reshape)
print(img_reshape.shape)
'''

img_small = Image.fromarray(img_array)
img_small = img_small.convert("L")
img_small.save('small.png')
img_small.show()

'''
# Save matrix as txt
img_mat = matrix(img_array)
with open('outfile.txt','wb') as f:
    for line in img_mat:
        savetxt(f, line, fmt='%.2f')
'''