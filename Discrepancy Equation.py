import numpy as np
#TO DO: Replace A, B and M with the picture matrixes when adding this part to the main code.

#Matrix of the picture without a mask applied
A=np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
#Matrix of the picture after it has been repaired
B=np.arange(9).reshape(3,3)
#Matrix of just the mask 
M=[[1, 0, 0], [0, 1, 0], [0, 0, 1]]


#Number of matrix elements
n=A.size
#Horizontal matrix size
h=A.shape[0]
#Vertical matrix size
v=A.shape[1]
#First sum in discrepancy equation
I1=0
#Second sum in discrepancy equation
I2=0

#Discrepancy score program

#This loop calculates the sum in the numerator of the discrepancy equation
for i in range(0,h):
    for j in range(0,v):
        #If statement ensures that only points that had the mask are considered in the equation
        if M[i][j] == 1:
            I1=I1+(A[i][j]-B[i][j])
        
#This loop calculates the sum in the denomenator of the discrepancy equation.        
for i in range(0,h):
    for j in range(0,v):
        #If statement ensures that only points that had the mask are considered in the equation
        if M[i][j] == 1:
            Imean=(A[i][j]+B[i][j])/2
            I2=I2+(B[i][j]-Imean)

#Define the denominator in the discrepancy equation
omega=np.sqrt(1/(n-1)*I2**2)

#Calculate chi value
X=np.sqrt(I1**2/(n*omega**2))

#Print chi squared
print(X**2)
