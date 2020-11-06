import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import sys

"""Author: Tom Howie 25317679"""

"""************************************************************************"""
"""Question 1 (Four Point Stencil) & general PDE setup- Function Defintions"""

"""This section sets up the Partial Differential Equation(PDE 1) in the matrix 
form Ax = b using code provided by Prof Simon J. Cox, Dr Samuel Sinayoko[1].
A = the laplacian operator
x = u -- which consists of the unknown variables we are trying to find
B = rho -- rho(0.5,0.5)=2 = (and zero elsewhere)
Boundary conditions are given: as u = u(x,y), and  where  0≤x,y,≤1, and 
u = 0 on the boundaries 

The size of the mesh used is defined to be of size n with n**2 unknowns
 to be found in it. So the matrix will be of size (n+2) with boundary 
 conditions included
 
 This section also sets up the embed adn check for rho(0.5,0.5)=2 functions 
 for use with Q1,2,4. Also defines the plot_graph function"""
 
n=5 # Defines Size, n must be odd to have centre point defined.
    #And n must be greater than 3 fo 4p or 5 for 8p stencil. 
h = 1/(n+1.0) #interval on  mesh for graph

#if (n % 2 == 0):
#    sys.exit("Mesh size is not odd process stopped") 

def A_mat(n):
    """This function returns the A matrix of the PDE in question when expressed
    in the matrix form Ax=b. It is of the size n**2 x n**2"""
    
    # Clear matrix and set it up
    a=np.zeros([n**2,n**2])
    #Build full matrix Interior"""
    for i in range(1,n-1):
        for j in range(1,n-1):
            north = (i-1)*n+j
            west = i*n+j-1
            index= i*n+j
            east = i*n+j+1
            south = (i+1)*n+j
            
            a[index,north]=1
            a[index,west] =1
            a[index,index]=-4
            a[index,east] =1
            a[index,south]=1
    
    #North/ Top (nothing further North)
    i=0 #First Row Number
    #Note that the range (1,n-1) means that we JUST middle ones
    #e.g. if n=5 then range(1,4) =[1,2,3]
    for j in range(1,n-1):
        #north = (i-1)*n+j
        west = i*n+j-1
        index= i*n+j
        east = i*n+j+1
        south = (i+1)*n+j
    
        #a[index,north]=1
        a[index,west] =1
        a[index,index]=-4
        a[index,east] =1
        a[index,south]=1
        
    
    #West/ Left (nothing further West)
    j=0 #First Column Number
    
    for i in range(1,n-1):
        north = (i-1)*n+j
        #west = i*n+j-1
        index= i*n+j
        east = i*n+j+1
        south = (i+1)*n+j
    
        a[index,north]=1
        #a[index,west] =1
        a[index,index]=-4
        a[index,east] =1
        a[index,south]=1
        
    
    #East/ Right (nothing further East)
    j=n-1 # Last Column number
    for i in range(1,n-1):
        north = (i-1)*n+j
        west = i*n+j-1
        index= i*n+j
        #east = i*n+j+1
        south = (i+1)*n+j
    
        a[index,north]=1
        a[index,west] =1
        a[index,index]=-4
        #a[index,east] =1
        a[index,south]=1
        #print i,j,index
        
    #South/ Bottom (nothing further South)
    i=n-1 # Last row number
    for j in range(1,n-1):
        north = (i-1)*n+j
        west = i*n+j-1
        index= i*n+j
        east = i*n+j+1
        #south = (i+1)*n+j
    
        a[index,north]=1
        a[index,west] =1
        a[index,index]=-4
        a[index,east] =1
        #a[index,south]=1
    
    #Top Left
    i=0
    j=0
    index= i*n+j
    east = i*n+j+1
    south = (i+1)*n+j
    
    a[index,index]=-4
    a[index,east] =1
    a[index,south]=1
    
    #Top Right
    i=0
    j=n-1
    west = i*n+j-1
    index= i*n+j
    south = (i+1)*n+j
    
    a[index,west] =1
    a[index,index]=-4
    a[index,south]=1
        
    #Bottom Left
    i=n-1
    j=0
    north = (i-1)*n+j
    index= i*n+j
    east = i*n+j+1
    
    a[index,north]=1
    a[index,index]=-4
    a[index,east] =1
    
    #Bottom Right
    i=n-1
    j=n-1
    north = (i-1)*n+j
    west = i*n+j-1
    index= i*n+j
    
    a[index,north]=1
    a[index,west] =1
    a[index,index]=-4
                   
    return a

def b_vec(n):
    """This function returns the vector b which represents the values of rho, which
    we know""" 
    b=np.zeros([n**2,1])
    mid_val = (n**2-1)/2.0
    b[mid_val] = 2.0
    
    return b

def A_mat_test(n):   
    """Test to check if A is set up correctly"""
    #Set up the printing of the array so it shows "nicely"
    np.set_printoptions(precision=0,linewidth=140)   
    A_check = A_mat(n)
    print ("")
    print 'A = '
    print (A_check)
    print ("")
    #Reset printing options
    np.set_printoptions(edgeitems=3,infstr='inf',linewidth=75, nanstr='nan', 
                       precision=8,suppress=False, threshold=1000)
    
    return "Test A (4 point) complete"

def b_vec_test(n):
    """Test to check if B is set up correctly"""
    b_check = b_vec(n)
    print ("")
    print 'b = '
    print b_check
    print ("")
    
    return "Test B complete"

def embed(soln):  
    """In order to check whether the code obeys rho(0.5,0.5)=2 the solution 
    must before wrapped(reshaped) into the correct form before finding the 
    value of rho using the stencil. The function wraps the column vector     
    solution into the correct form of a n*n matrix, and then embeds it with 
    the zero boundary conditions to give a n+2*n+2 matrix."""
    U_wrp=np.reshape(soln,[n,n])
    U_emb=np.zeros([n+2,n+2])
    for i in range(1,n+1):
        for j in range(1,n+1):
            U_emb[i,j]=U_wrp[i-1,j-1]
    return U_emb

def rho_check(U):
    """The functionchecks whether the solution value of rho(0.5,0.5)=2 as given
    in Q. Rather than creating a loop to iterate through all values, the
    equation works for the centre point, saving computational time."""
    i = (n+1)/2.0 #array index for midpoint on stencil mesh
    j = i #array index for midpoint on stencil mesh
    rho = (U[i-1,j]+U[i,j-1]+-4*U[i,j]+U[i,j+1]+U[i+1,j])
    return rho

def plot_graph(U,n):
    """Function to plot the graph of the b vector / rho values."""
    fig = plt.figure(3)
    ax = fig.gca(projection='3d')
    X = np.arange(0,n+2,1)*h
    Y = np.arange(0,n+2,1)*h
    X, Y = np.meshgrid(X,Y)
    ax.plot_wireframe(X, Y, U, rstride=1, cstride=1)
    ax.set_xlabel('x position')
    ax.set_ylabel('y position')
    ax.set_zlabel('u(x,y)')
    plt.show()

"""************************************************************************"""
"""Question 2 - Over-Relaxation Gauss-Seidel - Function Defintions"""

"""The solution is found using an over-relaxation technique, the Gauss-Seidel 
method. The function iter_eq defines the equations that are to be used by the 
function gauss_seide, in order to find the solutions to a given tolerance"""

def iter_eq(A,x,b,omega):
    """This function provides the equations for iteration in order to determine
    the x values of Ax=b. Omega is the relaxation parameter."""
    for i in xrange(n**2):
        t=0.0
        for j in xrange(n**2): 
            if j!=i:          
                t=t+(A[i,j]*x[j])
        x[i]=(omega/A[i,i])*(b[i]-t)+(1-omega)*x[i]
    return x

def gauss_seidel(A,x,b,tol = 1.0e-7):
    """Similar to that of J. Kiusalaas' python code [2] the function uses the
    Gauss-Seidel Method to use the iteration equations as defined by function 
    inter_eq in order to find the solutions of x to within the given tolerance. 
    The function compares the previous iteration to the current iteration
    to check if they are within the given tolerance. If they are, the current 
    iteration is the result/solution of the function, this result, the 
    number of iterations and value of omega is returned. However if after 1000 
    iterations the tolerance is not met, the function prints an error message
    and returns None"""
   
    omega = 1.0
    k = 10
    p = 1
    for i in xrange(1,1001):
        xOld = x.copy()
        x = iter_eq(A,x,b,omega) 
        dx = math.sqrt(np.dot(x-xOld,x-xOld))
        if dx < tol: 
            return x,i,omega
        # Compute of relaxation factor after k+p iterations
        if i == k: dx1 = dx
        if i == k + p:
            dx2 = dx
            omega = 2.0/(1.0 + math.sqrt(1.0 - (dx2/dx1)**(1.0/p)))
    print('Gauss-Seidel failed to converge')
    return None
    
"""************************************************************************"""
"""Question 3 - 8 Point Stencil - Function Defintions"""


def A_8p_mat(n):
    """This function returns the A matrix of the PDE in question when expressed
    in the matrix form Ax=b. It is of the size n**2 x n**2 and for use with the
    8 point stencil"""
    
    # Clear matrix and set it up
    a=np.zeros([n**2,n**2])
    #Build full matrix Interior"""
    for i in range(2,n-2):
        for j in range(2,n-2):
            N2 = (i-2)*n+j
            N1 = (i-1)*n+j
            W2 = i*n+j-2
            W1 = i*n+j-1
            Index= i*n+j
            E1 = i*n+j+1
            E2 = i*n+j+2
            S1 = (i+1)*n+j
            S2 = (i+2)*n+j
            a[Index,N2]=-1
            a[Index,N1]=16
            a[Index,W2]=-1
            a[Index,W1] =16
            a[Index,Index]=-60
            a[Index,E1] =16
            a[Index,E2]=-1
            a[Index,S1]=16
            a[Index,S2]=-1
   
    #North/Top inner edge N1
    i=1 
    for j in range(2,n-2):
        #N2 = (i-2)*n+j
        N1 = (i-1)*n+j
        W2 = i*n+j-2
        W1 = i*n+j-1
        Index= i*n+j
        E1 = i*n+j+1
        E2 = i*n+j+2
        S1 = (i+1)*n+j
        S2 = (i+2)*n+j
        #a[Index,N2]=-1
        a[Index,N1]=16
        a[Index,W2]=-1
        a[Index,W1] =16
        a[Index,Index]=-60
        a[Index,E1] =16
        a[Index,E2]=-1
        a[Index,S1]=16
        a[Index,S2]=-1
    
    #North/Top outer edge N2
    i=0
    for j in range(2,n-2):
        #N2 = (i-2)*n+j
        #N1 = (i-1)*n+j
        W2 = i*n+j-2
        W1 = i*n+j-1
        Index= i*n+j
        E1 = i*n+j+1
        E2 = i*n+j+2
        S1 = (i+1)*n+j
        S2 = (i+2)*n+j
        #a[Index,N2]=-1
        #a[Index,N1]=16
        a[Index,W2]=-1
        a[Index,W1] =16
        a[Index,Index]=-60
        a[Index,E1] =16
        a[Index,E2]=-1
        a[Index,S1]=16
        a[Index,S2]=-1
    
    #West/Left inner edge W1
    j=1 
    for i in range(2,n-2):
        N2 = (i-2)*n+j
        N1 = (i-1)*n+j
        #W2 = i*n+j-2
        W1 = i*n+j-1
        Index= i*n+j
        E1 = i*n+j+1
        E2 = i*n+j+2
        S1 = (i+1)*n+j
        S2 = (i+2)*n+j
        a[Index,N2]=-1
        a[Index,N1]=16
        #a[Index,W2]=-1
        a[Index,W1] =16
        a[Index,Index]=-60
        a[Index,E1] =16
        a[Index,E2]=-1
        a[Index,S1]=16
        a[Index,S2]=-1
        
    #West/Left outer edge W2
    j=0
    for i in range(2,n-2):
        N2 = (i-2)*n+j
        N1 = (i-1)*n+j
        #W2 = i*n+j-2
        #W1 = i*n+j-1
        Index= i*n+j
        E1 = i*n+j+1
        E2 = i*n+j+2
        S1 = (i+1)*n+j
        S2 = (i+2)*n+j  
        a[Index,N2]=-1
        a[Index,N1]=16
        #a[Index,W2]=-1
        #a[Index,W1] =16
        a[Index,Index]=-60
        a[Index,E1] =16
        a[Index,E2]=-1
        a[Index,S1]=16
        a[Index,S2]=-1
    
    #East/Right inner E1
    j=n-2
    for i in range(2,n-2):
        N2 = (i-2)*n+j
        N1 = (i-1)*n+j
        W2 = i*n+j-2
        W1 = i*n+j-1
        Index= i*n+j
        E1 = i*n+j+1
        #E2 = i*n+j+2
        S1 = (i+1)*n+j
        S2 = (i+2)*n+j
        a[Index,N2]=-1
        a[Index,N1]=16
        a[Index,W2]=-1
        a[Index,W1] =16
        a[Index,Index]=-60
        a[Index,E1] =16
        #a[Index,E2]=-1
        a[Index,S1]=16
        a[Index,S2]=-1
        
    #E1ast/Right outer edge E2
    j=n-1
    for i in range(2,n-2):
        N2 = (i-2)*n+j
        N1 = (i-1)*n+j
        W2 = i*n+j-2
        W1 = i*n+j-1
        Index= i*n+j
        #E1 = i*n+j+1
        #E2 = i*n+j+2
        S1 = (i+1)*n+j
        S2 = (i+2)*n+j
        a[Index,N2]=-1
        a[Index,N1]=16
        a[Index,W2]=-1
        a[Index,W1] =16
        a[Index,Index]=-60
        #a[Index,E1] =16
        #a[Index,E2]=-1
        a[Index,S1]=16
        a[Index,S2]=-1
    
    #South/Bottom inner edge S1
    i=n-2
    for j in range(2,n-2):
        N2 = (i-2)*n+j
        N1 = (i-1)*n+j
        W2 = i*n+j-2
        W1 = i*n+j-1
        Index= i*n+j
        E1 = i*n+j+1
        E2 = i*n+j+2
        S1 = (i+1)*n+j
        #S2 = (i+2)*n+j
        a[Index,N2]=-1
        a[Index,N1]=16
        a[Index,W2]=-1
        a[Index,W1] =16
        a[Index,Index]=-60
        a[Index,E1] =16
        a[Index,E2]=-1
        a[Index,S1]=16
        #a[Index,S2]=-1
    
    #South/Bottom Outer edge S2
    i=n-1
    for j in range(2,n-2):
        N2 = (i-2)*n+j
        N1 = (i-1)*n+j
        W2 = i*n+j-2
        W1 = i*n+j-1
        Index= i*n+j
        E1 = i*n+j+1
        E2 = i*n+j+2
        #S1 = (i+1)*n+j
        #S2 = (i+2)*n+j
        a[Index,N2]=-1
        a[Index,N1]=16
        a[Index,W2]=-1
        a[Index,W1] =16
        a[Index,Index]=-60
        a[Index,E1] =16
        a[Index,E2]=-1
        #a[Index,S1]=16
        #a[Index,S2]=-1
    
        
    #N1W1
    i=1
    j=1
    #N2 = (i-2)*n+j
    N1 = (i-1)*n+j
    #W2 = i*n+j-2
    W1 = i*n+j-1
    Index= i*n+j
    E1 = i*n+j+1
    E2 = i*n+j+2
    S1 = (i+1)*n+j
    S2 = (i+2)*n+j
    #a[Index,N2]=-1
    a[Index,N1]=16
    #a[Index,W2]=-1
    a[Index,W1] =16
    a[Index,Index]=-60
    a[Index,E1] =16
    a[Index,E2]=-1
    a[Index,S1]=16
    a[Index,S2]=-1
    
    #N2W1
    i=0
    j=1
    #N2 = (i-2)*n+j
    #N1 = (i-1)*n+j
    #W2 = i*n+j-2
    W1 = i*n+j-1
    Index= i*n+j
    E1 = i*n+j+1
    E2 = i*n+j+2
    S1 = (i+1)*n+j
    S2 = (i+2)*n+j
    #a[Index,N2]=-1
    #a[Index,N1]=16
    #a[Index,W2]=-1
    a[Index,W1] =16
    a[Index,Index]=-60
    a[Index,E1] =16
    a[Index,E2]=-1
    a[Index,S1]=16
    a[Index,S2]=-1
    
    #N1W2
    i=1
    j=0
    #N2 = (i-2)*n+j
    N1 = (i-1)*n+j
    #W2 = i*n+j-2
    #W1 = i*n+j-1
    Index= i*n+j
    E1 = i*n+j+1
    E2 = i*n+j+2
    S1 = (i+1)*n+j
    S2 = (i+2)*n+j
    #a[Index,N2]=-1
    a[Index,N1]=16
    #a[Index,W2]=-1
    #a[Index,W1] =16
    a[Index,Index]=-60
    a[Index,E1] =16
    a[Index,E2]=-1
    a[Index,S1]=16
    a[Index,S2]=-1  
    
    #N2W2
    i=0
    j=0
    #N2 = (i-2)*n+j
    #N1 = (i-1)*n+j
    #W2 = i*n+j-2
    #W1 = i*n+j-1
    Index= i*n+j
    E1 = i*n+j+1
    E2 = i*n+j+2
    S1 = (i+1)*n+j
    S2 = (i+2)*n+j
    #a[Index,N2]=-1
    #a[Index,N1]=16
    #a[Index,W2]=-1
    #a[Index,W1] =16
    a[Index,Index]=-60
    a[Index,E1] =16
    a[Index,E2]=-1
    a[Index,S1]=16
    a[Index,S2]=-1
    
    #N1E1
    i=1
    j=n-2
    #N2 = (i-2)*n+j
    N1 = (i-1)*n+j
    W2 = i*n+j-2
    W1 = i*n+j-1
    Index= i*n+j
    E1 = i*n+j+1
    #E2 = i*n+j+2
    S1 = (i+1)*n+j
    S2 = (i+2)*n+j
    #a[Index,N2]=-1
    a[Index,N1]=16
    a[Index,W2]=-1
    a[Index,W1] =16
    a[Index,Index]=-60
    a[Index,E1] =16
    #a[Index,E2]=-1
    a[Index,S1]=16
    a[Index,S2]=-1
    
    #N2E1
    i=0
    j=n-2
    #N2 = (i-2)*n+j
    #N1 = (i-1)*n+j
    W2 = i*n+j-2
    W1 = i*n+j-1
    Index= i*n+j
    E1 = i*n+j+1
    #E2 = i*n+j+2
    S1 = (i+1)*n+j
    S2 = (i+2)*n+j
    #a[Index,N2]=-1
    #a[Index,N1]=16
    a[Index,W2]=-1
    a[Index,W1] =16
    a[Index,Index]=-60
    a[Index,E1] =16
    #a[Index,E2]=-1
    a[Index,S1]=16
    a[Index,S2]=-1
    
    #N1E2
    i=1
    j=n-1
    #N2 = (i-2)*n+j
    N1 = (i-1)*n+j
    W2 = i*n+j-2
    W1 = i*n+j-1
    Index= i*n+j
    #E1 = i*n+j+1
    #E2 = i*n+j+2
    S1 = (i+1)*n+j
    S2 = (i+2)*n+j
    #a[Index,N2]=-1
    a[Index,N1]=16
    a[Index,W2]=-1
    a[Index,W1] =16
    a[Index,Index]=-60
    #a[Index,E1] =16
    #a[Index,E2]=-1
    a[Index,S1]=16
    a[Index,S2]=-1
       
    #N2E2
    i=0
    j=n-1
    #N2 = (i-2)*n+j
    #N1 = (i-1)*n+j
    W2 = i*n+j-2
    W1 = i*n+j-1
    Index= i*n+j
    #E1 = i*n+j+1
    #E2 = i*n+j+2
    S1 = (i+1)*n+j
    S2 = (i+2)*n+j
    #a[Index,N2]=-1
    #a[Index,N1]=16
    a[Index,W2]=-1
    a[Index,W1] =16
    a[Index,Index]=-60
    #a[Index,E1] =16
    #a[Index,E2]=-1
    a[Index,S1]=16
    a[Index,S2]=-1
    
    #S1W1
    i = n-2
    j = 1
    N2 = (i-2)*n+j
    N1 = (i-1)*n+j
    #W2 = i*n+j-2
    W1 = i*n+j-1
    Index= i*n+j
    E1 = i*n+j+1
    E2 = i*n+j+2
    S1 = (i+1)*n+j
    #S2 = (i+2)*n+j
    a[Index,N2]=-1
    a[Index,N1]=16
    #a[Index,W2]=-1
    a[Index,W1] =16
    a[Index,Index]=-60
    a[Index,E1] =16
    a[Index,E2]=-1
    a[Index,S1]=16
    #a[Index,S2]=-1
        
    #S2W1
    i=n-1
    j=1
    N2 = (i-2)*n+j
    N1 = (i-1)*n+j
    #W2 = i*n+j-2
    W1 = i*n+j-1
    Index= i*n+j
    E1 = i*n+j+1
    E2 = i*n+j+2
    #S1 = (i+1)*n+j
    #S2 = (i+2)*n+j
    a[Index,N2]=-1
    a[Index,N1]=16
    #a[Index,W2]=-1
    a[Index,W1] =16
    a[Index,Index]=-60
    a[Index,E1] =16
    a[Index,E2]=-1
    #a[Index,S1]=16
    #a[Index,S2]=-1
    
    #S1W2
    i=n-2
    j=0
    N2 = (i-2)*n+j
    N1 = (i-1)*n+j
    #W2 = i*n+j-2
    #W1 = i*n+j-1
    Index= i*n+j
    E1 = i*n+j+1
    E2 = i*n+j+2
    S1 = (i+1)*n+j
    #S2 = (i+2)*n+j
    
    a[Index,N2]=-1
    a[Index,N1]=16
    #a[Index,W2]=-1
    #a[Index,W1] =16
    a[Index,Index]=-60
    a[Index,E1] =16
    a[Index,E2]=-1
    a[Index,S1]=16
    #a[Index,S2]=-1
  
    #S12W12
    i=n-1
    j=0
    N2 = (i-2)*n+j
    N1 = (i-1)*n+j
    #W2 = i*n+j-2
    #W1 = i*n+j-1
    Index= i*n+j
    E1 = i*n+j+1
    E2 = i*n+j+2
    #S1 = (i+1)*n+j
    #S2 = (i+2)*n+j
    a[Index,N2]=-1
    a[Index,N1]=16
    #a[Index,W2]=-1
    #a[Index,W1] =16
    a[Index,Index]=-60
    a[Index,E1] =16
    a[Index,E2]=-1
    #a[Index,S1]=16
    #a[Index,S2]=-1
    
    #S11E11
    i=n-2
    j=n-2
    N2 = (i-2)*n+j
    N1 = (i-1)*n+j
    W2 = i*n+j-2
    W1 = i*n+j-1
    Index= i*n+j
    E1 = i*n+j+1
    #E2 = i*n+j+2
    S1 = (i+1)*n+j
    #S2 = (i+2)*n+j
    a[Index,N2]=-1
    a[Index,N1]=16
    a[Index,W2]=-1
    a[Index,W1] =16
    a[Index,Index]=-60
    a[Index,E1] =16
    #a[Index,E2]=-1
    a[Index,S1]=16
    #a[Index,S2]=-1
    
    #S11E12
    i=n-2
    j=n-1
    N2 = (i-2)*n+j
    N1 = (i-1)*n+j
    W2 = i*n+j-2
    W1 = i*n+j-1
    Index= i*n+j
    #E1 = i*n+j+1
    #E2 = i*n+j+2
    S1 = (i+1)*n+j
    #S2 = (i+2)*n+j
    a[Index,N2]=-1
    a[Index,N1]=16
    a[Index,W2]=-1
    a[Index,W1] =16
    a[Index,Index]=-60
    #a[Index,E1] =16
    #a[Index,E2]=-1
    a[Index,S1]=16
    #a[Index,S2]=-1
    
    #S12E11
    i=n-1
    j=n-2
    N2 = (i-2)*n+j
    N1 = (i-1)*n+j
    W2 = i*n+j-2
    W1 = i*n+j-1
    Index= i*n+j
    E1 = i*n+j+1
    #E2 = i*n+j+2
    #S1 = (i+1)*n+j
    #S2 = (i+2)*n+j
    a[Index,N2]=-1
    a[Index,N1]=16
    a[Index,W2]=-1
    a[Index,W1] =16
    a[Index,Index]=-60
    a[Index,E1] =16
    #a[Index,E2]=-1
    #a[Index,S1]=16
    #a[Index,S2]=-1
        
    #E12N12
    i=n-1
    j=n-1
    N2 = (i-2)*n+j
    N1 = (i-1)*n+j
    W2 = i*n+j-2
    W1 = i*n+j-1
    Index= i*n+j
    #E1 = i*n+j+1
    #E2 = i*n+j+2
    #S1 = (i+1)*n+j
    #S2 = (i+2)*n+j
    a[Index,N2]=-1
    a[Index,N1]=16
    a[Index,W2]=-1
    a[Index,W1] =16
    a[Index,Index]=-60
    #a[Index,E1] =16
    #a[Index,E2]=-1
    #a[Index,S1]=16
    #a[Index,S2]=-1
    
    return a
    
def b_8p_vec(n):
    """This function returns the vector b * by a factor of 12""" 
    b=np.zeros([n**2,1])
    mid_val = (n**2-1)/2.0
    b[mid_val] = 24.0
    
    return b
    
def A_mat_8p_test(n):   
    """Test to check if A_8p is set up correctly"""
    np.set_printoptions(precision=0,linewidth=120)
    A_8p_check = A_8p_mat(n)
    print ("")
    print "A_8p = "
    print A_8p_check
    print ("")
    np.set_printoptions(precision=0,linewidth=140)
    np.set_printoptions(edgeitems=3,infstr='inf',linewidth=75, nanstr='nan', 
                       precision=8,suppress=False, threshold=1000)
                   
    return "Test A (8 point) complete"
    
def embed_8p(soln):  
    """In order to check whether the code obeys rho(0.5,0.5)=2 the solution 
    must before wrapped(reshaped) into the correct form before finding the 
    value of rho using the stencil. The function wraps the colum vector     
    solution into the correct form of a n*n matrix, and then embeds it with 
    the two layers of zero boundary conditions on each side to give a n+4*n+4 
    matrix."""
    U_wrp=np.reshape(soln,[n,n])
    U_emb=np.zeros([n+4,n+4])
    for i in range(2,n+2):
        for j in range(2,n+2):
            U_emb[i,j]=U_wrp[i-2,j-2]
    return U_emb

def rho_check_8p(U):
    """checks whether the solution value of rho(0.5,0.5)=2 for the 8 point
    stencil as given in Q"""
    i = (n+3)/2 #as there is a double layer of zeros on each side mid index 
    j = i       #value needs to be changed.
    rho = (-U[i,j-2]+16*U[i,j-1]-U[i-2,j]+16*U[i-1,j]-60*U[i,j]+16*U[i+1,j] \
            -U[i+2,j]+16*U[i,j+1]-U[i,j+2])/12.0
    return rho

"""************************************************************************"""
"""Question 4 - Gauss-Seidel “Red-Black” Solver - Function Definitions"""

def red_black_iter_eq(A,x,b,omega):
    """This function is similar to the previous function in q2, providing the 
    equations for iteration, but offers the potential to decrease computational
    time as the loops are seperated, therefore could use parallel computing """
    for i in xrange(1,n**2-1,2):
        t=0.0
        for j in xrange(n**2):
            if j!=i:          
                t=t+(A[i,j]*x[j])
        x[i]=(omega/A[i,i])*(b[i]-t)+(1-omega)*x[i]
         
    for i in xrange(0,n**2,2):
        t=0.0
        for j in xrange(n**2):
            if j!=i:          
                t=t+(A[i,j]*x[j])
        x[i]=(omega/A[i,i])*(b[i]-t)+(1-omega)*x[i]  
    return x
    
def red_black_gauss_seidel(A,x,b,tol = 1.0e-7):
    """Similar to that of J. Kiusalaas' python code [2] the function uses the
    Gauss-Seidel Method to use the iteration equations as defined by function 
    inter_eq in order to find the solutions of x to within the given tolerance. 
    However now it utilises the red_black_iter_eq function."""
   
    omega = 1.0
    k = 10
    p = 1
    for i in xrange(1,1001):
        xOld = x.copy()
        x = red_black_iter_eq(A,x,b,omega) 
        dx = math.sqrt(np.dot(x-xOld,x-xOld))
        if dx < tol: 
            return x,i,omega
        # Compute of relaxation factor after k+p iterations
        if i == k: dx1 = dx
        if i == k + p:
            dx2 = dx
            omega = 2.0/(1.0 + math.sqrt(1.0 - (dx2/dx1)**(1.0/p)))
    print('Gauss-Seidel failed to converge')
    return None
    
"""************************************************************************"""
"""Tests to check that the Matrix A (for both four and five point stencil) and
Vector b has set up correctly."""
#print A_mat_test(n)
#print b_vec_test(n)
#print A_mat_8p_test(n)

"""************************************************************************"""
"""Question 1 (Four Point Stencil) - numpy Linalg Solver - Implementation"""
 
A1=A_mat(n) 
b1=b_vec(n)

#starttime=time.time() 
soln1 = np.linalg.solve(A1,b1)
#T1 = time.time()-starttime

#print("")
#print("Solution 1 =")
#print(soln1)
#print("")

U1 = embed(soln1)
#print("U1=")
#print(U1)
#print("")

rho1 = rho_check(U1)
#plot_graph(U1, n)

print ("Q1 solver (4 Point Stencil - np.linalg): rho(0.5,0.5) = %f" %rho1)
#print"Computational Time = "+ str(T1) + "s"
#print("")

"""************************************************************************"""
"""Question 2 - Over-Relaxation Gauss-Seidel - Implementation"""

A2=A_mat(n) 
b2=b_vec(n) 
soln2=np.zeros([n**2])

#starttime=time.time()
soln2,n_iter,omega = gauss_seidel(A2,soln2,b2)
#T2 = time.time()-starttime
#print("Solution 2 =")
#print(soln2)
#print("")

U2= embed(soln2)
#print("U2 =")
#print(U2)
#print("")

rho2 = rho_check(U2)
#plot_graph(U2, n)

print ("Q2 solver (Over-Relaxation Gauss-Seidel): rho(0.5,0.5) = %f" %rho2)
#print"Computational Time = "+ str(T2) + "s"
#print"Relaxation Factor =",omega
#print"Number of Iterations =",n_iter
#print("")

"""************************************************************************"""
"""Question 3 - 8 Point Stencil - Implementation"""

A3=A_8p_mat(n) 
b3=b_8p_vec(n) 

#starttime=time.time()
soln3 = np.linalg.solve(A3,b3)
#T3 = time.time()-starttime

#print("Solution 3 =")
#print(soln3)
#print("")

U3 = embed_8p(soln3)
#print("U3 =")
#print(U3)
#print("")

rho3 = rho_check_8p(U3)
#plot_graph(U3, n+2)

print("Q3 solver (8 point stencil - np.linalg) rho(0.5,0.5): %f" %rho3)
#print"Computational Time = "+ str(T3) + "s"
#print("")

"""************************************************************************"""
"""Question 4 - Gauss-Seidel “Red-Black” Solver - Function Implementation"""

A4=A_mat(n) 
b4=b_vec(n) 
soln4=np.zeros([n**2])

#starttime=time.time()
soln4,n_iter,omega = red_black_gauss_seidel(A4,soln4,b4)
#T4 = time.time()-starttime

#print("Solution 4 =")
#print(soln4)
#print("")

U4 = embed(soln4)
#print("U4 =")
#print(U4)
#print("")

rho2 = rho_check(U4)
#plot_graph(U4, n)

print("Q4 solver (Red-Black Guass-Seidel): rho(0.5,0.5) = %f" %rho2)
#print"Computational Time "+ str(T4) + "s"
#print"Relaxation Factor =",omega
#print"Number of Iterations =",n_iter
#print("")