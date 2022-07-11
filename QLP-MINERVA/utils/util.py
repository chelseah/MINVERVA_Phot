import numpy as np
import scipy as sp
import math
from scipy import linalg

def ls_polynxy(x,y,f,n):
    xi = sp.ones(len(x))
    A=(xi)[:,np.newaxis]

    for i in range(1,n+1):
      for j in range(0,i+1):
        A=np.hstack((A,(1./math.factorial(i-j)*1./math.factorial(j)*x**(i-j)*y**j)[:,np.newaxis]))
    c,resid,rank,sigma = sp.linalg.lstsq(A,f)
    return c

def poly_ordern_eval(x,y,m):
    order = int(np.sqrt(2*(len(m)-1)))-1
    #print len(m),order
    z=np.zeros(len(x))
    k=1
    for i in range(1,order+1):
      for j in range(0,i+1):
        #z+=m[k]*math.factorial(i-j)*math.factorial(j)*x**(i-j)*y**j
        z+=m[k]*1./math.factorial(i-j)*1./math.factorial(j)*x**(i-j)*y**j
        k+=1
    z+=m[0]
    return z


