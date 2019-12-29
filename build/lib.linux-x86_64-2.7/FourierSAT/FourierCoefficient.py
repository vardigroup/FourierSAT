import math
from decimal import *
from scipy import poly1d
import scipy.special
import numpy as np
"""
CardinalityFC(n,k)

inputs:
n: number of variables
k: the value of the threshold

CNF clauses are special cases of cardinality constraints

return:
the n*1 vector of Fourier coefficients of Card(x1,x2,...,xn)>=k, in other words, x1+x2+...+xk<=n-2k. Ordered by degree, thus ret[0] is the constant.
"""
def CardinalityFC(n,k):
    getcontext().prec = 50
    d = np.zeros(n + 1)
    if k==0:
        d[0]=-1
        return d
    templist = ([-1 for i in range(n-k)] + [1 for i in range(k-1)])
    templist = np.array(templist,dtype=object)
    temppoly = poly1d(templist,True) * (1-2.0*((k-1)%2))
    temparray = temppoly.c[::-1]/(2.0**(n-1))
    binocoef = [scipy.special.binom(n-1,k-1)/scipy.special.binom(n-1,i) for i in range(n)]
    temparray = (temparray * binocoef)
    for j in range(n):
        d[j+1] = temparray[j]

    templist3 = [scipy.special.binom(n,i) for i in range(k,n+1)]
    d[0] = 1 - sum(templist3) / (2.0 ** (n - 1))
    return d

# the Fourier Coefficients of Not-all-equal constraints
def NAEFC(n):
    d = np.zeros(n + 1)
    d[0] = (0.5)**(n-2)-1
    for i in range(2,n+1,2):
        d[i] = (0.5)**(n-2)
    return d

