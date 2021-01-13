import math
from decimal import *
from scipy import poly1d
import scipy.special
import numpy as np
from fractions import Fraction as F
import time
"""
CardinalityFC(n,k)

inputs:
n: number of variables
k: the value of the threshold

CNF clauses are special cases of cardinality constraints

return:
the n*1 vector of Fourier coefficients of Card(x1,x2,...,xn)>=k, in other words, x1+x2+...+xk<=n-2k. Ordered by degree, thus ret[0] is the constant.
"""
def binom(n,k):
    if k==0: return 1
    ans = F(1)
    for i in range(n):
        ans = ans * (i+1)
        if i < k:
            ans = ans/(i+1)
        else:
            ans = ans/(n-i)
    return ans

def verify_paserval(x):
    ans = 0
    for i in range(len(x)):
        ans += x[i]*x[i]*binom(len(x)-1,i)
    return ans

def CardinalityFC(n,k):
#    print('n='+repr(n)+' k='+repr(k))
    startall = time.time()
    if n==1 and k==1:
        return [0,1]
    d = np.zeros(n + 1)
    if k==0:
        d[0]=-1
        return d
    templist = ([-1 for i in range(n-k)] + [1 for i in range(k-1)])
    templist = np.array(templist,dtype=object)
#    print('before poly1d takes '+repr(time.time()-startall))
    
    start = time.time()
    temppoly = poly1d(templist,True) * (1-2*((k-1)%2))   # complexity: O(n(logn)^2)
#    print('poly1d takes '+repr(time.time()-start))
    temparray = temppoly.c[::-1]
    start = time.time()
    temparray = [F(temparray[i],(2**(n-1))) for i in range(len(temparray))]
#    print('temparray/2^n-1 takes '+repr(time.time()-start))
    start = time.time()
    binom_nk = binom(n-1,k-1)
    binom_ni_inv = [1 for i in range(n)]  # iteratively computing binom(n-1,i)
    for i in range(0,n-1):
        binom_ni_inv[i+1] = binom_ni_inv[i] * (i+1) / (n-i-1)
    binocoef = [binom_nk * binom_ni_inv[i] for i in range(n)]  #complexity: O(nk)
    #binocoef = [binom_nk / binom(n-1,i) for i in range(n)]  #complexity: O(nk)
 #   print('binom takes '+repr(time.time()-start))
    temparray = [temparray[i] * binocoef[i] for i in range(len(temparray))]
    for j in range(n):
        d[j+1] = F(temparray[j])
    start = time.time()
    templist3 = []
    templist3.append(binom(n,k))
    for i in range(k,n):
        templist3.append(templist3[i-k] * (n-i) / (i+1) )
  #  print('m1 uses '+repr(time.time()- start))
    #templist3 = [F(binom(n,i)) for i in range(k,n+1)]
    start = time.time()
    d[0] = F(1 - F(int(sum(templist3)),(2 ** (n - 1))))
  #  print('m2 uses '+repr(time.time()- start))
  #  print('CardFC takes '+repr(time.time()-startall))
    return d



# the Fourier Coefficients of Not-all-equal constraints
def NAEFC(n):
    d = np.zeros(n + 1)
    d[0] = (F(1,2))**(n-2)-1
    for i in range(2,n+1,2):
        d[i] = (F(1/2))**(n-2)
    return d

#total = 0
#for i in range(0,101):
#    total = abs(coef[i])*binom(100,i)
#    print(repr(total)+','),


