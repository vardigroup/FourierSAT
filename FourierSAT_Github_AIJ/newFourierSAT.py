import scipy.optimize as so
import numpy as np
from scipy import poly1d
import time
import math
import numpy
import random
from FourierCoefficient import CardinalityFC
from FourierCoefficient import NAEFC
from boolean_formula import Formula
from fractions import Fraction as F
import signal
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool
from functools import partial
import argparse
from bddtest import Bdd,gradient_on_bdd, fval_on_bdd

globalargs = ()
globalx = []

def norm(x):
    return sum([x[i]*x[i] for i in range(len(x))])

def prod(iterable):
    p= 1
    for n in iterable:
        p *= n
    return p

class Opt_result:
    def __init__(self):
        self.x = []
        self.fun = 0

def hamming(x,y):
    dis = 0
    for i in range(len(x)):
        if x[i]!=y[i]:
            dis += 1
    return dis

def distance(x,y):
    dis = 0
    for i in range(len(x)):
        dis += (x[i]-y[i])*(x[i]-y[i])
    return dis


def cardcons(x):
    start = time.time()
#    x = truncate(x)
    roots = np.array([-x[i] for i in range(len(x))],dtype=object)
    #poly3 = np.poly1d(roots,True).coeffs
    xlist = [-x[i] for i in range(len(x))]
    return poly1d(xlist,True).c
    
    end = time.time()
        
    
def rounding(y):
    x = list(y)
    for i in range(len(x)):
        if x[i]>0: x[i]=1
        else: x[i]=-1
    return x
def truncate(x):
    for i in range(len(x)):
        if x[i]>1-1e-2: x[i]=1
        if x[i]<-1+1e-2: x[i]=-1
    return x
def count(x):
    y = rounding(x)
    return (len(y)-sum(y))//2
def modify_weight(weight,formula,x,klist):
    x = rounding(x)
    return weight
def project(x):
    for i in range(len(x)):
        if x[i]>1:
            x[i] = 1
        elif x[i]<-1:
            x[i]=-1
    return x
def dot(x,y):
    return sum([x[i]*y[i] for i in range(len(x))])

def relax(x):
    for i in range(len(x)):
        if x[i]>0.8:
            x[i] = 0.8
        elif x[i]<-0.8:
            x[i]=-0.8
    return x

def LS(clauses,x,klist,weight,ctype):
    n = len(x)
    bestx = list(x)
    _, best_unsat_weight, _ = verify_sol(clauses,x,klist,weight,ctype)
    for i in range(n):
        tempx = list(x)
        tempx[i] = -1*tempx[i]
        _, unsat_weight_new, _ = verify_sol(clauses,tempx,klist,weight,ctype) 
        if unsat_weight_new < best_unsat_weight:
            best_unsat_weight = unsat_weight_new
            bestx = list(tempx)
    return bestx,sum(weight) - best_unsat_weight
        

def compute_FC_table(clauses,klist,ctype,coefs,comparator):
    FC_table = []
    for k in range(len(clauses)):
        if ctype[k]=='c':
            FC_table.append(CardinalityFC(len(clauses[k]),klist[k]))
        elif ctype[k]=='x':
            FC_table.append([0]) #place holder
        elif ctype[k]=='n':
            FC_table.append(NAEFC(len(clauses[k])))
        elif ctype[k] == 'p':
            bdd = Bdd()
            bdd.build(coefs[k],klist[k],len(clauses[k]),comparator[k])
       #     print('bdd size = '+repr(bdd.getSize()))
       #     print('k '+repr(klist[k]))
            FC_table.append(bdd)
    return FC_table

def compute_FC_vector(clauses,klist,ctype,weight,n):
    FC_vector = np.zeros(n)
    for k in range(len(clauses)):
        if ctype[k]=='c':
            coef =  CardinalityFC(len(clauses[k]),klist[k])
            for i in range(1,len(coef)):
                FC_vector[i-1] += coef[i]
        if ctype[k]=='x':
            FC_vector[n-1] += 1
        if ctype[k]=='n':
            coef = NAEFC(len(clauses[k]))
            for i in range(1,len(coef)):
                FC_vector[i-1] += coef[i] * weight[k]
    return FC_vector
 
def fast_fval_evaluate(clause,y,k,ctype):
    x = list(y)
    x = [i/abs(i)*x[abs(i)-1] for i in clause]
    n = len(x)
    nump1 = 0
    numn1 = 0
    for i in range(len(x)):
        if x[i]==1:
            nump1 += 1
        elif x[i]==-1:
            numn1 += 1
        else:
            x[i]=0
    num0 = n - numn1 - nump1
    if ctype == 'c':
        if k>0:
            if numn1>=k: return -1
            if nump1>n-k: return 1
        if k<0:
            if nump1>n-k: return -1
            if numn1>k: return 1
    if ctype == 'x':
        if num0==0:
            return prod(x)
    if ctype == 'n':
        if nump1>0 and numn1>0: return -1
        if num0==0 and (numn1==0 or nump1==0): return 1
    return "unknown"
            

def slow_fval_evaluate(clause,x,k,ctype,FC_table):
    if ctype=='c' or ctype=='n':
        cardcons_vec = cardcons([j/abs(j)*x[abs(j)-1] for j in clause])
        four_coef = FC_table
        return dot(cardcons_vec,four_coef)
    elif ctype == 'x':
        return prod(j/abs(j)*x[abs(j)-1] for j in clause)
    elif ctype == 'p':
        x_slice = [ x[abs(i)-1] * i / abs(i) for i in clause ] # assume no negative literal in PB constraints
        return fval_on_bdd(x_slice,FC_table)
def xor_grad(clause,x):
    for i in range(len(clause)):
        if x[abs(clause[i])-1] == 0:
            x[abs(clause[i])-1] = 1e-2
    total_prod = prod([xi/abs(xi)*x[abs(xi)-1] for xi in clause])
    return [total_prod  / x[abs(clause[i])-1] for i in range(len(clause))]

def fast_gradient(clause,x,k,ctype,i):
    x_pos_at_i = list(x)
    x_neg_at_i = list(x)
    x_pos_at_i[abs(i)-1] = 1
    x_neg_at_i[abs(i)-1] = -1
    fval_at_pos = fast_fval_evaluate(clause,x_pos_at_i,k,ctype)
    fval_at_neg = fast_fval_evaluate(clause,x_neg_at_i,k,ctype)
    if fval_at_pos!="unknown" and fval_at_neg!="unknown":
        return 0.5 * (fval_at_pos-fval_at_neg) 
    return "unknown"
   
def truncate_vector(v,k):
    return v[0:k]
    
def slow_gradient(clause,x,k,ctype,i,FC_table):
    x_rest = [xi/abs(xi)*x[abs(xi)-1] for xi in clause if xi!=i]
    l1 = [FC_table[ci] for ci in range(1,len(clause)+1)]
    l2 = cardcons(x_rest)
    return np.dot(l1,l2) *(abs(i)/i)

 
def grad_BP(clause,x,k,ctype,FC_table):
    grad = []
    nv = len(clause)
    coef = [FC_table[ci] for ci in range(1,len(clause)+1)]
    if nv < 40:
        for i in clause:
            x_rest = [xi/abs(xi)*x[abs(xi)-1] for xi in clause if xi!=i]
            esps = cardcons(x_rest)    
            grad.append(np.dot(coef,esps) * (abs(i)/i))
        return grad
    esps= [] # Elementary symetric Polymomials
    x_prime = [ -x[abs(clause[i])-1] * (clause[i]/abs(clause[i])) for i in range(len(clause))]
    forward_message = []
    backward_message = []
    forward_message.append([1])
    backward_message.append([1])
    for i in range(1,nv):
        forward_message.append(np.polymul(forward_message[i-1],poly1d(x_prime[i-1],True)))
        backward_message.append(np.polymul(backward_message[i-1],poly1d(x_prime[nv-i],True)))
    for i in range(1,nv+1):
        esps.append(np.polymul(forward_message[i-1],backward_message[nv-i]))
    for i in range(nv):
        grad.append(np.dot(coef,esps[i].c) * (abs(clause[i])/clause[i]))
    return grad 


def build_S(n,clauses,x,klist,ctype):
    S = np.zeros((n,n))
    for k in range(len(clauses)):
        clause = clauses[k]
        nv = len(clause)
        if nv < 40:
            for i in clause:
                 x_rest = [xi/abs(xi)*x[abs(xi)-1] for xi in clause if xi!=i]
                 esps = cardcons(x_rest) * (abs(i)/i)
                 S[abs(i)-1,:nv] += esps * abs(i)/i
        else:
            x_prime = [ -x[abs(clause[i])-1] * (clause[i]/abs(clause[i])) for i in range(len(clause))]
            forward_message = []
            backward_message = []
            forward_message.append([1])
            backward_message.append([1])
            for i in range(1,nv):
                forward_message.append(np.polymul(forward_message[i-1],poly1d(x_prime[i-1],True)))
                backward_message.append(np.polymul(backward_message[i-1],poly1d(x_prime[nv-i],True)))
            for i in range(1,nv+1):
                esps = np.polymul(forward_message[i-1],backward_message[nv-i])
                S[abs(clause[i-1])-1,:nv] += esps * abs(clause[i-1]/clause[i-1])
    return S 

def verify_sol(clauses,x,klist,weight,ctype,coefs,comparator):
    x = rounding(x)
    n = len(x)
    num_of_unsat_clauses = 0
    unsat_weight = 0
    i = 0
    ucnf = 0
    uxor = 0
    change_weight = 1
    for l_list in clauses:
        if ctype[i] == 'c':
            if sum(j/abs(j)*x[abs(j)-1] for j in l_list)>len(l_list)-2*klist[i]:
                num_of_unsat_clauses += 1
                ucnf+=1
                unsat_weight += weight[i]
                weight[i] *= change_weight
        elif ctype[i] == 'x':
            if prod(i/abs(i)*x[abs(i)-1] for i in l_list) > 0:
                num_of_unsat_clauses += 1
                unsat_weight += weight[i]
                uxor+=1
                weight[i] *= change_weight
        elif ctype[i] == 'n':
            i0 = l_list[0]
            flag = 1
            v0 = i0/(abs(i0))*x[abs(i0)-1]
            for j in l_list:
                if j/abs(j)*x[abs(j)-1] != v0:
                    flag = 0
                    break
            if flag == 1:
                num_of_unsat_clauses += 1
                unsat_weight += weight[i]
                weight[i] += change_weight * len(l_list)
        elif ctype[i] == 'p':
            s = 0
            for j in range(len(l_list)):
                s += (1-(abs(l_list[j])/l_list[j]) * x[abs(l_list[j])-1])/2 * coefs[i][j] 
            print('s='+repr(s)+' k= '+repr(klist[i]))
            if comparator[i] == '>=' and int(s) < klist[i]:
                num_of_unsat_clauses += 1
                unsat_weight += weight[i]
                weight[i] += change_weight * len(l_list)
            elif comparator[i] == '=' and int(s) != klist[i]:    
                num_of_unsat_clauses += 1
                unsat_weight += weight[i]
                weight[i] += change_weight * len(l_list)
        i+=1
    return num_of_unsat_clauses, unsat_weight, weight


def fun(x,*args):
    x = truncate(x)
    clauses, weight, klist, ctype, FC_table, FC_vector, n = args[0],args[1],args[2], args[3], args[4], args[5], args[6]
    start = time.time()
    fval = 0
    for k in range(len(clauses)):
        #signal = fast_fval_evaluate(clauses[k],x,klist[k],ctype[k])
        #if signal!="unknown":
        #    fval+=weight[k] * signal
        #else:
    #    print(clauses[k])
    #    print(klist[k])
    #    print(ctype[k])
    #    print(FC_table[k])
        fval+=weight[k] * slow_fval_evaluate(clauses[k],x,klist[k],ctype[k],FC_table[k])
        
    end = time.time()     
   # print('fun uses '+repr(end-start))
#    print('x count '+repr(count(x)))
    return fval

def abortable_worker(func,*args,**kwargs):
    timeout = kwargs.get('timeout',None)
    p = ThreadPool(1)
    res = p.apply_async(func,args=args)
    try:
           out = res.get(timeout)
           return out
    except multiprocessing.TimeoutError:
           res = so.OptimizeResult
           res.x = []
           return False,res,True

def gradient_matrix(x,*args):
    x = truncate(x)
    clauses, weight, klist, ctype, FC_table, FC_vector, n = args[0],args[1],args[2], args[3], args[4], args[5], args[6]
    counter = 0
    # matrix methods
    start = time.time()
    S = build_S(n,clauses,x,klist,ctype)
    #print(np.dot(S, FC_vector))
    grad = np.dot(S,FC_vector)   
    end = time.time()
    print('grad matrix takes '+repr(end-start))
    return grad

def gradient(x, *args):
    #print('before grad '+'x='+repr(x))
    switch = 0
    start = time.time()
    x = truncate(x)
    clauses, weight, klist, ctype, FC_table, FC_vector, n = args[0],args[1],args[2], args[3], args[4], args[5], args[6]
    grad = [0 for i in range(len(x))]
    counter = 0          
    signal = "unknown"
    for k in range(len(clauses)):
        signal = "unknown"
        if ctype[k] == 'x':
            xor_g = xor_grad(clauses[k],x)
            for i in range(len(clauses[k])):
                grad[abs(clauses[k][i])-1] += xor_g[i]
        elif ctype[k] == 'c':
            if switch == 0:
                for i in clauses[k]:
                  #  signal = fast_gradient(clauses[k],x,klist[k],ctype[k],i)
                    signal = "unknown"
                    if signal!="unknown":
                        grad[abs(i)-1] += weight[k] * signal
                    else:
                        grad[abs(i)-1] += weight[k] * slow_gradient(clauses[k],x,klist[k],ctype[k],i,FC_table[k])
            else: # belief propogation
                temp_grad = grad_BP(clauses[k],x,klist[k],ctype[k],FC_table[k])
                for i in range(len(clauses[k])):
                    clause = clauses[k]
                    grad[abs(clause[i])-1] += weight[k] * temp_grad[i]
        elif ctype[k] == 'p':
            x_slice = [ x[abs(i)-1] * abs(i)/i for i in clauses[k]] # assume no negative literal in PB constraints
            temp_grad = gradient_on_bdd(x_slice,FC_table[k])
            for i in range(len(clauses[k])):
                clause = clauses[k]
                grad[abs(clause[i])-1] += weight[k] * temp_grad[i] * abs(clause[i]) / clause[i]
    end = time.time() 
    #print('gradient uses '+repr(end-start))
    #print('gradient'+repr(grad))
   # print('x count '+repr(count(x)))
    #print('grad= '+repr(grad))
    return grad

def callback(xk):
    global globalargs
    clauses, weight, klist, ctype, FC_table, FC_vector, n = globalargs[0],globalargs[1],globalargs[2], globalargs[3], globalargs[4], globalargs[5], globalargs[6]
    fval = fun(xk,clauses,weight,klist,ctype,FC_table,FC_vector,n)
#    print('callback fval: '+repr(fval))
    #print('callback x:' +repr(xk))
    return

def use_scipy_kernel(x0,clauses,klist,param,weight,ctype,no,FC_table):
    # Random restart by setting x0
    global globalargs
    numofvars,numofclas = param
    np.random.seed(no+int(time.time()))
    x0 = 2 * np.random.rand(numofvars) - 1
    bnds = ()
    for j in range(numofvars):
        bnds += ((-1, 1),)
    FC_vector = compute_FC_vector(clauses,klist,ctype,weight,numofvars)
    args = (clauses,weight,klist,ctype,FC_table, FC_vector,numofvars)
    #args = (clauses,weight,klist,ctype,FC_table,numofvars)
    #print('gd='+repr(gradient(x0,clauses,weight,klist,ctype,FC_table)))
    opt = {'maxiter':50,'disp':False}
    globalargs = args
    res = so.minimize(fun, x0, method='SLSQP', bounds=bnds,jac=gradient,args=args,options=opt, callback=callback, tol=1e-50)  #Alternative method: L-BFGS-B tol=0.1
    #res = so.minimize(fun, x0, method='SLSQP', bounds=bnds,jac=gradient_matrix,args=args,options=opt, callback=callback, tol=1e-50)  #Alternative method: L-BFGS-B tol=0.1
    optimal = sum(weight)
    return False,res, False

def norm_of_projected_gradient(x,grad):
    newx = truncate([x[i] -grad[i] for i in range(len(grad))])
    norm = sum( (newx[i]-x[i])*(newx[i]-x[i]) for i in range(len(x)))
    return norm

def use_my_kernel(clauses,klist,param,weight,ctype,no,FC_table):
    numofvars,numofclas = param
    np.random.seed(int(time.time())+1000*no)
    x0 = 2 * np.random.rand(numofvars) - 1
    #x0 = [F(x0[i]).limit_denominator(100) for i in range(numofvars)]
    numofvars,numofclas = param
    grad_norm = 1
    x = x0
    res = Opt_result()
    res.x = []
    iter_num = 1
    iteration = 0
    delta = 1
    #while grad_norm>1e-4:
    while delta>1e-4:
        iteration += 1
    #    print('x='+repr(x))
        step_size = 0.1
        iter_num += 1 
        old_fval = fun(x,clauses,weight,klist,ctype,FC_table,numofvars)
        grad = gradient(x,clauses,weight,klist,ctype,FC_table,numofvars)
        grad_norm = norm_of_projected_gradient(x,grad)
     #   print('gnorm=' + repr(grad_norm))
        if grad_norm > 1: normalization = math.sqrt(grad_norm)
        else: normalization = 1
        x = truncate([x[i] -grad[i]/normalization*step_size for i in range(len(grad))])
        fval = fun(x,clauses,weight,klist,ctype,FC_table,numofvars)
        print('fval= '+repr(fval))
        print('old fval= '+repr(old_fval))
        print('count= '+repr(count(x)))
        print('grad norm= '+repr(grad_norm))
        delta = abs(old_fval - fval)
        _,fval,_ = verify_sol(clauses,x,klist,weight,ctype)
        print('delta = '+repr(delta))
        if fval == sum(-1*weight): break
    res.x = x
    res.fun = fval
    optimal = sum(weight)
    unsat_clauses_number,_,_ = verify_sol(clauses,x,klist,weight,ctype)
    print('unsat= ' + repr(unsat_clauses_number))
    print('iterations= ' + repr(iteration))
    if (res.fun + optimal) ** 2 < 1e-10:
        return True, res, False
    else:
        return False,res, False


def use_heavyball_kernel(clauses,klist,param,weight,ctype,no,FC_table):
    numofvars,numofclas = param
    np.random.seed((int(time.time()*1000)+1000*no)%(2**32))
    x0 = 2 * np.random.rand(numofvars) - 1
    #x0 = [F(x0[i]).limit_denominator(100) for i in range(numofvars)]
    numofvars,numofclas = param
    grad_norm = 1
    x = x0
    res = Opt_result()
    res.x = []
    iter_num = 1
    iteration = 0
    delta = 1
    x_last = x.copy()
    MAX_ITER = 100
    #while grad_norm>1e-4:
    while delta>1e-4:
        iteration += 1
    #    print('x='+repr(x))
        step_size = 0.5
        beta = 0.9
        old_fval = fun(x,clauses,weight,klist,ctype,FC_table,[],numofvars)
        grad = gradient(x,clauses,weight,klist,ctype,FC_table,[],numofvars)
     #   grad_norm = norm_of_projected_gradient(x,grad)
     #   print('gnorm=' + repr(grad_norm))
     #   if grad_norm > 1: normalization = math.sqrt(grad_norm)
     #   else: normalization = 1
        x_last_temp = x.copy()
        x = truncate([x[i] -grad[i] * step_size + (x[i]-x_last[i])*beta for i in range(len(grad))])
        x_last = x_last_temp.copy()
        fval = fun(x,clauses,weight,klist,ctype,FC_table,[],numofvars)
     #   print('fval= '+repr(fval))
        delta = abs(old_fval - fval)
        if iteration>MAX_ITER: break
    res.x = x
    res.fun = fval
    optimal = sum(weight)
    if (res.fun + optimal) ** 2 < 1e-10:
        return True, res, False
    else:
        return False,res, False



def use_nes_kernel(clauses,klist,param,weight,ctype,no,FC_table):
    numofvars,numofclas = param
    np.random.seed(int(time.time())+1000*no)
    x0 = 2 * np.random.rand(numofvars) - 1
    numofvars,numofclas = param
    grad_norm = 1
    x = x0
    res = Opt_result()
    res.x = []
    iteration = 0
    delta = 1
    MAXITER = 100
    y = x.copy()
    while delta>1e-4:
        iteration += 1
        step_size = 0.01
        beta = 0.8
        old_fval = fun(x,clauses,weight,klist,ctype,FC_table,[],numofvars)
        grad_y = gradient(y,clauses,weight,klist,ctype,FC_table,[],numofvars)
        x_tplusone = truncate([y[i] - grad_y[i]*step_size for i in range(len(y))])
        y = [x_tplusone[i] + beta*(x_tplusone[i]-x[i]) for i in range(len(y))]
        x = x_tplusone.copy()
        fval = fun(x,clauses,weight,klist,ctype,FC_table,[],numofvars)
        print('fval= '+repr(fval))
        print('old fval= '+repr(old_fval))
        print('count= '+repr(count(x)))
        delta = abs(old_fval - fval)
        _,fval,_ = verify_sol(clauses,x,klist,weight,ctype)
        print('delta = '+repr(delta))
        if fval == sum(-1*weight): break
        if iteration > MAXITER: break
    res.x = x
    res.fun = fval
    optimal = sum(weight)
    unsat_clauses_number,_,_ = verify_sol(clauses,x,klist,weight,ctype)
    print('unsat= ' + repr(unsat_clauses_number))
    print('iterations= ' + repr(iteration))
    if (res.fun + optimal) ** 2 < 1e-10:
        return True, res, False
    else:
        return False,res, False

def SA_generator(x,step_size):
    rand_direction = numpy.random.rand(len(x))
    rand_direction = [1-2*rand_direction[i] for i in range(len(x))]
    norm = math.sqrt(sum([rand_direction[i]*rand_direction[i] for i in range(len(x))]))
    rand_direction = [rand_direction[i]/norm for i in range(len(x))]
    newx = truncate([x[i] + step_size * rand_direction[i] for i in range(len(x))])
    return newx

def use_SA_kernel(clauses,klist,param,weight,ctype,no,FC_table):
    numofvars,numofclas = param
    np.random.seed(int(time.time())+1000*no)
    x0 = 2 * np.random.rand(numofvars) - 1
    #x0 = [F(x0[i]).limit_denominator(100) for i in range(numofvars)]
    numofvars,numofclas = param
    x = x0
    res = Opt_result()
    res.x = []
    iteration = 0
    delta = 1
    MAXITER = 30000
    gamma = 10000.0
    m0 = 100.0
    temperature = gamma / math.log2(m0)
    while iteration<MAXITER:
        print('iter= '+repr(iteration))
        iteration += 1
        step_size = 0.5
        decay = 1.0/(len(x))
        old_fval = fun(x,clauses,weight,klist,ctype,FC_table,numofvars)
        newx = SA_generator(x,step_size)
        fval = fun(newx,clauses,weight,klist,ctype,FC_table,numofvars)
        try:
            prob = min(1,math.exp((old_fval-fval)/(temperature)))
        except OverflowError:
            if fval<old_fval: prob =1-1e-3
            else: prob =1e-3
        print('temperature= '+repr(temperature))
        #prob = 1/(1+math.exp((-old_fval+fval)/(1e-5+temperature)))
        print('prob= '+repr(prob))
        if prob>random.random():
           x = newx
        #temperature = gamma/math.log2(iteration+m0)
        temperature = gamma*math.exp(-iteration*decay)
        print('old fval= '+repr(old_fval))
        print('fval= '+repr(fval))
        print('count= '+repr(count(x)))
        delta = abs(old_fval - fval)
        _,fval,_ = verify_sol(clauses,x,klist,weight,ctype)
        print('delta = '+repr(delta))
        if fval == sum(-1*weight): break
    res.x = x
    res.fun = fval
    optimal = sum(weight)
    unsat_clauses_number,_,_ = verify_sol(clauses,x,klist,weight,ctype)
    print('unsat= ' + repr(unsat_clauses_number))
    print('iterations= ' + repr(iteration))
    if (res.fun + optimal) ** 2 < 1e-10:
        return True, res, False
    else:
        return False,res, False



def use_adam_kernel(clauses,klist,param,weight,ctype,no,FC_table):
    numofvars,numofclas = param
    np.random.seed(int(time.time())+1000*no)
    x0 = 2 * np.random.rand(numofvars) - 1
    numofvars,numofclas = param
    grad_norm = 1
    x = x0
    res = Opt_result()
    res.x = []
    iteration = 0
    delta = 1
    MAXITER = 200
    m = np.zeros(numofvars)
    mp = np.zeros(numofvars)
    v = np.zeros(numofvars)
    vp = np.zeros(numofvars)
    while delta>1e-4:
        iteration += 1
        step_size = 0.1
        beta1 = 0.9
        beta2 = 0.999
        old_fval = fun(x,clauses,weight,klist,ctype,FC_table,[],numofvars)
        grad = gradient(x,clauses,weight,klist,ctype,FC_table,[],numofvars)

        m = [ beta1 * m[i] + (1-beta1) * grad[i] for i in range(numofvars) ]
        v = [ beta2 * v[i] + (1-beta2) * grad[i] * grad[i] for i in range(numofvars) ]
        mp = [ m[i]/(1-beta1**iteration) for i in range(numofvars) ]
        vp = [ v[i]/(1-beta2**iteration) for i in range(numofvars) ]
             
        x = truncate([x[i] - mp[i]*step_size/(1e-8+math.sqrt(vp[i])) for i in range(numofvars)])

        fval = fun(x,clauses,weight,klist,ctype,FC_table,[],numofvars)
        print('fval= '+repr(fval))
        print('old fval= '+repr(old_fval))
        print('count= '+repr(count(x)))
        delta = abs(old_fval - fval)
        _,fval,_ = verify_sol(clauses,x,klist,weight,ctype)
        print('delta = '+repr(delta))
        if fval == sum(-1*weight): break
        if iteration > MAXITER: break
    res.x = x
    res.fun = fval
    optimal = sum(weight)
    unsat_clauses_number,_,_ = verify_sol(clauses,x,klist,weight,ctype)
    print('unsat= ' + repr(unsat_clauses_number))
    print('iterations= ' + repr(iteration))
    if (res.fun + optimal) ** 2 < 1e-10:
        return True, res, False
    else:
        return False,res, False



def use_LS_kernel(clauses,klist,param,weight,ctype,no,table):
    numofvars,numofclas = param
    np.random.seed(int(time.time())+1000*no)
    x0 = 2 * np.random.rand(numofvars) - 1
    x0 = rounding(x0)
    numofvars,numofclas = param
    x = x0
    res = Opt_result()
    res.x = list(x)
    _,fval,_ = verify_sol(clauses,x,klist,weight,ctype)
    fval = sum(weight)-fval
    delta = 1
    while delta>0:
        old_fval = fval
        x,fval = LS(clauses,x,klist,weight,ctype)
        delta = fval - old_fval
    res.x = x
    res.fun = fval
    optimal = sum(weight)
    if (res.fun + optimal) ** 2 < 1e-10:
        return True, res, False
    else:
        return False,res, False


def solve_by_OptSat(filepath,timelimit,tolerance,cpus,verbose,ismaxsat):
        formula = Formula()
        formula = Formula.read_DIMACS(filepath)
        numofvars = len(formula._variables)
        numofcla = len(formula.clauses)
        param = numofvars, numofcla
        weight = formula._weight
    #    print('weight='+repr(weight))
        coefs = formula._coefs
        comparator = formula._comparator
        klist = formula._klist
        ctype = formula._ctype
        elapsed = 0
        bestx = []
        if ismaxsat == 1:
            weight = [1 for i in range(numofcla)]
        start = time.time()
        FC_table = compute_FC_table(formula.clauses,klist,ctype,coefs,comparator)
    #    print('compute Fourier table takes '+repr(time.time()-start))
        if verbose == 1:
            xlist = []
        best_unsat_num = numofcla
        total_time_limit = timelimit*0.8
        solved_flag = 0
        weight_group = []
        for _ in range(cpus):
            weight_group.append(weight)
        num_of_trials = 0
        start = time.time()
        np.random.seed(1)
        x0 = 2 * np.random.rand(numofvars) - 1
        while elapsed<total_time_limit:
            results = []
            pool = multiprocessing.Pool()
            for cpu_index in range(cpus):
                time_left = total_time_limit - elapsed 
                abortable_func = partial(abortable_worker,use_scipy_kernel,timeout = time_left)
                #abortable_func = partial(abortable_worker,use_my_kernel,timeout = time_left)
                #abortable_func = partial(abortable_worker,use_nes_kernel,timeout = time_left)
                #abortable_func = partial(abortable_worker,use_adam_kernel,timeout = time_left)
                #abortable_func = partial(abortable_worker,use_SA_kernel,timeout = time_left)
                #abortable_func = partial(abortable_worker,use_heavyball_kernel,timeout = time_left)
                #abortable_func = partial(abortable_worker,use_LS_kernel,timeout = time_left)
                #print('feed x0='+repr(x0))
                result = pool.apply_async(abortable_func, args=(x0,formula.clauses,klist,param,weight_group[cpu_index],ctype,cpu_index,FC_table))
                results.append(result)
                num_of_trials += 1
            pool.close()
            pool.join()
            
            for cpu_index in range(cpus):
                solved_flag, res,timeout_flag = results[cpu_index].get()
                x0 = res.x
                #print('before x0='+repr(x0))
                #x0 = relax(x0)
                #x0 = project(turblance  + x0)
                #print('after x0='+repr(x0))
                if timeout_flag == True or len(res.x)<numofvars:
                    res.x = [1 for _ in range(numofvars)]
                if verbose == 1:
                    xlist.append(rounding(res.x))
                num_of_unsat, weight_of_unsat, weight_group[cpu_index] = verify_sol(formula.clauses,res.x,klist,weight_group[cpu_index],ctype, coefs, comparator)
                
                if verbose == 1:
                #    print('x= '+repr(res.x))
                #    print('v '),
                    res.x = rounding(res.x)
                #    for i in range(len(res.x)):
                #        print(repr(int(-1*res.x[i]*(i+1)))+' '),
                #    print('\n')
                if num_of_unsat == 0:
                	bestx = rounding(res.x)
                	best_unsat_num = 0
                	break
                else:
                    if num_of_unsat < best_unsat_num:
                        bestx = rounding(res.x)
                        best_unsat_num = num_of_unsat
                        if ismaxsat == 1:
                            print('o '+repr(num_of_unsat))
                if verbose == 1:
                    print('current #UNSAT clauses: '+repr(num_of_unsat))
                    print('best #UNSAT clauses: '+repr(best_unsat_num))
            elapsed = time.time() - start
            bestx = rounding(bestx)
                 
            if best_unsat_num <= tolerance:
                solved_flag = 1
                print('s Solved')
                if verbose == 1:
                    print('t '+repr(elapsed)+' with '+repr(num_of_trials)+' trials. Average time per trial = '+repr(elapsed/num_of_trials))
                print('v '),
                for i in range(len(bestx)):
                    print(repr(int(-1*bestx[i]*(i+1)))+' '),
                print('\n')
                break
        if solved_flag==0 and ismaxsat==0:
            print('s Not Solved in ' + repr(timelimit) + ' seconds and ' + repr(num_of_trials) + ' trials. Minimum number of unSAT clauses = ' + repr(best_unsat_num))
            if verbose == 1:
                print('v '),
                for i in range(len(bestx)):
                    print(repr(int(-1*bestx[i]*(i+1)))+' '),
                print('\n')
           
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('filepath',type=str,help='set the file path')
    parser.add_argument('--timelimit', type=int, help='set the time limit')
    parser.add_argument('--tolerance',type=int, help='set the number of clauses allowed to be unsatisfied')
    parser.add_argument('--cpus',type=int,help='set the number of cores')
    parser.add_argument('--verbose',type=int,help='set verbose to 1 to output more information')
    parser.add_argument('--ismaxsat',type=int,help='set to 1 if you are solving an unweighted MAXSAT problem')
    args = parser.parse_args()
    if not args.timelimit:
        timelimit = 60
    else:
        timelimit = args.timelimit
    if not args.filepath:
        print('Please provide a filepath!')
        exit(0)
    if not args.tolerance:
        tolerance = 0
    else:
        tolerance = args.tolerance
    if not args.cpus:
        cpus = 1
    else:
        cpus = args.cpus
    if not args.verbose:
        verbose = 0
    else:
        verbose = 1
    if not args.ismaxsat:
        ismaxsat = 0
    else:
        ismaxsat = 1
    solve_by_OptSat(args.filepath,timelimit,tolerance,cpus,verbose,ismaxsat)

