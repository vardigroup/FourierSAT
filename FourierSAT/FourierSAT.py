import scipy.optimize as so
import numpy as np
from scipy import poly1d
import time
import math
import random
from FourierCoefficient import CardinalityFC
from FourierCoefficient import NAEFC
from boolean_formula import Formula
import signal
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool
from functools import partial
import argparse

def norm(x):
    return sum([x[i]*x[i] for i in range(len(x))])

def prod(iterable):
    p= 1
    for n in iterable:
        p *= n
    return p
def cardcons(x):
    templist = [-x[i] for i in range(len(x))]
    poly = poly1d(templist,True)
    return poly.c
def rounding(x):
    for i in range(len(x)):
        if x[i]>0: x[i]=1
        else: x[i]=-1
    return x
def count(x):
    x = rounding(x)
    return (len(x)-sum(x))//2
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
def verify_sol(formula,x,klist,weight,ctype):
    clauses = formula.clauses
    x = rounding(x)
    n = len(x)
    num_of_unsat_clauses = 0
    unsat_weight = 0
    i = 0
    ucnf = 0
    uxor = 0
    for l_list in clauses:
        if ctype[i] == 'c':
            if sum(j/abs(j)*x[abs(j)-1] for j in l_list)>len(l_list)-2*klist[i]:
                num_of_unsat_clauses += 1
                ucnf+=1
                unsat_weight += weight[i]
                weight[i]*=len(l_list)
        elif ctype[i] == 'x':
            if prod(i/abs(i)*x[abs(i)-1] for i in l_list) > 0:
                num_of_unsat_clauses += 1
                unsat_weight += weight[i]
                uxor+=1
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
        i+=1
    return num_of_unsat_clauses, unsat_weight, weight
def fun(x,*args):
    clauses, weight, klist, ctype = args[0], args[1], args[2], args[3]
    cnftotal = sum(weight[k]*np.dot(cardcons([j/abs(j)*x[abs(j)-1] for j in clauses[k]]),CardinalityFC(len(clauses[k]),klist[k])) for k in range(len(clauses)) if (ctype[k]=='c'))
    naetotal = sum(weight[k] * np.dot(cardcons([j / abs(j) * x[abs(j) - 1] for j in clauses[k]]),
                                      NAEFC(len(clauses[k]))) for k in range(len(clauses)) if
                   (ctype[k] == 'n'))
    xortotal = sum(weight[k] * prod(j/abs(j)*x[abs(j)-1] for j in clauses[k]) for k in range(len(clauses)) if (ctype[k]=='x'))
    fval = cnftotal+xortotal+naetotal
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
def gradient(x, *args):
    clauses, weight, klist, ctype = args[0],args[1],args[2], args[3]
    grad = [0 for i in range(len(x))]
    for k in range(len(clauses)):
        if ctype[k] == 'c':
            card_coef = CardinalityFC(len(clauses[k]),klist[k])
            for i in clauses[k]:
                x_rest = [xi/abs(xi)*x[abs(xi)-1] for xi in clauses[k] if xi!=i]
                grad[abs(i)-1]+=np.dot([card_coef[ci] for ci in range(1,len(clauses[k])+1)],\
                                                cardcons(x_rest)) *(abs(i)/i) * weight[k]
        elif ctype[k] == 'x':
            for i in clauses[k]:
                x_rest = [xi/abs(xi)*x[abs(xi)-1] for xi in clauses[k] if xi!=i]
                grad[abs(i)-1]+=prod(x_rest)*(abs(i)/i)*weight[k]
        if ctype[k] == 'n':
            nae_coef = NAEFC(len(clauses[k]))
            for i in clauses[k]:
                x_rest = [xi/abs(xi)*x[abs(xi)-1] for xi in clauses[k] if xi!=i]
                grad[abs(i)-1]+=np.dot([nae_coef[ci] for ci in range(1,len(clauses[k])+1)],\
                                                cardcons(x_rest)) *(abs(i)/i) * weight[k]
    return grad

def callback(xk):
    pass

def use_scipy_kernel(clauses,klist,param,weight,ctype,no):
    # Random restart by setting x0
    numofvars,numofclas = param
    np.random.seed(int(time.time())+1000*no)
    x0 = 2 * np.random.rand(numofvars) - 1
    numofvars,numofclas = param
    bnds = ()
    for j in range(numofvars):
        bnds += ((-1, 1),)
    opt = {'maxiter':100}#,'disp':True}
    args = (clauses,weight,klist,ctype)
    res = so.minimize(fun, x0, method='SLSQP', bounds=bnds,jac=gradient,args=args,options=opt, callback=callback, tol=1e-20)  #Alternative method: L-BFGS-B tol=0.1
    opt = sum(weight)
    if (res.fun + opt) ** 2 < 1e-10:
        return True, res, False
    else:
        return False,res, False

def solve_by_OptSat(filepath,timelimit,tolerance,cpus,verbose):
        formula = Formula()
        formula = Formula.read_DIMACS(filepath)
        numofvars = len(formula._variables)
        numofcla = len(formula.clauses)
        param = numofvars, numofcla
        weight = formula._weight
        klist = formula._klist
        ctype = formula._ctype
        elapsed = 0
        start = time.time()
        bestx = []
        best_unsat_num = numofcla
        total_time_limit = timelimit*0.8
        solved_flag = 0
        weight_group = []
        for _ in range(cpus):
            weight_group.append(weight)
        while elapsed<total_time_limit:
            results = []
            pool = multiprocessing.Pool()
            for cpu_index in range(cpus):
                time_left = total_time_limit - elapsed 
                abortable_func = partial(abortable_worker,use_scipy_kernel,timeout = time_left)
                #abortable_func = partial(abortable_worker,use_my_kernel,timeout = time_left)
                result = pool.apply_async(abortable_func, args=(formula.clauses,klist,param,weight_group[cpu_index],ctype,cpu_index))
                results.append(result)
            pool.close()
            pool.join()
            for cpu_index in range(cpus):
                solved_flag, res,timeout_flag = results[cpu_index].get()
                if timeout_flag == True or len(res.x)<numofvars:
                    res.x = [1 for _ in range(numofvars)]
                num_of_unsat, weight_of_unsat, weight_group[cpu_index] = verify_sol(formula,res.x,klist,weight_group[cpu_index],ctype)
                if num_of_unsat == 0:
                	bestx = res.x
                	best_unsat_num = 0
                	break
                else:
                	if num_of_unsat < best_unsat_num:
                    		bestx = res.x
                    		best_unsat_num = num_of_unsat
            if verbose == 1:
                print('#UNSAT clauses: '+repr(best_unsat_num))
            if best_unsat_num <= tolerance:
                solved_flag = 1
                print('s Solved')
                print('v '),
                for i in range(len(bestx)):
                    print(repr(int(-1*bestx[i]*(i+1)))+' '),
                print('\n')
                break
            elapsed = time.time() - start
        if solved_flag==0:
            print('s Not Solved in ' + repr(timelimit) + ' seconds. Minimum number of unSAT clauses = ' + repr(best_unsat_num))
            if verbose == 1:
                print('v '),
                for i in range(len(bestx)):
                    print(repr(int(-1*bestx[i]*(i+1)))+' '),
                print('\n')
                

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filepath',type=str,help='set the file path')
    parser.add_argument('--timelimit', type=int, help='set the time limit')
    parser.add_argument('--tolerance',type=int, help='set the number of clauses allowed to be unsatisfied')
    parser.add_argument('--cpus',type=int,help='set the number of cores')
    parser.add_argument('--verbose',type=int,help='set verbose to 1 to output more information')
    args = parser.parse_args()
    if not args.timelimit:
        timelimit = 10
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
    solve_by_OptSat(args.filepath,timelimit,tolerance,cpus,verbose)

if __name__ == "__main__":
    main()
