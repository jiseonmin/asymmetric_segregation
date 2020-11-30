# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 15:08:49 2020

@author: jim903
"""

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

LAM0 = 0.2
LAM1 = 1.2
N = 2
Slist = np.arange(0.05, 1.05, 0.05)

tol = 1e-3
num_I = 50
klist = range(50)

lam = lambda sigma : (LAM0 + LAM1 * sigma ** N) / (1 + sigma ** N)

def tau(k, S):
    I = lambda V : 1 / V / lam((k * S - S) / V + S)
    res = quad(I, 1, 2)
    return(res[0])

def tau_prime(k, S):
    dk = 0.001
    res = (tau(k + dk, S) - tau(k - dk, S)) / (2 * dk * S)
    return(res)
    
def I(lambdap_1):
    res = 0
    for i in range(num_I + 1):
        res += np.exp(-lambdap_1 * np.sum([tau(j, S) for j in range(i + 1)]))
    return(res)


C1list = []
for S in Slist:
    # First find Lambdap_1
    Lambdap_1 = 1
    n_it = 0
    while np.abs(I(Lambdap_1) - 1) > tol:
        if I(Lambdap_1) > 1:
            Lambdap_1 += 1e-3
            n_it += 1
        else:
            Lambdap_1 -= 1e-3
            n_it += 1
        if n_it > 1000:
            break
    
    def mean(n, S):
        mean0 = np.sum([(i + 1) / 2 * S 
                        * np.exp(-Lambdap_1 * 
                                 np.sum([tau(j, S) for j in range(i + 1)])) for i in range(num_I + 1)])
        if n == 0:
            return(mean0)
        else:
            mean = mean0 - S * n * (n + 1) / 4
            return(mean)
        
    C1_num = np.sum([np.exp(-Lambdap_1 
                            * np.sum([tau(j, S) for j in range(i + 1)]))
                    * np.sum([tau_prime(k, S) * mean(k, S) for k in range(i + 1)]) 
                    for i in klist])
    C1_den = np.sum([np.exp(-Lambdap_1 * np.sum([tau(j, S) for j in range(i + 1)]))
            * np.sum([tau(k, S) for k in range(i + 1)]) for i in klist])
    C1 = - C1_num / C1_den
    C1list.append(C1)

plt.figure()
plt.plot(Slist, C1list)
plt.plot(Slist, np.zeros(len(Slist)), color = 'k', linewidth = 0.5)
plt.xlabel('S')
plt.ylabel('$C_1$')
