# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 17:39:08 2020

@author: jim903
"""
# This file uses SymPy to find the coefficients of Lambda_p and save it to the file


import numpy as np
#Extend to 4, 6, ... orders. 
order = 6
from sympy import symbols, integrate, series, exp, Poly, Function, linsolve, log
#from sympy.utilites.lambdify import implemented_function
x, V, sigb, S, sig, Lambdap = symbols('x V sigb S sig Lamdap', positive = True)

tau = Function('tau')(sigb)

a = symbols('a')
Cs = [symbols('C%d' % i) for i in np.arange(2, order + 1, 2)]
Ds = [[] for _ in range(int(order / 2))]
for _ in range(int(order / 2)):
    Ds[_].extend([symbols('D%d_%d' %(_ * 2 + 2, i)) for i in range(_ * 2 + 2)])


Eqns = [[] for _ in range(int(order / 2))]
solns = [[] for _ in range(int(order / 2))]
for i in range(order + 1):
    print(i)
    Exp = series((sigb - S) ** i - (((sigb + S) * (1 + a) / 2 - S) ** i 
                  + ((sigb + S) * (1 - a) / 2 - S) ** i) 
    * exp(- Lambdap * tau) , x = sigb, x0 = S, n = order + 1)
    Exp = Exp.removeO()
    
    lam_exp = 1
    for j in range(len(Cs)):
        lam_exp += Cs[j] * a ** (2 * (j + 1))
    Exp = Exp.subs(Lambdap, log(2) / tau.subs(sigb, S) * lam_exp)
    const = series(Exp, x = sigb, x0 = S, n = 1).removeO()
    
    Exp2 = const
    for k in range(int(order / 2)):
        for m in range(k, int(order / 2)):
            Exp2 += Exp.coeff(sigb - S, 2 * k + 1) * Ds[m][2 * k] * a ** (2 * m + 2)
            Exp2 += Exp.coeff(sigb - S, 2 * k + 2) * Ds[m][2 * k + 1] * a ** (2 * m + 2)
    
    Exp2 = series(Exp2, x = a, x0 = 0, n = order + 1).removeO()
    for m in np.arange(int(order / 2)):
        coef = Exp2.coeff(a, m * 2 + 2)
        if coef != 0:
            p = Poly(coef, tuple([Cs[m]] + Ds[m]))
            pdoit = p.doit()
            pexp = pdoit.as_expr()
            Eqns[m].append(pexp)

for m in np.arange(int(order / 2)):
    print(m)
    solns[m] = list(linsolve(Eqns[m], tuple([Cs[m]] + Ds[m])).args[0])


for m in np.arange(1, int(order / 2)):
    for sol_ind in range(len(solns[m])):
        for l in range(m):
            solns[m][sol_ind] = solns[m][sol_ind].subs(Cs[l], solns[l][0])
            for k in range(len(Ds[l])):
                solns[m][sol_ind] = solns[m][sol_ind].subs(Ds[l][k], solns[l][1+k])

Lambdap_exp = log(2) / tau.subs(sigb, S) 
for i in range(len(Cs)):
    Lambdap_exp += log(2) / tau.subs(sigb, S) * solns[i][0].expand() * a ** (2 * i + 2)

lam0, lam1, n = symbols('lam0 lam1 n')
lam_prl = lambda t: (lam0 + lam1 * t ** n) / (1 + t ** n)
tau_prl = lambda x: integrate(1 / V / lam_prl((x - S) / V + S), (V , 1, 2))
Lambdap_exp_prl = Lambdap_exp.replace(Function('tau')(S), tau_prl(S))

p = Poly(Lambdap_exp_prl, a)
Clst = p.all_coeffs()
C6_S = Clst[-7].subs([(lam0, 1), (lam1, 0), (n, 2)]).doit().evalf()
C4_S = Clst[-5].subs([(lam0, 1), (lam1, 0), (n, 2)]).doit().evalf()
C2_S = Clst[-3].subs([(lam0, 1), (lam1, 0), (n, 2)]).doit().evalf()
C0_S = Clst[-1].subs([(lam0, 1), (lam1, 0), (n, 2)]).doit().evalf()

for s in np.arange(0.01, 1.01, 0.01):
    C6 = C6_S.subs([(S, s)]).doit().evalf()
    C4 = C4_S.subs([(S, s)]).doit().evalf()
    C2 = C2_S.subs([(S, s)]).doit().evalf()
    C0 = C0_S.subs([(S, s)]).doit().evalf()
    
    np.savetxt('C_damage_S={:.2f}.out'.format(s), [C0, C2, C4, C6])

C6_S = Clst[-7].subs([(lam0, 0.2), (lam1, 1.2), (n, 2)]).doit().evalf()
C4_S = Clst[-5].subs([(lam0, 0.2), (lam1, 1.2), (n, 2)]).doit().evalf()
C2_S = Clst[-3].subs([(lam0, 0.2), (lam1, 1.2), (n, 2)]).doit().evalf()
C0_S = Clst[-1].subs([(lam0, 0.2), (lam1, 1.2), (n, 2)]).doit().evalf()

for s in np.arange(0.01, 1.01, 0.01):
    C6 = C6_S.subs([(S, s)]).doit().evalf()
    C4 = C4_S.subs([(S, s)]).doit().evalf()
    C2 = C2_S.subs([(S, s)]).doit().evalf()
    C0 = C0_S.subs([(S, s)]).doit().evalf()
    
    np.savetxt('C_benefit_S={:.2f}.out'.format(s), [C0, C2, C4, C6])
    
            
## As an exmample, try plotting for different S's for damage case
#Slst = np.arange(0.3, 0.7, 0.1)
#alst = np.arange(0, 1.02, 0.02)
#C6_S = Clst[-7].subs([(lam0, 1), (lam1, 0), (n, 2)]).doit().evalf()
#C4_S = Clst[-5].subs([(lam0, 1), (lam1, 0), (n, 2)]).doit().evalf()
#C2_S = Clst[-3].subs([(lam0, 1), (lam1, 0), (n, 2)]).doit().evalf()
#C0_S = Clst[-1].subs([(lam0, 1), (lam1, 0), (n, 2)]).doit().evalf()
#
#import matplotlib.pyplot as plt
#plt.figure()
#Lambdap_mat_4th = np.zeros((len(Slst), len(alst)))
#Lambdap_mat_6th = np.zeros((len(Slst), len(alst)))
#
#
#for Sind in range(len(Slst)):
#    C6 = C6_S.subs([(S, Slst[Sind])]).doit().evalf()
#    C4 = C4_S.subs([(S, Slst[Sind])]).doit().evalf()
#    C2 = C2_S.subs([(S, Slst[Sind])]).doit().evalf()
#    C0 = C0_S.subs([(S, Slst[Sind])]).doit().evalf()
#
#    print(C6)
#
#    Lambdap0 = lam_prl(Slst[Sind]).subs([(lam0, 1), (lam1, 0), (n, 2)])
#    for aind in range(len(alst)):
#        Lambdap_val_4th = C0 + C2 * alst[aind] ** 2 + C4 * alst[aind] ** 4
#        Lambdap_val_6th = Lambdap_val_4th + C6 * alst[aind] ** 6
##        Lambdap_val = Lambdap_exp_prl.subs([(lam0, 1), (lam1, 0), (n, 2), (S, Slst[Sind]), (a, alst[aind])]).doit().evalf()
##        Lambdap_mat[Sind][aind] = Lambdap_val
#        Lambdap_mat_4th[Sind][aind] = Lambdap_val_4th
#        Lambdap_mat_6th[Sind][aind] = Lambdap_val_6th
#
#    plt.plot(alst, Lambdap_mat_4th[Sind] - Lambdap_mat_4th[Sind][0]
#    , label = '4th order, S = ' + str(np.round(Slst[Sind], 2))
#    , color=plt.cm.RdYlBu(Sind), linestyle = 'dashed')
#    plt.plot(alst, Lambdap_mat_6th[Sind] - Lambdap_mat_6th[Sind][0]
#    , label = '6th order, S = ' + str(np.round(Slst[Sind], 2))
#    , color=plt.cm.RdYlBu(Sind), alpha = 0.5)
#
#    max_ind_4th = np.argmax(Lambdap_mat_4th[Sind])
#    plt.scatter(alst[max_ind_4th], Lambdap_mat_4th[Sind][max_ind_4th] - Lambdap_mat_4th[Sind][0], c = 'r', marker = '*')
#    max_ind_6th = np.argmax(Lambdap_mat_6th[Sind])
#    plt.scatter(alst[max_ind_6th], Lambdap_mat_6th[Sind][max_ind_6th] - Lambdap_mat_6th[Sind][0], c = 'g', marker = '*')
#
#    print(str(alst[max_ind_6th]) + ', ' +str(alst[max_ind_6th]))
#
#plt.legend()
#plt.xlabel('a')
#plt.ylabel('$\Lambda_p - \Lambda_p(a = 0)$')
#plt.title('damage, order = ' + str(order))
#
#
## Plot a_c as function of S
#Slst = np.arange(0.01, 1.01, 0.01)
#a_lst = np.arange(0, 1.001, 0.001)
#ac_max_list = []
#for s in Slst:
#    C6 = float(C6_S.subs([(S, s)]).doit().evalf())
#    C4 = float(C4_S.subs([(S, s)]).doit().evalf())
#    C2 = float(C2_S.subs([(S, s)]).doit().evalf())
#    C0 = float(C0_S.subs([(S, s)]).doit().evalf())
#    Lambdap_lst = C0 + C2 * a_lst ** 2 + C4 * a_lst ** 4 + C6 * a_lst ** 6
#    ac_ind = np.argmax(Lambdap_lst)
#    ac = a_lst[ac_ind]
#    ac_max_list.append(ac)
#plt.figure()
#plt.plot(Slst, ac_max_list)
#plt.xlabel('S')
#plt.ylabel('$a_c$')
#plt.title('damage')
#
#
#
#
#Slst = np.arange(0.3, 0.7, 0.1)
#alst = np.arange(0, 1.01, 0.01)
#C6_S = Clst[-7].subs([(lam0, 0.2), (lam1, 1.2), (n, 2)]).doit().evalf()
#C4_S = Clst[-5].subs([(lam0, 0.2), (lam1, 1.2), (n, 2)]).doit().evalf()
#C2_S = Clst[-3].subs([(lam0, 0.2), (lam1, 1.2), (n, 2)]).doit().evalf()
#C0_S = Clst[-1].subs([(lam0, 0.2), (lam1, 1.2), (n, 2)]).doit().evalf()
#plt.figure()
#Lambdap_mat_4th = np.zeros((len(Slst), len(alst)))
#Lambdap_mat_6th = np.zeros((len(Slst), len(alst)))
#for Sind in range(len(Slst)):
#    C6 = C6_S.subs([(S, Slst[Sind])]).doit().evalf()
#    C4 = C4_S.subs([(S, Slst[Sind])]).doit().evalf()
#    C2 = C2_S.subs([(S, Slst[Sind])]).doit().evalf()
#    C0 = C0_S.subs([(S, Slst[Sind])]).doit().evalf()
#    print(C6)
#
#    Lambdap0 = lam_prl(Slst[Sind]).subs([(lam0, 0.2), (lam1, 1.2), (n, 2)])
#    for aind in range(len(alst)):
#        Lambdap_val_4th = C0 + C2 * alst[aind] ** 2 + C4 * alst[aind] ** 4
#        Lambdap_val_6th = Lambdap_val_4th + C6 * alst[aind] ** 6
#        Lambdap_mat_4th[Sind][aind] = Lambdap_val_4th
#        Lambdap_mat_6th[Sind][aind] = Lambdap_val_6th
#
#    plt.plot(alst, Lambdap_mat_4th[Sind] - Lambdap_mat_4th[Sind][0]
#    , label = '4th order, S = ' + str(np.round(Slst[Sind], 2))
#    , color=plt.cm.RdYlBu(Sind), linestyle = 'dashed')
#    plt.plot(alst, Lambdap_mat_6th[Sind] - Lambdap_mat_6th[Sind][0]
#    , label = '6th order, S = ' + str(np.round(Slst[Sind], 2))
#    , color=plt.cm.RdYlBu(Sind), alpha = 0.5)
#
#    max_ind_4th = np.argmax(Lambdap_mat_4th[Sind])
#    plt.scatter(alst[max_ind_4th], Lambdap_mat_4th[Sind][max_ind_4th] - Lambdap_mat_4th[Sind][0], c = 'r', marker = '*')
#    max_ind_6th = np.argmax(Lambdap_mat_6th[Sind])
#    plt.scatter(alst[max_ind_6th], Lambdap_mat_6th[Sind][max_ind_6th] - Lambdap_mat_6th[Sind][0], c = 'g', marker = '*')
#
#    print(str(alst[max_ind_6th]) + ', ' +str(alst[max_ind_6th]))
#
#plt.legend()
#plt.xlabel('a')
#plt.ylabel('$\Lambda_p - \Lambda_p(a = 0)$')
#plt.title('benefit, order = ' + str(order))
#
#
## Plot a_c as function of S
#Slst = np.arange(0.01, 1.01, 0.01)
#a_lst = np.arange(0, 1.001, 0.001)
#ac_max_list = []
#for s in Slst:
#    C6 = float(C6_S.subs([(S, s)]).doit().evalf())
#    C4 = float(C4_S.subs([(S, s)]).doit().evalf())
#    C2 = float(C2_S.subs([(S, s)]).doit().evalf())
#    C0 = float(C0_S.subs([(S, s)]).doit().evalf())
#    Lambdap_lst = C0 + C2 * a_lst ** 2 + C4 * a_lst ** 4 + C6 * a_lst ** 6
#    ac_ind = np.argmax(Lambdap_lst)
#    ac = a_lst[ac_ind]
#    ac_max_list.append(ac)
#
#plt.figure()
#plt.plot(Slst, ac_max_list)
#plt.xlabel('S')
#plt.ylabel('$a_c$')
#plt.title('benefit')
#
#
#sigA = 0
#sigS = 0
#sigT = 0
#DorB = 'B'
#Slst_numeric = np.arange(0.43, 0.58, 0.01)
#Slst_numeric = np.append(Slst_numeric, np.arange(0.05, 0.5, 0.1))
#Slst_numeric = np.append(Slst_numeric, np.arange(0.55, 1, 0.1))
#a_max_lst = []
#for Sind in range(len(Slst_numeric)):
#    S = np.round(Slst_numeric[Sind], 2)
#
#    lambda_file_ref = 'extended_Lambda_p_S={:.2f},sigA={:.3f},sigS={:.3f},sigT={:.3f},{}.out'.format(
#                S, sigA, sigS, sigT, DorB)
#    Lambda_p_ref = np.loadtxt(lambda_file_ref)
#    Lambda_a_0 = Lambda_p_ref[0]
#    a_file = 'extended_Alst_S={:.2f},sigA={:.3f},sigS={:.3f},sigT={:.3f},{}.out'.format(
#                S, sigA, sigS, sigT, DorB)
#    a_extended_lst = np.loadtxt(a_file)
#
#    max_ind = np.argmax(Lambda_p_ref)
#    a_max = a_extended_lst[max_ind]
#    a_max_lst.append(a_max)
#
#plt.scatter(Slst_numeric, a_max_lst)
