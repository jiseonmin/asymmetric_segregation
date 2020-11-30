"""
Simulate the population with a that depends on sigma_d. Specifically, we say
a[D] = tanh(k D). k is the slope at D = 0. Assume everything to be deterministic.
Scan S = 0.5 to 1.
"""

import numpy as np
import sys
import scipy.integrate as integrate
# use age instead of volume/protein to check when the cell divides
# don't save vblst since it is too expensive memorywise

N = 2
DT = 0.001
MAX_POP = 5e6
T_SAMPLE = 5
CELLS = 100 # number of cells at time zero

Slst = np.arange(0.02, 1.002, 0.02)
k = float(sys.argv[1])
DorB = str(sys.argv[2])


if DorB == 'D':
    Lambda0 = 1
    Lambda1 = 0
elif DorB == 'B':
    Lambda0 = 0.2
    Lambda1 = 1.2


def main():
    Lambdap_lst = []
    for S in Slst:
        print(S)
        Lambdap = build_ylst(S)
        print(Lambdap)
        Lambdap_lst.append(Lambdap)

    np.savetxt(
    'Lambda_p_adaptive_a_k={:.2f},{}.out'.format
    (k, DorB), Lambdap_lst
            )
    np.savetxt('Slst_k={:.2f},{}.out'.format
    (k, DorB), Slst
            )

def build_cells(s):
    return [Cell(s) for _ in range(CELLS)]

def step(cells, t_curr, s):
    for idx in range(len(cells)):
        cell = cells[idx]
        cell.grow()
        division = cell.divide(s)
        if division is not None:
            cells.append(division)
                

        
def vsum(cells):
    result = 0
    for cell in cells:
        result += 2 ** (cell.age / cell.lifetime)
    return result

def tau(s, sigmab):
    f = lambda V : 1 / V * (1 + (s + (sigmab - s) / V) ** N) \
    / (Lambda0 + Lambda1 * (s + (sigmab - s) / V) ** N)
    t0, err = integrate.quad(f, 1, 2)
    return t0

    
def build_ylst(s):
    t_curr = 0 
    cells = build_cells(s)
    logvlst = []
    tlst = []
    while len(cells) < MAX_POP:
        step(cells, t_curr, s)
        if t_curr > T_SAMPLE:
            logvlst.append(np.log(vsum(cells)))
            tlst.append(t_curr)
        t_curr += DT
    p = np.polyfit(tlst, logvlst, 1)
    return p[0]


       
###########
class Cell():
    def __init__(self, s, age = None, protein_b = None, 
                 lifetime = None):
        if age is None:
            self.protein_b = s
            self.lifetime = tau(s, s)
            # assume that the volume is distributed uniformly.
            self.age = self.lifetime * np.log(np.random.uniform(low=1.0, high = 2.0)) / np.log(2)
            
        else:
            self.age = age
            self.protein_b = protein_b
            self.lifetime = lifetime

    def grow(self):
        self.age += DT
        
    def divide(self, s):
        if self.age > self.lifetime - DT:
            protein_d = self.protein_b + s
            a = np.tanh(k * protein_d)
            self.protein_b = protein_d * (1 - a) / 2
            protein2 = protein_d * (1 + a) / 2
            self.age = 0
            age2 = 0
            self.lifetime = tau(s, self.protein_b)
            lifetime2 = tau(s, protein2)
            return Cell(s, age2, protein2, lifetime2)
        else:
            return None

if __name__ == '__main__':
    main()

