import sys
from random import uniform
import numba as nb
import numexpr as ne
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy import signal


def f_t2(F0,t,tc):
    return F0/pow(1-t/tc,3./8.)
def A(Mc,F0,t,D,tc) : 
    return 4*pow(Mc,5./3.)*pow(pi*f_t2(F0,t,tc),2./3.)/D
def phi(t,Mc,F0,tc):
    return 10.05*F0*tc*pow(1-t/tc,5./8.)

def generation(M1,M2,d,F_ech): 
    #M1 et M2 en masse solaires
    #D en pc
    #F_ech en Hz

    G = 6.674e-8 #cm3.g-1.s-2
    c = 2.99e10 #cm.s-1
    M_sun = 1.99e33 #g
    pc = 3.086e18 #cm

    temp_t = []
    temp_hc = []
    F = []
    D = d*pc
    m1 = M1*M_sun
    m2 = M2*M_sun
    mtot = m1+m2
    eta = m1*m2/pow(mtot,2.)
    Mc = mtot*pow(eta,3./5.)
    F0= 10.
    tc = 9.23e-4*pow(c,5.)/(pow(F0,8./3.)*pow(Mc,5./3.)*pow(G,5./3.))
    t = 0
    i=0
    window = signal.hann(int(F_ech*0.1*2), sym = False)
    while t<=tc- 1/F_ech: 
        t+=1/F_ech
        temp_t.append(t)
        if t<0.1 : 
            tapering = window[i]
        elif (tc-t)<0.1 : 
            tapering =window[i-int(F_ech*tc)]
        else : 
            tapering = 1.
        temp = tapering*pow(G,5./3.)*A(Mc,F0,t,D,tc)*cos(phi(t,Mc,F0,tc))/pow(c,4.)
        temp_hc.append(temp)
        i+=1
    print("nombre de points : ",i," durée analytique tc (s): ",tc)
    df = pd.DataFrame(temp_hc)
    df.to_csv("template.txt", sep = ' ',index = False, header = None)
    return ( temp_hc)

    
    
    
    
