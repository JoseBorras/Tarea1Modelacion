#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 18:46:43 2018

@author: luiggi
"""

import numpy as np
from Coefficients import Coefficients

class Advection1D(Coefficients):
    
    def __init__(self, nvx = None, rho = None, dx = None):
        super().__init__(nvx)
        self.__nvx = nvx
        self.__rho = rho
        self.__dx = dx
        self.__u = np.zeros(nvx-1)

    def __del__(self):
        del(self.__nvx)
        del(self.__rho)
        del(self.__dx)
        del(self.__u)

    def setU(self, u):
        if type(u) == float:
            self.__u.fill(u)
        else:
            self.__u = u

    def u(self):
        return self.__u
    
    def calcCoef(self,metodo='Upwind2'):
        aE = self.aE()
        aW = self.aW()
        aEE = self.aEE()
        aWW = self.aWW()
        aP = self.aP()
        u = self.__u
        rho = self.__rho

        for i in range(1,self.__nvx-1):
            if metodo=='DifCentrales':
                # Diferencias Centrales
                CE = - rho * u[i] * 0.5
                CW =   rho * u[i-1] * 0.5 
                aE[i] += CE
                aW[i] += CW
                aP[i] += CE + CW + rho * (u[i] - u[i - 1])
            elif metodo == 'Upwind1':
                # ------------- Upwind-------------------
                CE = max((-rho*u[i],0)) 
                CW = max((rho*u[i-1],0))
                aE[i] += CE 
                aW[i] += CW
                aP[i] += CE + CW + rho * (u[i] - u[i-1])
            elif metodo == 'Upwind2':                
                #-------------Second Order Upwind -----------------
                CE = 1.5*max((-rho*u[i],0)) + 0.5*max((-rho*u[i-1],0))
                CW = 1.5*max((rho*u[i-1],0)) + 0.5*max((rho*u[i],0))
                CEE = -0.5*max((-rho*u[i],0))     
                CWW = -0.5*max((rho*u[i-1],0))
                aE[i] += CE 
                aW[i] += CW
                aEE[i] += CEE 
                aWW[i] += CWW
                aP[i] += CE + CW+ CEE + CWW + rho * (u[i] - u[i-1]) 
            elif metodo == 'Quick':
                #----------------QUICK---------------------
                CE = (-3./8.)*max((rho*u[i],0)) + (6./8.)*max((-rho*u[i],0)) + (1./8.)*max((-rho*u[i-1],0))
                CW = (1./8.)*max((rho*u[i],0)) + (6./8.)*max((rho*u[i-1],0)) + (-3./8.)*max((-rho*u[i-1],0))
                CEE = -(1./8.)*max((-rho*u[i],0))     
                CWW = -(1./8.)*max((rho*u[i-1],0))
                aE[i] += CE 
                aW[i] += CW
                aEE[i] += CEE 
                aWW[i] += CWW
                aP[i] += CE + CW+ CEE + CWW + rho * (u[i] - u[i-1])                            

if __name__ == '__main__':
    
    nx = 5
    u = np.sin(np.linspace(0,1,nx))
#    u = np.ones(nx)
    print('-' * 20)  
    print(u)
    print('-' * 20)  

    af1 = Advection1D(6, 1, 1)
    af1.alloc(6)
    af1.setU(u)
    print(af1.u())
    print('-' * 20)  

    af1.calcCoef()
    print(af1.aP(), af1.aE(), af1.aW(), af1.Su(), sep = '\n')
    print('-' * 20)  

    af1.bcDirichlet('LEFT_WALL', 2)
    af1.bcDirichlet('RIGHT_WALL', 1)
    print(af1.aP(), af1.aE(), af1.aW(), af1.Su(), sep = '\n')
    print('-' * 20)  



