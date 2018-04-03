#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 18:15:28 2018

@author: luiggi
"""
import numpy as np

class Matrix():
    
    def __init__(self, nvx = None):
        self.__N = nvx - 2 
        self.__A = np.eye(self.__N)

    def __del__(self):
        del(self.__N)
        del(self.__A)
        
    def mat(self):
        return self.__A
    
    def build(self, coefficients = None):
# nx = 5, nvx = 6
# 0     1     2     3     4     5  <-- Volumes 
# o--|--x--|--x--|--x--|--x--|--o
#       0     1     2     3        <-- Unknowns    
#
#        0   1   2   3
# --+---------------------
# 0 | [[12. -4.  0.  0.]
# 1 | [ -4.  8. -4.  0.]
# 2 | [  0. -4.  8. -4.]
# 3 | [  0.  0. -4. 12.]]

        aP = coefficients.aP()
        aE = coefficients.aE()
        aW = coefficients.aW()
        aEE = coefficients.aEE()
        aWW = coefficients.aWW()
        A = self.__A
        A[0][0] = aP[1]
        A[0][1] = -aE[1]
        A[0][2] = -aEE[1]
        
        A[1][0] = -aW[2]        
        A[1][1] = aP[2]
        A[1][2] = -aE[2]
        A[1][3] = -aEE[2]
        for i in range(2,self.__N-2): # range(1,N-3)  <-- (1,2)
            A[i][i] = aP[i+1]
            A[i][i+1] = -aE[i+1]
            A[i][i-1] = -aW[i+1]
            A[i][i+2] = -aEE[i+1]
            A[i][i-2] = -aWW[i+1]
        A[-1][-1] = aP[-2]
        A[-1][-2] = -aW[-2]
        A[-1][-3] = -aWW[-2]
        
        A[-2][-1] = -aE[-3]        
        A[-2][-2] = aP[-3]
        A[-2][-3] = -aW[-3]
        A[-2][-4] = -aWW[-3]

if __name__ == '__main__':

    a = Matrix(6)
    print('-' * 20)  
    print(a.mat())
    print('-' * 20)  
    
    from Diffusion import Diffusion1D        
    df1 = Diffusion1D(6, 1, 0.25)
    df1.alloc(6)
    df1.calcCoef()
    df1.setSu(100)
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  

    df1.bcDirichlet('LEFT_WALL', 2)
    df1.bcDirichlet('RIGHT_WALL', 1)
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  
    
    a.build(df1)
    print(a.mat())
    print('-' * 20)  
