#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Coefficients import Coefficients

class Diffusion1D(Coefficients):
    """
    Clase que se encarga de calcular los coeficientes difusivos y actualizar los coeficientes generales
    (aP, aW,aE). Esta clase hereda de la clase Coefficients.
    
    Métodos:
        constructor(nvx,Gamma,dx): inicia los atributos nvx y dx de acuerdo a la clase padre además del
                                   atributo Gamma
        destructor(): delete atributes
        calcCoef(): Calcula los coeficientes difusivos y actualiza los coficientes generales (aP, aW,aE).

        
    Atributos:
        aP: Coeficiente central
        aE: Coeficiente de nodo siguiente
        aW: Coeficiente de nodo anterior
        Gamma: Coeficiente Difussivo
        dx: tamaño del los volumenes
        
    """
    
    def __init__(self, nvx = None, Gamma = None, dx = None):
        super().__init__(nvx, dx)
        self.__nvx = nvx
        self.__Gamma = Gamma
        self.__dx = dx

    def __del__(self):
        del(self.__Gamma)
        del(self.__dx)
    
    def calcCoef(self):
        #obtiene los coeficientes usando los geters
        aE = self.aE()
        aW = self.aW()
        aP = self.aP()
        
        #actualiza los coeficientes (Si gamma no es constante se deben descomentar las siguientes líneas y comentar estas)
        aE += self.__Gamma / self.__dx
        aW += self.__Gamma / self.__dx
        aP += aE + aW
 
#        for i in range(self.__nvx):
#            aE[i] += self.__Gamma / self.__dx
#            aW[i] += self.__Gamma / self.__dx
#            aP[i] += aE[i] + aW[i]

if __name__ == '__main__':
    
    df1 = Diffusion1D(5, 5, 1)
    df1.alloc(5)
    df1.calcCoef()
    df1.setSu(100)

    print('-' * 20)  
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  

    df1.bcDirichlet('LEFT_WALL', 2)
    df1.bcDirichlet('RIGHT_WALL', 1)
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  
