#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class Coefficients():
    """
    Esta clase define los arreglos principales para los coeficientes del
    metodo de Volumen Finito. Los arreglos son definidos como variables de
    clase para que sean compartidos por todos los objetos de esta clase.
    
        Métodos:
        constructor(nvx,delta): inicia los atributos nvx y delta
        destructor(): delete atributes
        alloc(n): asigna arreglos con ceros a los atributos de coeficientes (aP, aE, etc.) con
                  el objetivo de reservar memoria
        setVolumes(nvx): set atributo nvx
        setDelta(delta): set atributo delta
        aP():get aP
        aW(): get aW
        aWW(): get aWW
        aE(): get aE
        aEE(): get aEE
        sU(): get sU                
        bcDirichlet(wall,phi): ajusta los coeficientes de la frontera 'wall' (puede ser 'LEFT_WALL' o 'RIGHT_WALL')
                               de acuerdo a la condición de frontera con valor 'phi'
        bcNeumman(wall,flux): ajusta los coeficientes de la frontera 'wall' (puede ser 'LEFT_WALL' o 'RIGHT_WALL')
                               de acuerdo a la condición de frontera con flujo de valor 'flux'
        setsU(q): set atributo Su
        setSp(Sp): set atributo Sp

        
    Atributos:
        aP: Coeficiente central
        aE: Coeficiente de nodo siguiente
        aW: Coeficiente de nodo anterior
        aEE: Coeficiente del segundo nodo siguiente
        aWW: Coeficiente de segundo nodo anterior
        Gamma: Coeficiente Difussivo
        delta: tamaño del los volumenes
        u: velocidad
        sU: coeficientes que representan los términos independientes del sistema de ecuaciones a
               resolver y que son consecuencia de las fuentes
        sP: Corrección en los coeficientes aP debido a fuentes.
        
    """    
    
    #atributos
    __aP = None
    __aE = None
    __aW = None
    __aEE = None
    __aWW = None
    __Su = None
    __nvx = None
    __delta = None

    def __init__(self, nvx = None, delta = None):
        Coefficients.__nvx = nvx
        Coefficients.__delta = delta


    @staticmethod
    def alloc(n):
        if Coefficients.__nvx:
            nvx = Coefficients.__nvx
        else:
            nvx = n
        Coefficients.__aP = np.zeros(nvx)
        Coefficients.__aE = np.zeros(nvx)
        Coefficients.__aW = np.zeros(nvx)
        Coefficients.__aEE = np.zeros(nvx)
        Coefficients.__aWW = np.zeros(nvx)
        Coefficients.__Su = np.zeros(nvx)        
    
    def setVolumes(self, nvx):
        Coefficients.__nvx = nvx
        
    def setDelta(self, delta):
        Coefficients.__delta = delta
        
    def aP(self):
        return Coefficients.__aP

    def aE(self):
        return Coefficients.__aE
    
    def aW(self):
        return Coefficients.__aW
    
    def aEE(self):
        return Coefficients.__aEE
    
    def aWW(self):
        return Coefficients.__aWW
    
    def Su(self):
        return Coefficients.__Su

    @staticmethod
    def bcDirichlet(wall, phi):
        aP = Coefficients.__aP
        aE = Coefficients.__aE
        aW = Coefficients.__aW
        aEE = Coefficients.__aEE
        aWW = Coefficients.__aWW
        Su = Coefficients.__Su

        if wall == 'LEFT_WALL':
            aP[1] += aW[1] + 3*aWW[1]
            Su[1] += (2 *aW[1] + 4*aWW[1]) * phi
            aW[2] -= aWW[2] #condición del segundo nodo (requerida en métodos de orden mayor a 1)
            Su[2] += (2 *aWW[2] ) * phi
        elif wall == 'RIGHT_WALL':
            aP[-2] += aE[-2] + 3*aEE[1]
            Su[-2] += (2 *aE[-2] + 4*aEE[1] )* phi
            aE[2] -= aEE[2] #condición del penúltimo nodo (requerida en métodos de orden mayor a 1)
            Su[-3] += (2*aEE[1] )* phi
            
#        if wall == 'LEFT_WALL':
#            Su[1] += (0.25*max((rho*u[i-1],0))*1)
#            aP[1] += 2*aW[1] + 9*aWW[1]
#            aE[1] += (1./3.)*aW[1] + 2*aWW[1]
#            Su[1] += ( (8./3.)*aW[1] + 8*aWW[1]) * phi
#            aP[2] -= (1./3.)*aWW[2] #condición del segundo nodo (requerida en métodos de orden mayor a 1)
#            aW[2] -= 2*aWW[2]
#            Su[2] += ( (8./3.)*aWW[2] ) * phi
#        elif wall == 'RIGHT_WALL':
#            aP[1] += 2*aE[1] + 9*aEE[1]
#            aW[1] += (1/3.)*aE[1] + 2*aEE[1]
#            Su[1] += ( (8/3.)*aE[1] + 8*aEE[1]) * phi
#            aP[2] -= (1/3.)*aEE[2] #condición del segundo nodo (requerida en métodos de orden mayor a 1)
#            aE[2] -= 2*aEE[2]
#            Su[2] += ( (8/3.)*aEE[2] ) * phi
            
    @staticmethod
    def bcNeumman(wall, flux):
        aP = Coefficients.__aP
        aE = Coefficients.__aE
        aW = Coefficients.__aW
        Su = Coefficients.__Su
        dx = Coefficients.__delta

        if wall == 'LEFT_WALL':
            aP[1] -= aW[1]
            Su[1] -= aW[1] * flux * dx
        elif wall == 'RIGHT_WALL':
            aP[-2] -= aE[-2]
            Su[-2] += aE[-2] * flux * dx  
            
    def setSu(self, q):
        Su = Coefficients.__Su
        dx = Coefficients.__delta
        Su += q * dx
        
    def setSp(self, Sp):
        aP = Coefficients.__aP
        dx = Coefficients.__delta
        aP -= Sp * dx
        

if __name__ == '__main__':
    
    coef1 = Coefficients(6, 0.25)
    coef1.alloc(6)
    coef1.setSu(100)
    coef1.setSp(-2)
    
    print('-' * 20)  
    print(coef1.aP(), coef1.aE(), coef1.aW(), coef1.Su(), sep = '\n')
    print('-' * 20)  

    ap = coef1.aP()
    ap[2] = 25
    print(ap, coef1.aP(),sep='\n')
    print('-' * 20)  

    ae = coef1.aE()
    aw = coef1.aW()
    su = coef1.Su()
    ae.fill(5)
    aw.fill(5)
    ap.fill(10)
    coef1.setSp(-2)
    coef1.bcDirichlet('LEFT_WALL', 2)
    coef1.bcNeumman('RIGHT_WALL', 1)
    print(coef1.aP(), coef1.aE(), coef1.aW(), coef1.Su(), sep = '\n')
    print('-' * 20)  

