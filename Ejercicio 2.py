# -*- coding: utf-8 -*-
"""
Tarea 2 

@author: María Elena Esquivel Murillo
"""

import scipy.integrate as spint
import numpy as np

def DefinirFunción(y,p):
    '''
    Permiten obtener los valores de la función f=dp(y)/dy a partir de los
    cuales se desean encontrar los valores de p(y)

    Parámetros
    ----------
    y : altura, variable independiente.
    p : presión atmosférica, variable dependiente.

    Devuelve
    -------
    valorFunción : valor de la derivada dp(y)/dy en un y específico

    '''
    M=0.0289647
    R=8.314462
    g=9.8
    valorFunción=(200*M*g*p)/(R*(-58.6*10**3+y))
    return valorFunción

y0=0
p0=101325
y=3000
separaciónAlturas=100
cantidadPuntos=int((y-y0)/separaciónAlturas+1)
alturas=np.linspace(y0,y,cantidadPuntos)

métodoRK='RK45'
valoresPresión=spint.solve_ivp(DefinirFunción,[y0,y],[p0],method=métodoRK,t_eval=alturas)
resultados=valoresPresión.y
print(resultados)
