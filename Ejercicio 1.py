# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 20:29:24 2021

@author: yisus
"""
import numpy as np
import matplotlib.pyplot as plt

interespaciadoH=100
valorInicialY=0
valorFinalY=3000
presionInicial=101325

cantidadDePuntos=int((valorFinalY-valorInicialY)/interespaciadoH+1)

valores_de_alturas=np.linspace(valorInicialY,valorFinalY,cantidadDePuntos)

valores_de_presiones=np.zeros(cantidadDePuntos)
valores_de_presiones[0]=presionInicial



def DefinirFuncion(y,p):
    '''
    Definición de la función f de la ecuación diferencial tal que dp/dy=f(y,p)

    Parameters
    ----------
    y : valor de altura al que se quiere calcular la función f
    p : valor de presión a la que se quiere calcular f

    Returns
    -------
    f : Valor numérico de la función evaluada en (y,p)

    '''
    M=0.0289647
    R=8.314462
    g=9.8
    f=(200*M*g*p)/(R*(-58.6*10**3+y))
    return f
    

def Calcular_RK4(valores_de_alturas,valores_de_presiones,interespaciadoH):
    '''
    Se calcula por el método de Ruggen-Kutta 4 la  EDO con función definida F 
    antes.

    Parameters
    ----------
    valores_de_alturas : Array con las alturas a las que se quiere calcular la
    presión, el primero valor del array debe ser la altura inicial
    valores_de_presiones : Array con ceros donde se van a colocar los valores 
    de presiones calculados, el primer valor del array sebe ser la presión inicial
    interespaciadoH : interespaciado entre puntos denominado h

    Returns
    -------
    list
        Se retorna un array de dos entradas donde en la primera entrada están
        los valores de las alturas y en el segundo las presiones a esas alturas

    '''
    for i in range(1,len(valores_de_alturas)):
        
        k1=interespaciadoH*DefinirFuncion(valores_de_alturas[i-1],valores_de_presiones[i-1])   
        k2=interespaciadoH*DefinirFuncion(valores_de_alturas[i-1]+k1/2,valores_de_presiones[i-1]+interespaciadoH/2)
        k3=interespaciadoH*DefinirFuncion(valores_de_alturas[i-1]+k2/2,valores_de_presiones[i-1]+interespaciadoH/2)
        k4=interespaciadoH*DefinirFuncion(valores_de_alturas[i-1]+k3,valores_de_presiones[i-1]+interespaciadoH)
        valores_de_presiones[i]=valores_de_presiones[i-1]+1/6*(k1+2*k2+2*k3+k4)
        
    return [valores_de_alturas,valores_de_presiones]

valores_de_alturas,valores_de_presiones=Calcular_RK4(valores_de_alturas, valores_de_presiones, interespaciadoH)

fig, ax= plt.subplots()
ax.plot(valores_de_alturas,valores_de_presiones,label='RK4 Programado')
ax.set_xlabel('Alturas (m)')
ax.set_ylabel('Presiones (Pa)')
ax.set_title('Presión atmosférica vs altura sobre el nivel del mar')




'''
Esta sección no pertence al código del ejercicio 1 si no solo se utilizó para 
graficar y se considero pertinerte conservarlo.

Solución Análitica

valoresAnalitico=2.8230*10**(-28)*(58.6*10**3-valores_de_alturas)**(200*0.0289647*9.8/8.314462)




Solución de Elena 
    también presente en el documento ejercicio 2 aquí solo se
    colocó para graficar
    
# -*- coding: utf-8 -*-
"""
Tarea 2 

@author: María Elena Esquivel Murillo
"""

import scipy.integrate as spint
import numpy as np

def DefinirFunción(y,p):
    '''
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



fig, ax= plt.subplots()
ax.plot(valores_de_alturas,valores_de_presiones,label='RK4 Programado')

#descomentar si se desea ver los resultados de RK45 descomentar antes el código
#llamado Solución de Elena
#ax.plot(valores_de_alturas,resultados[0],label='RK45 SciPy')

#descomentar si se quiere graficar la solución analítica
#ax.plot(valores_de_alturas,valoresAnalitico,label='Solución Analítica')

ax.set_xlabel('Alturas (m)')
ax.set_ylabel('Presiones (Pa)')
ax.set_title('Presión atmosférica vs altura sobre el nivel del mar')
ax.legend()
'''



        
        