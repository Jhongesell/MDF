# -*- coding: utf-8 -*-
"""
@author: omar jonathan mendoza bernal

 Upwind scheme for
 u_t + a(x) u_x = 0

 Ejemplo
 u_t +  a u_x = 0
 Si u(x,0) =  f(x) es la condición inicial
 la solución analitica es
 u(x,t) = f(x-at
"""

from math import exp
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt

# Ecuación
# u_t + 2 u_x = 0
coef_a = 2

def condicion_inicial (x):
    """
    Condición inicial de la ecuación
    """
    y = exp (-(x - 0.2)*(x - 0.02))
    return y

def sol_analitica (x, t):
    """
    Solución analítica de la ecuación
    """
    y = exp(-(x-2*t)*(x-2*t))

    return y

############
# Dominio
############
a = -2.0
b = 8.0

# Partición del dominio
nNodos = 401
h = (b - a)/(nNodos - 1)

# Intervalo de tiempo
dt = 0.012

# creo el vector w donde se guradará la solución para cada tiempo
# B matriz del lado derecho
w = np.zeros((nNodos,1))
B = np.zeros((nNodos, nNodos))
B = np.matrix(B)
espacio = np.zeros((nNodos,1))

for i in xrange(nNodos):
    xx_ = a + i*h
    espacio[i] = xx_
    w[i] = condicion_inicial(xx_)

print "Espacio"
print espacio

print "Condición Inicial"
print w

mu = coef_a * dt / h

if mu <= 1:
    print "mu ", mu
    print "Buena aproximación"
else:
    print "mu ", mu
    print "Mala aproximación"

if coef_a >= 0:
    B[0,0] = 1 - mu
    for i in xrange (1, nNodos):
        B[i,i-1] = mu
        B[i,i] = 1 - mu
else:
    B[0,0] = 1 + mu;
    B[0,1] = -mu;

    # se completa las matrices desde el renglon 2 hasta el m-2
    for i in xrange(1, nNodos-1):
        B[i,i] =  1 + mu;
        B[i, i+1] = -mu;

    B[nNodos-1,nNodos-1] = 1 + mu

# para guardar la soulción analitica
xx = np.zeros((nNodos,1));

iteraciones = 201

# Matriz sparse csr
Bs = sparse.csr_matrix(B);

# Resolvemos el sistema iterativamente
for j in xrange (1, iteraciones):
    t = j*dt;
    w = Bs*w;

	# Imprimir cada 20 iteraciones
    if  j%20 == 0:
        print "t", t
        print "w", w
        #for k_ in xrange (nNodos):
         #   xx[k_] = sol_analitica(espacio[k_], t)
          #  print xx

        plt.plot(espacio, w)
        #plt.plot(espacio, xx)

#    print "XX"
#    print xx
#    print "W"
#    print w