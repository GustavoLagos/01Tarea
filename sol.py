# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as const
import scipy as sci

#P1: grafico del espectro solar

sun_AMO = np.genfromtxt('sun_AM0.dat', dtype='str')
sun_AMO_num = sun_AMO.astype(np.float)
wavelength = np.log10(10*sun_AMO_num[:,0]) #en angstrom
flujo = np.log10((10**10)*sun_AMO_num[:,1]) #en unidades cgs

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Espectro del Sol')
ax.set_xlabel('log10 [Longitud de onda ($\AA$)]')
ax.set_ylabel('log10 [Flujo ($gs^{-3}cm$)]')

plt.gcf()
plt.ion()
plt.plot(wavelength, flujo)
plt.show()
plt.savefig('espectro_solar.png')

#P2: Luminosidad del sol

flujo_inicial = (10**10)*sun_AMO_num[0,1]
flujo_final = (10**10)*sun_AMO_num[len(sun_AMO_num[:,1])-1,1]
delta_x = (np.amax(10*sun_AMO_num[:,0]) - np.amin(10*sun_AMO_num[:,0]))/(len(sun_AMO_num[:,0])-1) 

i = 1
suma = 0.0
while i < len(sun_AMO_num[:,1]):
    suma += (10**10)*sun_AMO_num[i,1]
    i += 1

integral1 = (delta_x/2.0)*(flujo_inicial + 2.0*suma + flujo_final) #metodo del trapecio

luminosidad = 4*np.pi*integral1*(10**(-12))*(const.au**2)

#P3: Integracion de la funcion de Planck y radio efectivo del Sol con el metodo del trapecio

def f(x):
    return (np.tan(x)**3)/(np.exp(np.tan(x))-1) 

n = 2000
x= np.linspace(0, np.pi/2.0, n+1) #arreglo de 0 a pi/2, de dimension n+1
h = (np.pi/2.0)/(n)

i = 0
a = np.zeros(shape=(n))
while i < n:
    a[i] = f(i+1)
    i+=1

s = 2.0*np.sum(a) + f(n) #f(0) = 0
integral = s*h/2.0

T = 5777 #temperatura efectiva del Sol
P = ((2*np.pi*const.h)/const.c)*((const.k_B*T/const.h)**4)*integral


#P4: Comparacion de scipy con integrales anteriores

I1 = sci.integrate.trapz(flujo, wavelength)
I2 = sci.integrate.quad(lambda x: (np.tan(x)**3)/(np.exp(np.tan(x))-1),0,np.pi/2.0)







































