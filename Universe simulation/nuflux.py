import numpy as np
import numpy.random
import scipy as sp
import scipy.integrate
import matplotlib.pyplot as plt
import astropy.cosmology
from astropy.cosmology import WMAP9 as wmap9
from astropy import units as u

def rho(z, rho0):
	
	if z < 1:
		ni = 3.4
	elif 1 <= z < 4:
		ni = -0.3
	else:
		ni = -3.5
	return rho0 * (1 + z)**ni

def Q(Nnu, rho0, d):
	
	
	fun1 = lambda z: rho(z, rho0) / wmap9.H(z).to('1/s').value
	fun2 = lambda z: rho(z, rho0)/rho(0, rho0) / (wmap9.H(z).value/wmap9.H(0).value)
	fun3 = lambda z: rho(z, rho0)/rho(0, rho0) / (wmap9.H(z).value/wmap9.H(0).value) * (1+z)**(-2)
	
	xi = [sp.integrate.quad(fun1, 0, np.inf)[0],
		  sp.integrate.quad(fun2, 0, np.inf)[0],
		  sp.integrate.quad(fun3, 0, np.inf)[0]]
	xi = (xi * u.s * (u.Mpc)**(-3)).to('s/cm3').value
	
	c = 29979245800 #In cm/s
	Aeff = (1 * (u.km)**2).to('cm2').value
	Tlive = 5 # in years
	Tlive = Tlive * 3600 * 24 * 365.25
	
	Phi = Nnu / (Aeff * Tlive * (4*np.pi))
	Q0 = Phi * 4*np.pi / (c * xi)
	
	F = Q0 / (4*np.pi*d**2)
	
	return F

print Q(1000, 10**(-5), (100 * u.Mpc.to('cm')))

"""

def Jnu(E, H0, d):
	
	H0m5 = H0 / (10**(-5))
	xi24 = 2.4
	d1 = d / 10 # Must be given in Mpc
	
	return 0.9 * 10**(-12) / (xi24 * H0m5 * d1**2) / E**2

E = np.logspace(0,18,19)
H0 = 10**(-5)
d = 50

plt.loglog(E, Jnu(E, H0, d))
plt.show()

H0 = np.logspace(-7,-4, 4)
t = 5 * 3600 * 24 * 365.25
plt.loglog(H0, Jnu(1, H0, d)*t)
plt.show()
"""