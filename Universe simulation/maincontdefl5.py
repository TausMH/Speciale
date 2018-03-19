# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 21:04:59 2017
 
@author: Taus
"""
 
"""
Assumptions
 
Homogeneity meaning constant density at rho = 10**(-5) / Mpc**3
 
Isotropic emission and constant luminosity, L(z) = L0 & Lcr = Lnu
Furthermore, n_pix(ne) is proportional to ne^-alpha.
                                             
"""
 
import numpy as np
import numpy.random as npr
import scipy as sp
import scipy.optimize
import scipy.stats
import healpy as hp
import matplotlib.pyplot as plt
import time
 
# Setup
# Primary variables
rho = 10**(-7)
Lcr = 1. * 10**(0) #MUST be a float
Lnu = Lcr
maxtime = 10 # in years
nside = 2**4 #resolution for skymap
rmaxCR = 400.
Ecr = 4.5*10**19
Lambda = 10**(-3.5)
B = 10**(-12)
 
# Secondary variables
Hp = 14.4 * 10**3 #Assuming expanding universe, in Mpc
Hp = 14.4 * 10**3 / 10.#Assuming expanding universe, in Mpc
#nestbool = True #skymap style
nestbool = False #skymap style
sortlim = 5000
Lnormbool = True
LnormvalCR = 10000
Lnormvalnu = 10000


maxtime = maxtime * 3600 * 24 * 365.25
 
 
def GenUniverse(rho, Hp):
	nsource = int(4/3 * Hp**3 * np.pi * rho)

	# For reference: https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
	# Generating radius distribution:
	r = Hp * npr.rand(nsource)**(1/3.)

	# Generating angular distribution
	z = npr.uniform(-1,1,size=(nsource))
	theta = 2*np.pi * npr.rand(nsource)
	x = np.multiply((1 - z**2)**(1/2.), np.cos(theta))
	y = np.multiply((1 - z**2)**(1/2.), np.sin(theta))

	xf = np.multiply(r, x)
	yf = np.multiply(r, y)
	zf = np.multiply(r, z)

	print "Universe generated"
	return np.array([xf, yf, zf])
 
def Numnu(coords, Lnu, maxtime, Lnormbool, Lnormval):
	norm = np.linalg.norm(coords,axis=0)
	F = Lnu / (4 * np.pi * norm**2)
	mF = F * maxtime
	if Lnormbool==1:
		mF = float(Lnormval) / sum(mF) * mF
	print "Neutrino flux generated"
	return np.random.poisson(mF)
 
 
def Numcr(coords, Lcr, maxtime, rmax, Lnormbool, Lnormval):
	norm = np.linalg.norm(coords,axis=0)
	mask = (norm < rmax)*1
	print "Number of CR sources: ", sum(mask)
	F = Lcr / (4 * np.pi * norm**2)
	mF = F*mask * maxtime
	if Lnormbool==1:
		mF = float(Lnormval) / sum(mF) * mF
	print "Cosmic ray flux generated"
	return np.random.poisson(mF)

def MagScatteringNear(Ncr, coords, sortlim, Lambda, B, E):
	nonzeromask =  Ncr >= sortlim
	mNcr = (Ncr[nonzeromask]).astype(int)
	mCoords = coords[:, nonzeromask]
	Ncrlen = len(mNcr)
	psum = sum(mNcr)
	sctcoords = np.zeros([3,psum])
	
	[thetaarr, phiarr] = hp.pixelfunc.vec2ang(np.transpose(mCoords)) #Might have to be (N,3) instead of (3,N)
	D = np.linalg.norm(mCoords ,axis=0)
	sigma = 0.025 * (D/Lambda)**(1/2.) * (B/10**(-11)) / (E/10**20) * np.pi / 180.	
	
	count = 0
	
	for i in xrange(Ncrlen):

		x = npr.uniform(0, 1, mNcr[i])
		t = npr.uniform(0, 2*np.pi, mNcr[i])
		cosalpha = 1 + (sigma[i])**2 * (np.log(1 - x * (1 - np.exp(-2/(sigma[i])**2))))
		
		alpha = np.arccos(cosalpha)
		theta = thetaarr[i]
		phi = phiarr[i]
		
		c = np.array([[np.cos(phi)*np.cos(theta)*np.sin(alpha)*np.cos(t) - np.sin(phi)*np.sin(alpha)*np.sin(t) + np.cos(phi)*np.sin(theta)*np.cos(alpha)],
                      [np.sin(phi)*np.cos(theta)*np.sin(alpha)*np.cos(t) + np.cos(phi)*np.sin(alpha)*np.sin(t) + np.sin(phi)*np.sin(theta)*np.cos(alpha)], 
                      [-np.sin(theta)*np.sin(alpha)*np.cos(t) + np.cos(theta)*np.cos(alpha)]])
		
		sctcoords[:,count:count+mNcr[i]] = c[:,0,:]
		
		count = count + mNcr[i]

	return sctcoords


def MagScatteringFar(Ncr, coords, sortlim, Lambda, B, E):
	nonzeromask = (sortlim > Ncr) * (Ncr > 0)
	mNcr = (Ncr[nonzeromask]).astype(int)
	mCoords = coords[:, nonzeromask]

	[theta, phi] = hp.pixelfunc.vec2ang(np.transpose(mCoords)) #HAS to be (n,3)
	
	D = np.linalg.norm(mCoords ,axis=0)
	sigma = 0.025 * (D/Lambda)**(1/2.) * (B/10**(-11)) / (E/10**20) * np.pi / 180.
	sigma = np.transpose(np.array([sigma]))
	theta = np.transpose(np.array([theta]))
	phi   = np.transpose(np.array([phi]))

	Ncrlen = len(mNcr)
	maxlen = max(mNcr)
	x = npr.uniform(0,1,[Ncrlen,maxlen])
	t = npr.uniform(0,2*np.pi,[Ncrlen,maxlen])
	cosalpha = 1 + sigma**2 * (np.log(1 - x * (1 - np.exp(-2/sigma**2))))
	alpha = np.arccos(cosalpha)
	
	c = np.array([[np.cos(phi)*np.cos(theta)*np.sin(alpha)*np.cos(t) - np.sin(phi)*np.sin(alpha)*np.sin(t) + np.cos(phi)*np.sin(theta)*np.cos(alpha)],
				  [np.sin(phi)*np.cos(theta)*np.sin(alpha)*np.cos(t) + np.cos(phi)*np.sin(alpha)*np.sin(t) + np.sin(phi)*np.sin(theta)*np.cos(alpha)], 
				  [-np.sin(theta)*np.sin(alpha)*np.cos(t) + np.cos(theta)*np.cos(alpha)]])
	c = c[:,0,:,:]
	
	#Create maskarr
	binmask = np.zeros([Ncrlen,maxlen])
	for i in xrange(Ncrlen):
		binmask[i, 0:mNcr[i]] = 1
	
	pc  = np.multiply(binmask, c)
	cnorm = np.linalg.norm(pc ,axis=0)
	mcnorm = cnorm != 0
	
	sctcoords  = pc[:,mcnorm]
	
	return sctcoords
	
def MagScattering(Ncr, coords, sortlim, Lambda, B, E):
	
	nearstart = time.time()
	coordsNear = MagScatteringNear(Ncr, coords, sortlim, Lambda, B, E)
	nearend = time.time()
	coordsFar = MagScatteringFar(Ncr, coords, sortlim, Lambda, B, E)
	farend = time.time()
	
	print "Time spent on near-scattering:", nearend - nearstart
	print "Time spent on far-scattering:",  farend  - nearend
	
	appended = np.append(coordsNear, coordsFar, axis=1)
	
	return appended

def Likelihood(ns, S, NCR, B):
	

	LH = ns/NCR * S + (NCR-ns)/NCR * B
	LH = 2*np.sum(np.log(LH/B), axis=0) # Same as LH/LH0

	# Negative in order to utilize minimize.
	return -LH
	

def MaxLikelihood(eta, cr4plot, crcoords, nside, nestbool):
	
	print "maximum start"
	
	npix = hp.nside2npix(nside)
	NCR  = crcoords.shape[1]
	pixcoords = np.asarray(hp.pix2vec(nside, np.arange(npix))) #np.repeat(np.arange(cr4plot.size))
	
	angdiff = np.arccos(np.dot(np.transpose(crcoords), pixcoords)) #Has to be shape n,3 x 3,m
	S = np.exp(-abs(angdiff / (2*np.pi))*50) #Arbitrary x50 factor to make maxima appear
	B = np.zeros([NCR,npix]) + 1./npix
	
	x0 = np.zeros(npix) + 2
	solutionsx = np.zeros(npix)
	solutionsf = np.zeros(npix)
	#bnds = ((0,None),(0,None)) * (npix/2) #Creates tuples for bounds making values positive
	bnds = ((0,10000),)
			   
	
	fl1 = time.time()
	for i in xrange(npix):
		res = sp.optimize.minimize(Likelihood, x0[i], args = (S[:,i], NCR, B[:,i]), bounds = bnds)
		solutionsx[i] = res.x
		solutionsf[i] = -res.fun #Since negative from minimization
	fl2 = time.time()
	print "Time elapsed for for loop:", fl2-fl1
	
	"""
	lc1 = time.time()
	lstcomp = [sp.optimize.minimize(Likelihood, x0[i], args = (S[:,i], NCR, B[:,i]), bounds = bnds) for i in xrange(npix)]
	#TSarray1, nsarray1 = zip(*lstcomp)
	lc2 = time.time()
	print "Time elapsed for list comprehension:", lc2-lc1
	print "Difference between methods:", (fl2-fl1) - (lc2-lc1)
	"""
	
	return solutionsx, solutionsf

#Delta based maxlikelihood

def MaxTS(N, Nin, B):
	Nout = N - Nin
	nsmax = (N * (Nin * (B - 1) + Nout * B)) / ((B - 1) * (Nin + Nout))
	TS = 2 * (Nin*np.log(nsmax/(B*N) + (1 - nsmax/N)) + Nout*np.log(1-nsmax/N))
	return TS, nsmax
							   
def MaxLikelihoodDelta(nside, pix):
	
	print "max likelihooding"
	
	npix = hp.nside2npix(nside)
	B = 1./npix
	N = np.sum(pix)
	
	TSarray = np.array([])
	nsarray = np.array([])
	
	for i in xrange(npix):
		Nin = pix[i]
		(TS,nsmax) = MaxTS(N, Nin, B)
		TSarray = np.append(TSarray,TS)
		nsarray = np.append(nsarray,nsmax)
	
	return TSarray, nsarray


ScriptStart = time.time()
# Generating data
universe = GenUniverse(rho, Hp)
Ncr = Numcr(universe, Lcr, maxtime, rmaxCR, Lnormbool, LnormvalCR)
Nnu = Numnu(universe, Lnu, maxtime, 		Lnormbool, Lnormvalnu)


# New scatter
ScatterStart2 = time.time()
sctCoords = MagScattering(Ncr, universe, sortlim, Lambda, B, Ecr)
ScatterEnd2 = time.time()


print "Number of CR events:", sctCoords[0,:].size

# Generating skymap using healpix
npix = hp.nside2npix(nside)
cr4plot = np.zeros(npix)
nu4plot = np.zeros(npix)

scnu4plot = hp.vec2pix(nside, universe[0,:], universe[1,:], universe[2,:], nest=nestbool)
sccr4plot = hp.vec2pix(nside, sctCoords[0,:], sctCoords[1,:], sctCoords[2,:], nest=nestbool)
sccr4plotundef = hp.vec2pix(nside, universe[0,:],universe[1,:],universe[2,:], nest=nestbool)

#Method 1
cr4plot = cr4plot + np.bincount(sccr4plot, minlength = npix)
nu4plot = nu4plot + np.bincount(scnu4plot, weights = Nnu, minlength = npix)
cr4plotundef = np.bincount(sccr4plotundef, weights = Ncr, minlength = npix) + 1

background = np.zeros(npix) + 1

cr4plot = cr4plot + background
nu4plot = nu4plot + background



# Plot skymaps

#hp.mollview(cr4plot, nest=nestbool, title="Cosmic ray nested distribution", norm='log')
#hp.graticule()
"""
hp.mollview(cr4plot, coord=['E','G'], nest=nestbool, title="Cosmic ray nested distribution", norm='log')
hp.graticule(coord=['E'])

hp.mollview(cr4plotundef, coord=['E','G'], nest=nestbool, title="Cosmic ray no deflection", norm='log')
hp.graticule(coord=['E'])
hp.mollview(cr4plotundef, nest=nestbool, title="Cosmic ray no deflection", norm='log')
hp.graticule()
"""
#hp.mollview(nu4plot, nest=nestbool, title="Neutrino nested distribution", norm='log')
#hp.graticule()


"""
# Plot power spectra and cross spectrum
pwrCR = hp.anafast(cr4plot)
pwrnu = hp.anafast(nu4plot)
crs   = hp.anafast(cr4plot,nu4plot)

ellCR = np.arange(len(pwrCR))
ellnu = np.arange(len(pwrnu))
ellcrs = np.arange(len(crs))

plt.figure()
plt.plot(ellCR[1:], pwrCR[1:])
plt.xlabel('ell')
plt.ylabel('c1')
plt.grid()
plt.title('CR power spectrum')

plt.figure()
plt.plot(ellnu[1:], pwrnu[1:])
plt.xlabel('ell')
plt.ylabel('c1')
plt.grid()
plt.title('Neutrino power spectrum')

plt.figure()
plt.plot(ellcrs[1:], crs[1:])
plt.xlabel('ell')
plt.ylabel('c1')
plt.grid()
plt.title('Cross spectrum')
"""

# Plot histogram of npix with ne
"""
plt.figure()
maxx = max(cr4plot[cr4plot != 1])
plt.hist(cr4plot[cr4plot != 1], bins = np.logspace(0.1,np.log10(maxx),50))
plt.gca().set_xscale("log",nonposx='clip')
plt.gca().set_yscale("log",nonposy='clip')
plt.xlabel('n_event')
plt.ylabel('n_pix')
plt.title('Cosmic rays')

plt.figure()
maxx = max(nu4plot[nu4plot != 1])
plt.hist(nu4plot[nu4plot != 1], bins = np.logspace(0.1,np.log10(maxx),50))
plt.gca().set_xscale("log",nonposx='clip')
plt.gca().set_yscale("log",nonposy='clip')
plt.xlabel('n_event')
plt.ylabel('n_pix')
plt.title('Neutrinos')
"""


# Likelihood testing
eta = 0.2
nsmax, TSmax = MaxLikelihood(eta, cr4plot, sctCoords, nside, nestbool)
hp.mollview(nsmax, nest=nestbool, title="Cosmic ray ns map") #, norm='log'
hp.graticule()
hp.mollview(TSmax, nest=nestbool, title="Cosmic ray TS map")
hp.graticule()

#Delta Likelihood testing
"""
TSarrayCR, nsarrayCR = MaxLikelihoodDelta(nside, cr4plot)
TSarraynu, nsarraynu = MaxLikelihoodDelta(nside, nu4plot)
TSarrayCRundef, nsarrayCRundef = MaxLikelihoodDelta(nside, cr4plotundef)

TSarrs = (TSarrayCR, TSarraynu, TSarrayCRundef)
plottitles = ('TS histogram for deflected CRs with scaled chi-squared dist','TS histogram for neutrinos with scaled chi-squared dist','TS histogram for undeflected CRs with scaled chi-squared dist')

for i in xrange(len(TSarrs)):		
	plt.figure()
	hst = plt.hist(TSarrs[i], bins=50)
	#histmax = max(hst[0])
	#x = np.linspace(0,100,400)
	#plt.plot(x, sp.stats.chi2.pdf(x,1)*histmax)
	plt.gca().set_yscale("log",nonposy='clip')
	plt.xlabel('TS')
	plt.ylabel('count')
	plt.title(plottitles[i])

#hp.mollview(TSarrayCR, nest=nestbool, title="Cosmic ray TS map", norm='log')
#hp.graticule()
"""

# Plotting the 3D circle
"""
from mpl_toolkits.mplot3d import Axes3D

stepsize = 100
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(sctCoords[0,0::stepsize], sctCoords[1,0::stepsize], sctCoords[2,0::stepsize], s=20, depthshade=True)
plt.xlim(-1,1)
plt.ylim(-1,1)
ax.set_zlim(-1,1)
"""

ScriptEnd = time.time()

print "Total time elapsed:", ScriptEnd-ScriptStart 
print "Time spent on new scattering:", ScatterEnd2-ScatterStart2

plt.show()