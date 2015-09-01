#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
import dispar3_sub as dispar 
import idlsave


class saved_grid(object):

	'''
	class which loads then stores variables from a pickle 
	save file. Similar funcitonality to readgrid in IDL version.
	'''

	def __init__(self, fname=None):

		# default is the user doesn't know about the object save file
		import pickle
		if fname == None:
			import os
			fname = os.path.dirname(__file__)	# should reside in same directory as main module
			fname += "objs.pickle"				# actual fname

		# read a number of variables from the module
		with open(fname, 'r') as f:
			self.iconv, self.mlab, self.inter, self.ml, self.pol, self.wls, self.qlab, self.tl, self.tlab, self.ql = pickle.load(f)


class agnspec_model(object):

	'''
    Create an Agnspec model
    
    inputs:
        i) basic physical parameters:

		mass - mass of the black hole (in the solar masses)
		       Default: MASS=1.e9
		mdot - mass accretion rate (in M_sun/year)
		       Default: MDOT=1
		angm - angular momentum of the black hole in geometrized units
		       Default: ANGM=0.998 (i.e. a maximum-rotation Kerr hole)
		alpha - Shakura-Sunyaev viscosity parameter alpha
		       Default: alpha=0.01
		mu   - cosine of the inclination angle (MU=1 would be a face-on disk)
		       Default: mu=0.6

		ii) parameters for setting up annuli:

		tmin - minimum Teff of an annulus which is taken from the grid
		       (all other, cooler, annuli are represented by black-bodies).
		       TMIN should be close to the minimum Teff of the grid, 10,000 K.
		       Default: tmin=9000
		deltar - if set to a non-zero value, the radii are set up equidistant
		       in log R, DELTAR is then delta log R (that is, 
		       log R_{i+1} = log R_i + DELTAR
		       Default: deltar=0.1
		rcut - if set to a non-zero value, it is the cutoff radius of
		       the overall disk
		       Default: rcut=0. (i.e. the cutoff radius is specified through
		       the limiting Teff - see below)
		tlim - if set to a non-zero value, the cutoff radius is specified
		       as that at which Teff(R) = TLIM
		       Default: tlim=1000

		iii) parameters for computing the integrated spectrum

		frmin - minimum frequency for the integrated spectrum in the
		        observer's frame
		        Default: frmin=1.e14
		frmax - maximum frequency for the integrated spectrum in the
		        observer's frame
		        Default: frmax=1.e17
		nfobs - number of frequencies of the integrated spectrum
		        Default: nfobs=300
    '''

	def __init__(self, mass, mdot, angm, mu, alpha=0.01, tlim=1000.0, frmax=1.e17, frmin = 1.e14, 
		         tmin = 9000.0, rout = 0., rcut = None, nfobs = 300):

		'''
		Initialise the model with user parameters
		'''

		self.mass = mass
		self.mdot = mdot
		self.angm = angm
		self.mu = mu
		self.alpha = alpha
		self.tlim = tlim
		self.tmin = tmin
		self.rcut = rcut
		self.frmax = frmax
		self.frmin = frmin
		self.nfobs = nfobs
		self.rout = rout
		self.__name__ = 'Agnspec_Model'

		self.nre = 0
		self.deltar = 0.1


	def get_dispar(self):

		'''
		run the FORTRAN submodule dispar3_sub using the variables defined in the model.
		populates the arrays

		tr 
		mr 
		qr 
		r 
		'''

		self.tr, self.mr, self.qr, self.r, self.nr = \
		        dispar.dispar3 (self.mass, self.mdot, self.angm, 
		        	            self.alpha, self.rout, self.tmin,
		        	            self.deltar, self.nre)

		# truncate arrays
		self.tr = self.tr[:self.nr]
		self.mr = self.mr[:self.nr]
		self.qr = self.qr[:self.nr]
		self.r = self.r[:self.nr]

		if self.rcut == None:
			self.rcut = self.r[self.nr-1] * (self.tlim / 10. ** self.tr[self.nr-1]) ** (-1.333)

		return (0)


# def ringspec(grid, teff,dm0,qgrav,int0):

# 	nt = len(grid.tl)
# 	it = findgen(nt)
# 	ii = interpol(it,grid.tl,teff) 
# 	i=fix(ii)+1
# 	i=i>1<(nt-1)
# 	i1=i-1
# 	a=(teff-tl(i1))/(tl(i)-tl(i1))
# 	a1=1.-a
# 	;
# 	nm=n_elements(ml)
# 	im=findgen(nm)
# 	ii=interpol(im,ml,dm0) 
# 	j=fix(ii)+1
# 	j=j>1<(nm-1)
# 	j1=j-1
# 	b=(dm0-ml(j1))/(ml(j)-ml(j1))
# 	b1=1.-b
# 	;
# 	nq=n_elements(ql)
# 	iq=findgen(nq)
# 	ii=interpol(iq,ql,qgrav) 
# 	k=fix(ii)+1
# 	k=k>1<(nq-1)
# 	k1=k-1
# 	c=(qgrav-ql(k1))/(ql(k)-ql(k1))
# 	c1=1.-c



mass = 1e9
mdot = 5
angm = 0
alpha = 0.01
rout = 0.
tmin = 9000.
deltar=0.1
nre=0 
mu = 0

grid = saved_grid()

# ----------------------
# Basic input parameters
# ----------------------
#
#     MSTAR      - M(star), either in M(Sun), or in grams;
#     MDOT       - M(dot), either in M(Sun)/year; or in g/s
#     AA         - angular momentum of black hole
#     ALPHA      - Shakura-Sunyaev viscosity parameter
#     NRE        - number of radial annuli

model = agnspec_model(mass, mdot, angm, mu)
model.get_dispar()

for i in range(len(model.r)):
	print model.r[i], model.tr[i], model.mr[i], model.qr[i]

#
# -----------------------------------------------
# 1c spectra for individual annuli
# -----------------------------------------------
#
tmpin = np.zeros(20)
totr = np.zeros(len(grid.wls))


