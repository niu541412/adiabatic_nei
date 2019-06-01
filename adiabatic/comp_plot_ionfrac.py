"""
Purpose: make a comparison plot: NEI v.s. CIE.

In this file, create the ionic fraction plot

Written by sw, Jun 02, 2019
"""

# Import used modules
import pyatomdb, pylab
import pickle, os
import numpy as np
import astropy.io.fits as pyfits
from astropy.io import ascii
import matplotlib as mpl

#system parameters
rootpath = os.getcwd()+'/'
# Zlist = [1,2,6,7,8,10,12,14,16,18,20,26,28]

#Read the condition file
confile = rootpath+'adia.exp_phy.info'
conditions = ascii.read(confile)
ncondi = len(conditions)

#1.Compare the ionic fraction========================================
cie_ionfrac = pickle.load(open(rootpath+'cie_case/tionfrac_cie.pkl','rb'))
nei_ionfrac = pickle.load(open(rootpath+'nei_case/tionfrac_nei.pkl','rb'))

# Only compare the ionic fraction of C, N, O, Ne, Mg, and Si because of
# the color table limit
Zlist = [6,7,8,10,12,14]

radius = conditions[0:ncondi]['R']/3.0856780e+18
mpl.style.use('default')
for Z in Zlist:
  elsymb = pyatomdb.atomic.Ztoelsymb(Z)
  fig, ax = pylab.subplots(1, 1)
  fig.show()
  for i in range(5,Z+1):
    ax.plot(radius, nei_ionfrac[Z][i,0:ncondi], color='C'+"%s"%(i-5), \
      linewidth=2, label=elsymb+' '+pyatomdb.atomic.int2roman(i+1))
    ax.plot(radius, cie_ionfrac[Z][i,0:ncondi], color='C'+"%s"%(i-5), \
      linewidth=2.5, linestyle='dashed')
  ax.plot(radius, np.zeros(ncondi), linewidth=2.5, \
    linestyle='dashed', label='CIE')
  ax.set_xlabel('Radius (pc)')
  ax.set_ylabel('Ionic Fraction')
  ax.set_ylim(1e-7,1.5)
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.legend(loc=0)
  pylab.draw()
  fig.savefig(rootpath+'figures/comp_ionfrac/ionfrac_comp_%s.png' %elsymb)
