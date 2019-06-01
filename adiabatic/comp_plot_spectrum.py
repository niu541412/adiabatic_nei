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

#Read the condition file
confile = rootpath+'adia.exp_phy.info'
conditions = ascii.read(confile)
ncondi = len(conditions)

Zlist = [1,2,6,7,8,10,12,14,16,18,20,26,28]
cie_spec = pickle.load(open(rootpath+'cie_case/tcspec_cie.pkl','rb'))
nei_spec = pickle.load(open(rootpath+'nei_case/tcspec_nei.pkl','rb'))
ebins = nei_spec['ebins']
nbins = len(cie_spec[0,:])

nei_tspec = np.zeros([ncondi,nbins], dtype=float)
for Z in Zlist:
  nei_tspec += nei_spec[Z]

condi_index = [26, 100, 125, 200, 250, 283, 342]
for i in condi_index:
  fig, ax = pylab.subplots(1, 1)
  fig.show()
  ax.loglog(ebins, nei_tspec[i,:]*100, drawstyle='steps', label='NEI')
  ax.loglog(ebins, cie_spec[i,:]*100, drawstyle='steps', label='CIE', \
    linestyle='dashed')
  ax.set_xlabel('Energy (keV)')
  ax.set_ylabel('Cts s$^{-1}$ cm$^3$ bin$^{-1}$')
  ax.legend(loc=0)
  ax.set_xlim([0.1,2.0])
  # ax.set_ylim([1e-27,1e-18])
  pylab.draw()
  radius=conditions[i]['R']/3.0856780e+18
  fig.savefig(rootpath+'figures/comp_spec/testmodel_comp_r%4.2fpc.png' % radius)
