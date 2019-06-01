#Purpose: derive the spectrum

# Import used modules
import nei_physcalc
import pyatomdb
import pickle, os
from datetime import datetime
from astropy.io import ascii
import numpy as np

#system parameters
rootpath = os.getcwd()+'/'
Zlist = [1,2,6,7,8,10,12,14,16,18,20,26,28]

#Read the condition file
confile = rootpath+'adia.exp_phy.info'
conditions = ascii.read(confile)
ncondi = len(conditions)
condi_index = range(0,ncondi)

# set up spectral bins
mineng = 0.1
maxeng = 2.0
nbins = 190
ebins = np.linspace(mineng,maxeng,nbins)

# pre-open the emissivity files
linefile  = \
  pyatomdb.pyfits.open(os.path.expandvars('$ATOMDB/apec_nei_line.fits'))
cocofile  = \
  pyatomdb.pyfits.open(os.path.expandvars('$ATOMDB/apec_nei_comp.fits'))

# spectrum calculation
now1 = datetime.now().hour*3600. + datetime.now().minute*60. + \
       datetime.now().second + datetime.now().microsecond/1e6
for Z in Zlist:
  nei_physcalc.calc_nei_spectrum([Z], condifile=conditions, \
                 condi_index=condi_index, ebins=ebins, \
                 linefile=linefile, cocofile=cocofile, \
                 ionfracfile='tionfrac_'+ \
                   pyatomdb.atomic.Ztoelsymb(Z)+'.pkl')
now2 = datetime.now().hour*3600. + datetime.now().minute*60. + \
       datetime.now().second + datetime.now().microsecond/1e6
print("Time Consuming:%7.2f sec." % (now2-now1))

# Combine the spectrum file
comb_spec = {}
for Z in Zlist:
  zspec = pickle.load(open(rootpath+'tspec_'+ \
                 pyatomdb.atomic.Ztoelsymb(Z)+'.pkl','rb'))
  comb_spec[Z] = zspec[Z]
comb_spec['ebins'] = zspec['ebins']

tmp = open(rootpath+'tcspec_nei.pkl','wb')
pickle.dump(comb_spec,tmp)
tmp.close()

