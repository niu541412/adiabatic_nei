#Purpose: derive the CIE ionic fraction at the condition positions.

# Import used modules
import cie_physderi
import pyatomdb
import pickle, os
from datetime import datetime
from astropy.io import ascii
import multiprocessing as mp

#system parameters
rootpath = os.getcwd()+'/'
Zlist = [1,2,6,7,8,10,12,14,16,18,20,26,28]

#Read the condition file
confile = rootpath+'adia.exp_phy.info'
conditions = ascii.read(confile)
ncondi = len(conditions)
condi_index = range(0,ncondi)

#1.derive the CIE ionic fraction in parallel way-----------
now1 = datetime.now().hour*3600. + datetime.now().minute*60. + \
       datetime.now().second + datetime.now().microsecond/1e6
pool = mp.Pool(mp.cpu_count())
res  = pool.starmap(cie_physderi.deri_cie_ionfrac, \
                    [([Z], conditions, condi_index) for Z in Zlist])
pool.close()
now2 = datetime.now().hour*3600. + datetime.now().minute*60. + \
       datetime.now().second + datetime.now().microsecond/1e6
print("Time Consuming:%7.2f sec." % (now2-now1))
## Time Consuming:  69.84 sec.
#Combine the ionic fraction file
comb_ionfrac = {}
for Z in Zlist:
  zionfrac = pickle.load(open(rootpath+'tciefrac_'+ \
                 pyatomdb.atomic.Ztoelsymb(Z)+'.pkl','rb'))
  comb_ionfrac[Z] = zionfrac[Z]

tmp = open(rootpath+'tionfrac_cie.pkl','wb')
pickle.dump(comb_ionfrac,tmp)
tmp.close()
