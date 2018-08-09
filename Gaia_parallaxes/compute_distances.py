###############################################################################
#                                                                             #
#                   MeerKAT 1-Million star target list:                       #
#                            compute_distances.py                             #
#                                                                             #
#                        Written by Logan A. Pearce (2018)                    #
#                                                                             #
###############################################################################
#
# This script takes in csv results from batched Gaia archive queries of "download_large_query.py
# and computes the distance
# estimates for each using the methods of Bailer-Jones 2015.  It computes the FWHM of
# the distance estimate distribution, and computes the standard deviation of the distance
# as sd = fwhm/2.355 (thus this estimate is only valid for distance posteriors that are
# well approximated as Gaussian - Gaia sources with fractional parallax uncertainty
# less than 0.1; see Gaia distance estimation jupyter notebook for detailed description).
# It then interpolates some values of Teff for which APSIS was unable to compute a value,
# and estimates a spectral type for all sources with Teff values using the model of
# Eric Mamajek
# http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
# Version 2018.05.24
#
# Requires:
#   python packages pandas, scipy, numpy
#   text file EEM_dwarf_UBVIJHK_colors_Teff.txt from above url
#
# Input:
#   csv files output from download_large_queries.py
#
# Output:
#   the same csv files with appended distance, teff, and sptype columns
#
# Useage:
# compute_distances.py 

from scipy.optimize import brentq
import numpy as np
import os
import time
import warnings
warnings.filterwarnings("ignore")
import pandas as pd

# Make a list of all Gaia job result files:
os.system('ls async*.csv > list')

# Open the relevant files
with open('list') as f:
    z = f.read().splitlines()

s = pd.read_table('EEM_dwarf_UBVIJHK_colors_Teff.txt',comment='#',delim_whitespace=True,index_col=False)

# Teff interpolation coefficients
c5gb = [-0.17129205,0.94192658,-1.93876506,2.03051831,-1.3019204,4.03152045]

count=0
for line in z:
    count=count+1
    print 'Computing number ',count,' of ',len(z)
    # read in csv file:
    k = pd.read_csv(line)
    p = k['parallax'].values
    pe = k['parallax_error'].values
    temp = k['teff_val'].values
    L=1350
    rmax=1e6
    # Initialize new column arrays:
    d,sd=np.zeros(len(k)),np.zeros(len(k))
    fwhm_lo,fwhm_hi = np.zeros(len(k)),np.zeros(len(k))
    tgb = np.zeros(len(k))
    spt = np.zeros(len(k),dtype=object)

    start1=time.time()
    print 'Staring distances'

    for i in range(len(k)):
        ############### compute distance:
        # establish the coefficients of the mode-finding polynomial:
        coeff = np.array([(1./L),(-2),((p[i]/1000.)/((pe[i]/1000.)**2)),-(1./((pe[i]/1000.)**2))])
        # use numpy to find the roots:
        g = np.roots(coeff)
        # Find the number of real roots:
        reals = np.isreal(g)
        realsum = np.sum(reals)
        # If there is one real root, that root is the  mode:
        if realsum == 1:
            gd = np.real(g[np.where(reals)[0]])
        # If all roots are real:
        elif realsum == 3:
            if p[i] >= 0:
                # Take the smallest root:
                gd = np.min(g)
            elif p[i] < 0:
                # Take the positive root (there should be only one):
                gd = g[np.where(g>0)[0]]
        d[i] = gd
    end=time.time()
    print 'That took ',(end-start1)/60/60,' hrs'
    
    print 'Computing SD'
    start=time.time()
    for i in range(len(k)):
        rmode = d[i]
        p1,pe1=p[i]/1000.,pe[i]/1000.
        M = (rmode**2*np.exp(-rmode/L)/pe1)*np.exp((-1./(2*(pe1)**2))*(p1-(1./rmode))**2)
        lo = brentq(lambda x: 2*np.log(x)-(x/L)-(((p1-(1./x))**2)/(2*pe1**2)) \
                    +np.log(2)-np.log(M)-np.log(pe1), 0.001, rmode)
        hi = brentq(lambda x: 2*np.log(x)-(x/L)-(((p1-(1./x))**2)/(2*pe1**2)) \
                    +np.log(2)-np.log(M)-np.log(pe1), rmode, rmax)
        fwhm_lo[i],fwhm_hi[i] = lo,hi

        sd[i] = (hi-lo)/2.355
    end=time.time()
    print 'That took ',(end-start)/60/60,' hrs'

    print 'Estimating Teff where needed, and SpType for all Teffs.'
    start = time.time()
    for i in range(len(k)):
        # Estimate teff's where applicable:
        if np.isnan(temp[i]) and k['bp_g'].values[i]>0.05 and k['bp_g'].values[i]<1.8:
            gb = k['bp_g'].values[i]
            t = c5gb[0]*gb**5 + c5gb[1]*gb**4 + c5gb[2]*gb**3 + c5gb[3]*gb**2 + c5gb[4]*gb + c5gb[5]
            tgb[i] = 10**t
        else:
            pass
        # Estimate spectral type for each source with a teff value
        if np.invert(np.isnan(temp[i])):
            index = np.argmin(np.abs(s['Teff']-temp[i]))
            spt[i] = s['SpT'][index]
            count = count+1
        elif tgb[i] != 0:
            count1=count1+1
            index = np.argmin(np.abs(s['Teff']-tgb[i]))
            spt[i] = s['SpT'][index]
        else:
            count2=count2+1
    end=time.time()
    print 'That took ',(end-start)/60/60,' hrs'

    print 'Writing the result out to file'
    k['dist.c'] = d
    k['fwhm_lo.c'] = fwhm_lo
    k['fwhm_hi.c'] = fwhm_hi
    k['sd.c'] = sd
    k['teff_val.c'] = tgb
    k['sptype.c'] = spt

    k.to_csv(line.strip('.')[0]+'_c.csv',index=False)
    end=time.time()
    print 'The whole file took ',(end-start1)/60/60,' hrs'
    print ''

print 'Donezo'
    
    
    
