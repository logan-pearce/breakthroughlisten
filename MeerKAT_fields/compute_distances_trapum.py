import MySQLdb
import time
from scipy.optimize import brentq
import numpy as np

db = MySQLdb.connect(host='104.154.94.28',db='loganp',\
                     read_default_file="~/.my.cnf",\
                     autocommit=True,\
                     local_infile = 1)
c=db.cursor()
cursor=db.cursor()

L=1350
rmax = 1e6

#######################################################################################################
db = MySQLdb.connect(host='104.154.94.28',user='loganp',passwd='bl13579',db='loganp',autocommit=True,\
                     local_infile = 1)
c=db.cursor()
cursor=db.cursor()

print 'Next is TRAPUM'
start=time.time()
print 'Adding rows to table...'
string = "ALTER TABLE meerkat_gaiadr2_trapum \
            ADD `dist.c` FLOAT NOT NULL,\
            ADD `fwhm_lo.c` FLOAT NOT NULL,\
            ADD `fwhm_hi.c` FLOAT NOT NULL,\
            ADD `sd.c` FLOAT NOT NULL"
c.execute(string)

c.execute("""SELECT source_id,parallax,parallax_error \
                FROM meerkat_gaiadr2_trapum\
                WHERE 1 """)
numrows = c.rowcount

print 'Starting computation! Hang on this will be a while!'

for x in xrange(0,numrows):
    row = c.fetchone()
    # Compute distance:
    coeff = np.array([(1./L),(-2),((row[1]/1000.)/((row[2]/1000.))**2),\
                      -(1./((row[2]/1000.)**2))])
    g = np.roots(coeff)
    # Find the number of real roots:
    reals = np.isreal(g)
    realsum = np.sum(reals)
    # If there is one real root, that root is the  mode:
    if realsum == 1:
        d = np.real(g[np.where(reals)[0]])
    # If all roots are real:
    elif realsum == 3:
        if row[1] >= 0:
            # Take the smallest root:
            d = np.min(g)
        elif row[1] < 0:
            # Take the positive root (there should be only one):
            d = g[np.where(g>0)[0]]
    # Compute FWHM: 
    rmode = d
    M = (rmode**2*np.exp(-rmode/L)/(row[2]/1000.))*\
        np.exp((-1./(2*(row[2]/1000.)**2))*(row[1]/1000.-(1./rmode))**2)
    lo = brentq(lambda x: 2*np.log(x)-(x/L)-((((row[1]/1000.)-(1./x))**2)/(2*(row[2]/1000.)**2)) \
               +np.log(2)-np.log(M)-np.log(row[2]/1000.), 0.001, rmode)
    hi = brentq(lambda x: 2*np.log(x)-(x/L)-((((row[1]/1000.)-(1./x))**2)/(2*(row[2]/1000.)**2)) \
               +np.log(2)-np.log(M)-np.log(row[2]/1000.), rmode, rmax)
    # Compute std dev:
    sigma = (hi-lo)/2.355
    
    string = """UPDATE `meerkat_gaiadr2_trapum` SET `dist.c`="""+str(d[0])+""" \
                    WHERE `source_id` ="""+str(row[0])
    cursor.execute(string)
    string = """UPDATE `meerkat_gaiadr2_trapum` SET `fwhm_lo.c`="""+str(lo)+""" \
                    WHERE `source_id` ="""+str(row[0])
    cursor.execute(string)
    string = """UPDATE `meerkat_gaiadr2_trapum` SET `fwhm_hi.c`="""+str(hi)+""" \
                    WHERE `source_id` ="""+str(row[0])
    cursor.execute(string)
    string = """UPDATE `meerkat_gaiadr2_trapum` SET `sd.c`="""+str(sigma)+""" \
                    WHERE `source_id` ="""+str(row[0])
    cursor.execute(string)
end=time.time()
print 'Whew.  That took ',(end-start)/60/60,' hrs'
