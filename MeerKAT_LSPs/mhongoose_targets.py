import pandas as pd
import numpy as np
import MySQLdb
import os
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

db = MySQLdb.connect(host='104.154.94.28',db='loganp',\
                     read_default_file="~/.my.cnf",\
                     autocommit=True,\
                     local_infile = 1)
c=db.cursor()

def get_targets(ra,dec,ang_diam=0.8,dist=5000,con=db,database='master_gaia_database'):
    """
    Given an RA/Dec array of any size, returns a Pandas databse of
    the targets from the master database that fall within 
    a circle of specified angular diameter of the given RA/Decs
    All angles in degrees
    Args:
        ra,dec (array, float [deg]): arrays of pointing coordinates in decimal degrees
        ang_diam (float [deg]): angular size of the diameter of the beam you are simulating.  Default of 0.8 deg
                          is slightly smaller than the MeerKAT beam in L-band to provide a conservative estimate
        dist (float, [pc]): depth of the desired query in parsecs.  Default is set larger than the largest distance
                          in the master database to return essentially no distance cut
        con (MySQL connection object): the connection you would like to use for the retrieval
        database (str): the database to query.  Default is the master database of all Gaia targets in our program
    Returns:
        Pandas dataframe of every object in the database meeting criteria
    """
    index = range(len(ra))
    appended_data = [] #make a list to store dataframes
    for r,d,i in zip(ra,dec,index):
        string = 'SELECT * FROM '+str(database)+' \
                    WHERE POWER((ra-('+str(r)+')),2) + POWER((decl - ('+str(d)+')),2) < '+str((ang_diam/2.)**2)+' \
                    AND `dist.c` <= '+str(dist)+';'
        dataframe = pd.read_sql(string, con=con)
        # store DataFrame in list
        appended_data.append(dataframe)
        print "I've done ",i+1," of ",len(ra)," total pointings"
    targets = pd.concat(appended_data, axis=0)
    return targets

# Import the mhongoose target list, in dms, ref epoch 2000, coord system ICRS (verified from NED database)
k = pd.read_table('mhongoose_webtable.txt',delim_whitespace=True)

# Convert coords to decimals using astropy sky coordinate objects:
sc = SkyCoord(ra=k['RA'].values, dec=k['Dec'].values, equinox='J2000.0')
k['RA_deg'],k['Dec_deg'] = sc.ra.degree,sc.dec.degree

# Get all targets within each pointing:
targets = get_targets(k['RA_deg'].values[0:5],k['Dec_deg'].values[0:5],ang_diam=0.8)
targets.to_csv('mhongoose.csv')
targets = get_targets(k['RA_deg'].values[5:10],k['Dec_deg'].values[5:10],ang_diam=0.8)
targets.to_csv('mhongoose.csv', mode='a', header=False)
targets = get_targets(k['RA_deg'].values[10:15],k['Dec_deg'].values[10:15],ang_diam=0.8)
targets.to_csv('mhongoose.csv', mode='a', header=False)
targets = get_targets(k['RA_deg'].values[15:20],k['Dec_deg'].values[15:20],ang_diam=0.8)
targets.to_csv('mhongoose.csv', mode='a', header=False)
targets = get_targets(k['RA_deg'].values[20:25],k['Dec_deg'].values[20:25],ang_diam=0.8)
targets.to_csv('mhongoose.csv', mode='a', header=False)
targets = get_targets(k['RA_deg'].values[25:],k['Dec_deg'].values[25:],ang_diam=0.8)
targets.to_csv('mhongoose.csv', mode='a', header=False)
