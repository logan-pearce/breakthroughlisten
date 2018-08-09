###############################################################################
#                                                                             #
#                   MeerKAT 1-Million star target list:                       #
#                           complete_sample_plots.py                          #
#                                                                             #
#                        Written by Logan A. Pearce (2018)                    #
#                                                                             #
###############################################################################
#
# This script generates png histograms of the distribution of distances and Teff's
# (where available) for all objects in the 1M star target list.  Blue = all targets
# from specific pointings, Gold = all objects in the volume complete sample (all objects
# in the BL database out to 160 pc).
#
# Requires:
#   python packages MySQLdb (or MySQL Connector/Python), pandas, matplotlib, astropy,
#   a credentials file "~/.my.cnf" with login credentials to BL MySQL server
#
# Input:
#   csv file "1_million_sample_complete.csv"
#      optionally, get that file from the BL SQL database.
#
# Output:
#   2 png plots of histograms of distance and Teff
#
# Useage:
# complete_sample_plots.py [-h] 

import pandas as pd
import matplotlib
# Set the backend to save plots without interacting with them:
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from astropy import units as u
# Import SQL client:

# Load the csv file of the 1M star target list:
k=pd.read_csv('1_million_sample_complete.csv')

# Or alternately get the list from the BL SQL database (uncomment to use):
# Using MySQLdb
'''
import MySQLdb
db = MySQLdb.connect(host='104.154.94.28',db='loganp',\
                     read_default_file="~/.my.cnf",\
                     autocommit=True,\
                     local_infile = 1)
# Using mysql connector/python:
# import mysql.connector
# db = mysql.connector.connect(option_files='/Users/loganpearce/.my.cnf',host='104.154.94.28',\
#                              database='loganp',autocommit=True,allow_local_infile=True)

string = 'SELECT * FROM `1M_target_list` WHERE 1'
k=pd.read_sql(string,con=db)
''' 

m=k.loc[np.where(k['project']=='Volume Complete')[0]]
n=k.loc[np.where(k['project']!='Volume Complete')[0]]

# Turn off interactive plotting to save plots without displaying them:
plt.ioff()
# Make plot figure and axes:
fig, ax = plt.subplots(ncols=1, figsize=(7, 7))
plt.subplot(211)
plt.hist(n['dist.c'][np.isfinite(n['dist.c'])],bins=100)
plt.xlabel('Distance (pc)')
plt.subplot(212)
plt.hist(n['teff_val'][np.where(n['teff_val']!=0)[0]][np.isfinite(n['teff_val'][np.where(n['teff_val']!=0)[0]])],bins=100)
plt.xlabel('Teff (K)')
plt.tight_layout()
plt.show()
plt.savefig('pointings_dist_teff.png',format='png')
plt.close(fig)

fig, ax = plt.subplots(ncols=1, figsize=(7, 7))
plt.subplot(211)
plt.hist(m['dist.c'][np.isfinite(m['dist.c'])],bins=100,color='goldenrod')
plt.xlabel('Distance (pc)')
plt.subplot(212)
plt.hist(m['teff_val'].values[np.where(m['teff_val'].values>2000.0)[0]],bins=100,color='goldenrod')
plt.xlabel('Teff (K)')
plt.tight_layout()
plt.show()
plt.savefig('vol_compl_dist_teff.png',format='png')
plt.close(fig)
