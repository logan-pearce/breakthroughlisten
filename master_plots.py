###############################################################################
#                                                                             #
#                   MeerKAT 1-Million star target list:                       #
#                                master_plots.py                              #
#                                                                             #
#                        Written by Logan A. Pearce (2018)                    #
#                                                                             #
###############################################################################
#
# This script generates png plots of a color-magnitude diagram (CMD) and 2d histogram
# on a Mollweide project for all objects in the BL catalog (subset of all Gaia objects
# which meet data quality criteria).
#
# Requires:
#   python packages MySQLdb (or MySQL Connector/Python), pandas, matplotlib, astropy,
#   a credentials file "~/.my.cnf" with login credentials to BL MySQL server
#
# Input:
#   None
#
# Output:
#   png plot of a CMD of all objects in the SQL database (not corrected for extinction)
#   png plot of a histogram of all objects in the database in Mollweide projection
#
# Useage:
# mater_plots.py [-h]  

import pandas as pd
import matplotlib
# Set the backend to save plots without interacting with them:
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from astropy import units as u
# Import SQL client:
import MySQLdb
# import mysql.connector


print 'Establishing connection'
# Connect via MySQLdb:
db = MySQLdb.connect
(
                     read_default_file="~/.my.cnf",\
                     autocommit=True,\
                     local_infile = 1
)
# Connect via MySQL Connector/Python:
# db = mysql.connector.connect(option_files='/Users/loganpearce/.my.cnf',autocommit=True,allow_local_infile=True)

# Retrieve some columns from all objects in the SQL database:
f = pd.read_sql("SELECT source_id,parallax,bp_rp,phot_g_mean_mag,ra,decl FROM master_gaia_database\
                WHERE 1", con=db)
print "Found ",f.shape,' objects'

# GBP - GRP color:
color = f['bp_rp'].values
# Apparent G band magnitude:
gmag = f['phot_g_mean_mag'].values 

# Convert to absolute magnitude (excluding extinction for now):
GMag = gmag + 5 + 5*np.log10(f['parallax'].values/1000)

print 'Making CMD'

# Turn off interactive plotting to save plots without displaying them:
plt.ioff()
# Make plot figure and axes:
fig, ax = plt.subplots(ncols=1, sharey=True, figsize=(7, 7))
# plot 2d histogram:
hb = ax.hexbin(color,GMag,gridsize=500,cmap='inferno',bins='log',mincnt=1)
# configure plots:
ax.set_title("CMD of all Gaia targets")
ax.set_ylabel('M$_{G}$')
ax.set_xlabel('G$_{BP}$-G$_{RP}$')
ax.invert_yaxis()
cb = fig.colorbar(hb, ax=ax)
cb.set_label('log10(count)')
plt.annotate('{0} stars'.format(color.shape[0]),xy=(0.55,0.85),xycoords='figure fraction')
# write out plot:
plt.savefig('master_cmd.png',format='png')
plt.close(fig)

print 'Making Mollweide'

# Make Sky Coord objects:
pra = coord.Angle(f['ra'].values*u.degree)
pra = pra.wrap_at(180*u.degree)
pdec = coord.Angle(f['decl'].values*u.degree)
# Plot 2d histogram of objects on Mollweide projection:
fig2 = plt.figure(2,figsize=(8,6))
ax = fig2.add_subplot(111, projection="mollweide")
hb = ax.hexbin(pra.radian, pdec.radian,gridsize=500,cmap='inferno',bins='log',mincnt=1)
ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'],color='grey')
ax.set_title('All Gaia targets \n')
cb = fig2.colorbar(hb, ax=ax,orientation="horizontal", pad=0.1)
cb.set_label('log10(count)')
ax.grid(True)
plt.savefig('master_mollweide.png',format='png')
plt.close(fig2)
