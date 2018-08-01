###############################################################################
#                                                                             #
#                   MeerKAT 1-Million star target list:                       #
#                                get_targets.py                               #
#                                                                             #
#                        Written by Logan A. Pearce (2018)                    #
#                                                                             #
###############################################################################
#
# This script takes in a set of telescope pointings and queries the BL SQL server
# for all objects on the target list falling within the beam during all of the pointings.
# It also checks for any confirmed exoplanet hosts and pulsars in the ATNF catalog which
# might fall within the beam.
#
# Requires:
#   python packages MySQLdb, pandas, astropy, astroquery, datetime
#   a credentials file "~/.my.cnf" with login credentials to BL MySQL server
#
# Input:
#   csv file with pairs of RA/Dec pointings, first column is RA, second column is Dec,
#        in decimal degrees
#
# Output:
#   csv file of all database objects within telescope beam
#   csv file of known exoplanet hosts within beam (if any)
#   csv file of pulsars within beam (if any)
#
# Useage:
# get_targets.py [-h] [-o OUTPUT_FILENAME] [-d ANG_DIAM] [-c CONDITIONS]
#                      [-t TABLE]
#                      filename
#
# positional arguments:
#   filename             the path to a csv file containing tuples of RA,Dec in
#                        decimal degrees

# optional arguments:
#  -h, --help            show this help message and exit
#  -o OUTPUT_FILENAME, --output_filename OUTPUT_FILENAME
#                        csv name for file to output results. default is
#                        targets_datetime.csv ex:
#                        targets-2018-08-01T14:32:59.txt
#  -d ANG_DIAM, --ang_diam ANG_DIAM
#                        angular diameter of the beam in degrees. default=0.8
#                        deg, approx diameter of MeerKAT L-beam
#  -c CONDITIONS, --conditions CONDITIONS
#                        selection criteria. ex: "dist_c<100" returns objects
#                        with distance less than 100 pc (quotes are required
#                        around conditional argument). default=no condition
#  -t TABLE, --table TABLE
#                        name of table to pull results from. default =
#                        1M_target_list
#
# examples:
#    python get_targets.py test.txt   <- Get all objects from the 1M target list that fall within all the
#                            pointings for L-band sized beams
#    python get_targets.py test.txt -c "dist_c < 100"   <- Get all objects out to a distance of 100 pcs
#    python get_targets.py test.txt -c "sptype_c LIKE '%G%'" -d 0.4 -t master_gaia_database   <- Get all SpType G stars
#                            that fall within an S-band sized beam (3 GHz) from the master list of all high quality
#                            Gaia sources


############### Definitions:
def get_targets(ra,dec,ang_diam,conditions,database):
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
    appended_data = [] 
    for r,d,i in zip(ra,dec,index):
        string = 'SELECT * FROM '+str(database)+' \
                    WHERE POWER((ra-('+str(r)+')),2) + POWER((decl - ('+str(d)+')),2) < '+str((ang_diam/2.)**2)+' \
                    AND '+str(conditions)+';'
        dataframe = pd.read_sql(string, con=db)
        # store DataFrame in list
        appended_data.append(dataframe)
        print "I've done ",i+1," of ",len(ra)," total pointings"
    targets = pd.concat(appended_data, axis=0)
    return targets

def find_exoplanets(ra,dec,ang_diam):
    """
    Given an RA/Dec array of any size, queries the NASA Exoplanet Archive and returns a Pandas databse of
    any known exoplanet hosts that fall within a circle of specified angular diameter of the given RA/Decs
    All angles in degrees
    Args:
        ra,dec (array, float [deg]): arrays of pointing coordinates in decimal degrees
        ang_diam (float [deg]): angular size of the diameter of the beam you are simulating.  Default of 0.8 deg
                          is slightly smaller than the MeerKAT beam in L-band to provide a conservative estimate
    Returns:
        Astropy table object of every confirmed exoplanet host within the circle.
    """
    from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive
    from astropy.table import QTable, Table
    # Get exoplanet archive table of all confirmed exoplanets:
    result_table = NasaExoplanetArchive.get_confirmed_planets_table()
    # Convert QTable to regular astropy Table object for easy referencing:
    e = Table(result_table)
    # Convert to pandas df
    cols = e.colnames[0:len(e.colnames)-1]
    ee = e[cols].filled(-9999)
    e = ee.to_pandas()
    
    # Create array for results:
    appended_data = []
    for r,d in zip(ra,dec):
        circle = ((r)-(e['ra']))**2 + ((d)-(e['dec']))**2 <= ang_diam**2
        appended_data.append(e[circle])
    exoplanets = pd.concat(appended_data, axis=0)
    return exoplanets

def find_pulsars(ra,dec,ang_diam):
    """
    Given an RA/Dec array of any size, queries our SQL server table version of the ATNF pulsar catalog.
    Args:
        ra,dec (array, float [deg]): arrays of pointing coordinates in decimal degrees
        ang_diam (float [deg]): angular size of the diameter of the beam you are simulating.  Default of 0.8 deg
                          is slightly smaller than the MeerKAT beam in L-band to provide a conservative estimate
    Returns:
       Pandas dataframe of any pulsars in the ATNF catalog which fall in any of the beams.
    """
    index = range(len(ra))
    appended_data = []
    for r,d,i in zip(ra,dec,index):
        string = 'SELECT * FROM ATNF_psr_cat \
                    WHERE POWER((RAJD-('+str(r)+')),2) + POWER((DECJD - ('+str(d)+')),2) < '+str((ang_diam/2.)**2)+';'
        dataframe = pd.read_sql(string, con=db)
        # store DataFrame in list
        appended_data.append(dataframe)
    pulsars = pd.concat(appended_data, axis=0)
    return pulsars

############### Define arguments:
import argparse
parser = argparse.ArgumentParser()
# Required positional arguments:
parser.add_argument("filename", help="the path to a csv file containing tuples of RA,Dec in decimal degrees")
# Optional positional arguments"
parser.add_argument("-o","--output_filename", help="csv name for file to output results.  default is targets_datetime.csv ex: targets-2018-08-01T14:32:59.txt")
parser.add_argument("-d","--ang_diam", help="angular diameter of the beam in degrees.  default=0.8 deg, approx diameter of MeerKAT L-beam",type=float)
parser.add_argument("-c","--conditions", help='selection criteria.  ex: "dist_c<100" returns objects with distance less than 100 pc \
    (quotes are required around conditional argument).  default=no condition',type=str)
parser.add_argument("-t","--table", help="name of table to pull results from.  default = 1M_target_list",type=str)

args = parser.parse_args()
filename=args.filename

if args.output_filename:
    output_filename=args.output_filename
else:
    import datetime
    output_filename='targets-{date:%Y-%m-%dT%H-%M-%S}.csv'.format( date=datetime.datetime.now() )

if args.ang_diam:
    ad=args.ang_diam
else:
    ad=0.8
if args.conditions:
    conditions = args.conditions
else:
    conditions = '1'
if args.table:
    database=table
else:
    database="1M_target_list"

############### Input pointings:
import pandas as pd
print 'Reading in RA/Dec...'
filein=pd.read_csv(filename)
RA = filein.ix[:,0]
Dec = filein.ix[:,1]

############### Connect to SQL server:
print 'Connecting to SQL server...'
import MySQLdb
db = MySQLdb.connect(host='104.154.94.28',db='loganp',\
                     read_default_file="~/.my.cnf",\
                     autocommit=True,\
                     local_infile = 1)

############### Get all targets in the database within beam pointings:
print 'Retrieving targets...'
targets = get_targets(RA,Dec,ad,conditions,database)
print 'I found ',len(targets),' total objects within those pointings.'

print 'Writing out targets to file...'
targets.to_csv(output_filename,index=False)

############### Search for known exoplanet hosts within the beam:
print 'Looking for confirmed exoplanet hosts...'
exoplanets=find_exoplanets(RA,Dec,ad)
print 'I found ',len(exoplanets),' known exoplanet hosts within these pointings.'
if len(exoplanets)!=0:
    print 'Writing exoplanet hosts to file...'
    exoplanets.to_csv(output_filename.split('.')[0]+'_exoplanet_hosts.csv')
else:
    print 'Zero is not enough to write to a file, sadly...'

############### Search for pulsars with in the beam:
print 'Looking for pulsars...'
pulsars=find_pulsars(RA,Dec,ad)
print 'I found ',len(pulsars),' pulsars within these pointings.'
if len(exoplanets)!=0:
    print 'Writing pulsars to file...'
    exoplanets.to_csv(output_filename.split('.')[0]+'_pulsars.csv')
else:
    print 'Zero is not enough to write to a file, sadly...'

############### Write comment with details of query:
comment='# Query performed on '+str(datetime.datetime.now())+' on input file '+str(filename)+' with beam size '+str(ad)+', with conditions '+str(conditions)+' \
from table '+str(database)+'\n'

with open(output_filename,'r') as contents:
      save = contents.read()
with open(output_filename,'w') as contents:
      contents.write(comment)
with open(output_filename,'a') as contents:
      contents.write(save)

print 'Donezo.'


