import os
import MySQLdb

print 'Establishing MySQL connection'
db = MySQLdb.connect(host='104.154.94.28',user='loganp',passwd='bl13579',db='loganp',autocommit=True)
cursor=db.cursor()

# Make a list of the output csv files and read in the list
os.system('ls *.csv > list')
with open('list') as f:
    z = f.read().splitlines()

for line in z:
    print 'Loading '+line
    # Replace all empty cells with NaNs:
    os.system('sh replace_nan '+line+' NaN')
    # Remove the first row of strings
    os.system('tail --lines=+2 '+line+' > '+line.split('.')[0]+'_2.csv')
    # Push the database to MySQL
    sqlstring = """LOAD DATA LOCAL INFILE '"""+line.split('.')[0]+"""_2.csv' INTO TABLE master_gaia_database \
        FIELDS TERMINATED BY ',';"""
    cursor.execute(sqlstring)
    print 'I have done ',line
