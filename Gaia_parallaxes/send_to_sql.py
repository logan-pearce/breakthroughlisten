import os
import MySQLdb

print 'Establishing MySQL connection'
db = MySQLdb.connect(host='104.154.94.28',db='loganp',\
                     read_default_file="~/.my.cnf",\
                     autocommit=True,\
                     local_infile = 1)
cursor=db.cursor()

# Make a list of the output csv files and read in the list
os.system('ls *.csv > list')
with open('list') as f:
    z = f.read().splitlines()

for line in z:
    print 'Loading '+line
    # Remove the first row of strings
    os.system('tail --lines=+2 '+line+' > '+line.split('.')[0]+'_2.csv')
    # Push the database to MySQL
    sqlstring = """LOAD DATA LOCAL INFILE '"""+line.split('.')[0]+"""_2.csv' INTO TABLE master_gaia_database \
        FIELDS TERMINATED BY ',';"""
    cursor.execute(sqlstring)
    print 'I have done ',line
    

'''for line in z:
    print 'Loading '+line
    # Remove the first row of strings
    os.system('tail --lines=+2 '+line+' > '+line.split('.')[0]+'_2.csv')
    # Push the database to MySQL through the "import_csv" bash script
    os.system('./import_csv '+line.split('.')[0]+'_2.csv')
    print 'I have done ',line'''
    
