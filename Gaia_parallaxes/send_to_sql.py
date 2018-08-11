###############################################################################
#                                                                             #
#                   MeerKAT 1-Million star target list:                       #
#                                send_to_sql.py                               #
#                                                                             #
#                        Written by Logan A. Pearce (2018)                    #
#                                                                             #
###############################################################################
#
# This script pushes .csv tables to the Breakthrough Listen MySQL server.
#
# Requires:
#   python packages MySQLdb
#   a credentials file "~/.my.cnf" with login credentials to BL MySQL server
#
# Input:
#   list of csv files you wish to push to MySQL, named "list", which all have
#       a 1-row header.  You must be sure that the tables match the table structure
#       on the SQL server or they will not be loaded correctly.
#


import os
import MySQLdb

print 'Establishing MySQL connection'
db = MySQLdb.connect
(
                     read_default_file="~/.my.cnf",\
                     autocommit=True,\
                     local_infile = 1
)
cursor=db.cursor()

# Connect via MySQL Connector/Python:
# db = mysql.connector.connect(option_files='/Users/loganpearce/.my.cnf')

# Make a list of the output csv files and read in the list
#os.system('ls *_c.csv > list')
with open('list') as f:
    z = f.read().splitlines()

for line in z:
    print 'Loading '+line
    # Remove the first row of strings
    os.system('tail -n+2 '+line+' > '+line.split('.')[0]+'_2.csv')
    # Push the database to MySQL
    sqlstring = """LOAD DATA LOCAL INFILE '"""+line.split('.')[0]+"""_2.csv' INTO TABLE `1M_target_list` \
        FIELDS TERMINATED BY ',';"""
    cursor.execute(sqlstring)
    print 'I have done ',line

os.system('rm list')
