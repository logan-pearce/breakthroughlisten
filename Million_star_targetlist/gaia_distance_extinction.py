###############################################################################
#                                                                             #
#                            Million Star Target List:                        #
#              Computing distance uncertainties from Gaia parallaxes          #
#                                                                             #
#                        Written by Logan A. Pearce (2018)                    #
#                                                                             #
###############################################################################
#
# The script takes in a list of Gaia source galactic lat/long, and distances computed from
# Gaia parallaxes, to determine the extinction in E(B-V) for the source.  The extinction
# is estimated using interpolation from the 3D extinction map StilISM (http://stilism.obspm.fr).
# The interpolation script is adapted from the StilISM group's python scripts found at:
# (https://gitlab.obspm.fr/gplum/stilism/blob/5742f93272a637c0d4bde4d128cc0c72334ff685/web/stilism/stilism.py)
# For the script to generate the input file see the Jupyter Notebook
# "Properties of the Nearest 1M stars".  It runs in parallel process to speed calculation.  
#
# Inputs:
#  - CSV file of original dataframe index, Gaia galactic lat/long, a computed distance
#  - The stilism data cube called stilism_cube.h5 found here: http://stilism.obspm.fr/about
# 
# Output:
#  - csv arrays of each process' result array titled "rank_finalarray_mpi.csv" in the same directory as
#        this script
#
# From the terminal, execute as follows:
#   mpiexec -n xx python gaia_distance_extinction.py  
# Where xx is number of cores.
#
# If running on TACC, write an sbatch script with the following parameters:
#   #SBATCH -N 1
#   #SBATCH -n 48
#   #SBATCH -p normal
#   #SBATCH -t 48:00:00  (<- the max time allowed for Lonestar 5.  Will probably not take that long,
#                          depends on value of nsamples)
#   And in the command line area, write:
#      ibrun python gaia_distance_extinction.py



import numpy as np
import os
import time
import warnings
import csv
from astropy.table import Table
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
import scipy.interpolate as spi
import h5py

warnings.filterwarnings("ignore")

filename = 'input_to_extinction.csv'
hdf5file = 'stilism_cube.h5'

#######################################################################################################
#################################### Definititions ####################################################
#######################################################################################################

def init(hdf5file):
    """Load hdf5, calculate axes values corresponding to data.

    Args:
        hdf5file (str): full path for STILISM HDF5 file.

    Returns:
        dict: headers contains in HDF5 file.
        :func:`np.array`: 3D array which contains E(B-V).
        tuple: (x, y, z) where x,y,z contains array of axes
            corresponding to cube values.
        array: value min for x, y, z axes.
        array: value max for x, y, z axes.

    """
    cube = None
    cubeXerr = None
    cubeYerr = None
    headers = None

    # read hdf5 file
    with h5py.File(hdf5file, 'r') as hf:
        cube = hf['stilism/cube_datas'][:]
        cubeXerr = hf['stilism/cube_err_distance'][:]
        cubeYerr = hf['stilism/cube_err_magnitudemax'][:]
        dc = hf['stilism/cube_datas']

        # Less method call are done with this version:
        headers = {k: v for k, v in dc.attrs.items()}

    sun_position = headers["sun_position"]
    gridstep_values = headers["gridstep_values"]

    # Calculate axes for cube value, with sun at position (0, 0, 0)
    min_axes = -1 * sun_position * gridstep_values
    max_axes = np.abs(min_axes)
    axes = (
        np.linspace(min_axes[0], max_axes[0], cube.shape[0]),
        np.linspace(min_axes[1], max_axes[1], cube.shape[1]),
        np.linspace(min_axes[2], max_axes[2], cube.shape[2])
    )

    # S. Ferron variable for map
    step = np.array(headers["gridstep_values"])
    hw = (np.copy(cube.shape) - 1) / 2.
    points = (
        np.arange(0, cube.shape[0]),
        np.arange(0, cube.shape[1]),
        np.arange(0, cube.shape[2])
    )
    s = hw * step

    return (
        headers,
        cube,
        cubeXerr,
        cubeYerr,
        axes,
        min_axes,
        max_axes,
        step,
        hw,
        points,
        s
    )


def reddening2(vlong, vlat, distance, ulong=u.deg, ulat=u.deg, frame='galactic', step_pc=5):
    """Calculate Reddening versus distance from Sun for one particular distance, lat and long.

    Args:
        vlong (str or double): Longitude value.
        ulong (str): Longitude unit used in :class:`SkyCoord`.
        vlat (str or double): Latitude value.
        ulat (str): Latitude unit used in :class:`SkyCoord`.
        frame (str): Galactic, icrs ... values supported by :class:`SkyCoord`.

    Kwargs:
        step_pc (int): Incremental distance in parsec

    Returns:
        array: Parsec values.
        array: E(B-V) value obtain with integral of linear extrapolation.

    """
    # Calculate the position for 1pc
    sc = SkyCoord(
        vlong,
        vlat,
        distance=1 * u.pc,
        unit=(ulong, ulat),
        frame=frame)
    coords_xyz = sc.transform_to('galactic').represent_as('cartesian').get_xyz().value

    # Find the number of parsec I can calculate before go out the cube
    # (exclude divide by 0)
    not0 = np.where(coords_xyz != 0)
    max_pc = np.amin(
        np.abs(
            np.take(max_axes, not0) / np.take(coords_xyz, not0)))

    # Calculate all coordinates to interpolate (use step_pc)
    distances = np.arange(0, max_pc, step_pc)
    sc = SkyCoord(
        vlong,
        vlat,
        distance=distance,
        unit=(ulong, ulat, 'pc'),
        frame=frame)
    sc = sc.transform_to('galactic').represent_as('cartesian')
    #coords_xyz = np.array([coord.get_xyz().value for coord in sc])
    sc_xyz = sc.get_xyz().value

    # linear interpolation with coordinates
    interpolation = spi.interpn(
        axes,
        cube,
        sc_xyz,
        method='linear')
    xvalue = np.around(distance)
    yvalue = interpolation * step_pc

    # errors
    '''xerrors = spi.interpn(
        axes,
        cubeXErr,
        coords_xyz,
        method='linear')
    yerrors = spi.interpn(
        axes,
        cubeYErr,
        coords_xyz,
        method='linear')'''

    return (xvalue,
        np.around(yvalue, decimals=3))
        #np.around(xerrors, decimals=0),
        #np.around(yerrors, decimals=3))

#######################################################################################################
#################################### Establish parallel process #######################################
#######################################################################################################

from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

# define the communicator
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
ncor = size

################################# Load data and initialize: ##########################################

if rank==0:
    print "Importing data..."

# Initialisation
headers, cube, cubeXErr, cubeYErr, axes, min_axes, max_axes, step, hw, points, s = init(hdf5file)

k = np.loadtxt(filename,skiprows=1,delimiter=',')
index = k[:,0]
b = k[:,1]
l = k[:,2]
dist = k[:,3]


######################### Divy up the data to give each rank some to chew on: #######################
#print index.shape[0]
# How many chunks to divide the data into:
number = float(index.shape[0]/ncor)
#print number

# Create an array called "mylist" which is an array of just a piece of the indicies.  So, rank 0 gets
# indicies index[0*number:1*number], rank 1 gets index[1*number:2*number], etc.
# For ~ 200,000 sources, the number is 4166.6875, so each process gets 4166-4167 sources to chew on.
#print rank,index[number*rank:number*(rank+1)]

# Now each rank pulls out their piece of the puzzle:

b = b[number*rank:number*(rank+1)]
l = l[number*rank:number*(rank+1)]
dist = dist[number*rank:number*(rank+1)]
myindex = index[number*rank:number*(rank+1)]


########### Start the calculation for loop: ################

extinction = np.array([])
count=0
start = time.time()

if rank==0:
    print 'Starting Computation...'

for i in range(len(l))[0:5]:
    print 'Rank ',rank,i
    count=count+1
    x,y = reddening2(l[i], b[i], dist[i])
    extinction = np.append(extinction,y[0])

    mod=count%10000
    if mod==0:
        print 'Rank ',rank,' has done ',count,' loops'
        os.system('touch ext_count_'+str(count)+'.txt.')

        # Write it out.
        t = Table([myindex, dist, b, l, extinction], \
                names=('index', 'dist', 'b','l', 'extinction'))
        t.write(str(rank)+'_finalextinction_mpi.csv', format='csv',overwrite=True)
print extinction

end = time.time()

# Write it out done.
t = Table([myindex, dist, b, l, extinction], \
                  names=('index', 'dist', 'b','l', 'extinction'))
t.write(str(rank)+'_finalextinction_mpi.csv', format='csv',overwrite=True)

end = time.time()
print 'Rank ',rank,' took ',(end - start),' s'

################################ Collect all the arrays into one file ########################################
if rank ==0:
    finalarray = np.loadtxt('0_finalextinction_mpi.csv',delimiter=',',skiprows=1)
    for i in range(1,ncor):
        print i
        a = np.loadtxt(str(i)+'_finalextinction_mpi.csv',delimiter=',',skiprows=1)
        finalarray = np.vstack ([finalarray,a])
    print finalarray.shape
    t = Table([myindex, dist, b, l, extinction], \
                  names=('index', 'dist', 'b','l', 'extinction'))
    t.write('final_extinction_output.csv', format='csv',overwrite=True)

    #os.system('rm *_finalextinction_mpi.csv')
    
