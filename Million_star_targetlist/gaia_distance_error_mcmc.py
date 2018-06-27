###############################################################################
#                                                                             #
#                            Million Star Target List:                        #
#              Computing distance uncertainties from Gaia parallaxes          #
#                                                                             #
#                        Written by Logan A. Pearce (2018)                    #
#                                                                             #
###############################################################################
#
# The script takes in a list of Gaia source parallaxes, and parallax uncertainties,
# and uses MCMC to determine the 95% CI on the distance estimates by sampling the
# posterior.  For a description of the posterior see the Jupyter Notebook "Computing
# Gaia Distances".  For the script to generate the input file see the Jupyter Notebook
# "Nearest 1M Stars".  It runs in parallel process to speed calculation.  
#
# Inputs:
#  - CSV file of original dataframe index, Gaia parallaxes and uncertainties, a computed distance
# 
# Output:
#  - csv arrays of each process' result array titled "rank_finalarray_mpi.csv" in the same directory as
#        this script
#
# From the terminal, execute as follows:
#   mpiexec -n xx python grand_step3.py path_to_image_file
# Where xx is number of cores and path to image file is the absolute path beginning at the home directory.
#
# If running on TACC, write an sbatch script with the following parameters:
#   #SBATCH -N 1
#   #SBATCH -n 48
#   #SBATCH -p normal
#   #SBATCH -t 48:00:00  (<- the max time allowed for Lonestar 5.  Will probably not take that long,
#                          depends on value of nsamples)
#   And in the command line area, write:
#      ibrun python gaia_distance_error_mcmc.py

#######################################################################################################
#################################### User defined settings: ###########################################
#######################################################################################################
#                                                                                                     #
#       Change the value of these variables to set the MCMC to desired settings:                      #
#           -filename: input file csv                                                                 #
#           -nsamples: number of samples desired from MCMC                                            #
#                                                                                                     #
#######################################################################################################



import numpy as np
import os
import time
import warnings
import csv
from scipy.stats import norm
from scipy.stats.mstats import mquantiles
from astropy.table import Table
import pandas as pd

warnings.filterwarnings("ignore")


filename = 'input_to_mcmc.csv'
nsamples = 7000

#######################################################################################################
#################################### Definititions ####################################################
#######################################################################################################

def prob(r,omega,sigma):
    p = (r**2*np.exp(-r/L)/(sigma))*np.exp((-1/(2*(sigma)**2))*(omega-(1/r))**2)
    if not np.isfinite(p):
        return -np.inf 
    return p

def sampler(lnprob,mu_init=0,omega=0,sigma=0,proposal_width=0.5,nsamples=50,prob=[0.05, 0.5, 0.95]):
    mu_current = mu_init
    posterior = [mu_current]
    yes_accept = 0
    for i in range(nsamples):
        mu_proposal = norm(mu_current,proposal_width).rvs()
        prob_of_proposal = lnprob(mu_proposal,omega,sigma)
        prob_current = lnprob(mu_current,omega,sigma)
        p_accept = prob_of_proposal/prob_current
        dice = np.random.rand()
        accept = dice < p_accept
        if accept:
            mu_current = mu_proposal
            yes_accept = yes_accept+1 #for tracking acceptance rate
        posterior.append(mu_current)
    return posterior,mquantiles(posterior, prob=prob),np.float(yes_accept)/np.float(nsamples)


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

if rank==0:
    print "Importing data..."

k = np.loadtxt(filename,skiprows=1,delimiter=',')
index = k[:,0]
parallax = k[:,1]
para_sigma = k[:,2]
dist = k[:,3]

L = 1350 # scale length in parsecs

######################### Divy up the data to give each rank some to chew on: #######################
#print index.shape[0]
# How many chunks to divide the data into:
number = float(index.shape[0]/ncor)
#print number

# Create an array called "mylist" which is an array of just a piece of the indicies.  So, rank 0 gets
# indicies index[0*number:1*number], rank 1 gets index[1*number:2*number], etc.
# For ~ 200,000 sources, the number is 4166.6875, so each process gets 4166-4167 sources to chew on.
#print rank,index[number*rank:number*(rank+1)]

exec 'mylist = index[%s:%s]' %(number*rank,number*(rank+1))

# Now each rank pulls out their piece of the puzzle:
#parallax = parallax[np.int_(mylist)]
#para_sigma = para_sigma[np.int_(mylist)]
#dist = dist[np.int_(mylist)]

parallax = parallax[number*rank:number*(rank+1)]
para_sigma = para_sigma[number*rank:number*(rank+1)]
dist = dist[number*rank:number*(rank+1)]
myindex = index[number*rank:number*(rank+1)]

'''if rank==0:
    print 'Rank 0 is doing ',myindex
elif rank==1:
    print 'Rank 1 is doing ',myindex
elif rank==2:
    print 'Rank 2 is doing ',myindex
elif rank==3:
    print 'Rank 3 is doing ',myindex

if rank==0:
    print 'Rank 0 is doing ',myindex.shape
elif rank==1:
    print 'Rank 1 is doing ',myindex.shape
elif rank==2:
    print 'Rank 2 is doing ',myindex.shape
elif rank==3:
    print 'Rank 3 is doing ',myindex.shape'''


########### Start the MCMC for loop: ################
quant_median,quant_lo,quant_hi = np.array([]),np.array([]),np.array([])
count = 0
start = time.time()

if rank==0:
    print 'Starting MCMC...'

for i in range(len(mylist))[0:5]:
    count=count+1
    post,quant,accept_rate = sampler(prob,mu_init=dist[i],omega=parallax[i]/1000,sigma=para_sigma[i]/1000,nsamples=nsamples,proposal_width=1.0)
    quant_lo=np.append(quant_lo,quant[0])
    quant_median=np.append(quant_median,quant[1])
    quant_hi=np.append(quant_hi,quant[2])

    mod=count%100
    if mod==0 and rank==0:
        print 'Rank 0 has done ',count,' loops'


# Write it out done.
t = Table([myindex[0:5], dist[0:5], parallax[0:5], para_sigma[0:5], quant_lo, quant_median, quant_hi], \
                  names=('index', 'dist', 'parallax','para_sigma', 'quant_lo', 'quant_median', 'quant_hi'))
t.write(str(rank)+'_finalarray_mpi.csv', format='csv',overwrite=True)

end = time.time()
print 'Rank ',rank,' took ',(end - start),' s'

################################ Collect all the arrays into one file ########################################
if rank ==0:
    finalarray = np.loadtxt('0_finalarray_mpi.csv',delimiter=',',skiprows=1)
    for i in range(1,ncor):
        print i
        a = np.loadtxt(str(i)+'_finalarray_mpi.csv',delimiter=',',skiprows=1)
        finalarray = np.vstack ([finalarray,a])
    print finalarray.shape
    t = Table([finalarray[:,0], finalarray[:,1], finalarray[:,2], finalarray[:,3], finalarray[:,4], finalarray[:,5], finalarray[:,6]], \
                  names=('index', 'dist', 'parallax','para_sigma', 'quant_lo', 'quant_median', 'quant_hi'))
    t.write('final_mcmc_output.csv', format='csv',overwrite=True)

    os.system('rm *_finalarray_mpi.csv')
    
