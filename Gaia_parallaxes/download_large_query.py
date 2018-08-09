###############################################################################
#                                                                             #
#                   MeerKAT 1-Million star target list:                       #
#                         download_large_query.py                             #
#                                                                             #
#                        Written by Logan A. Pearce (2018)                    #
#                                                                             #
###############################################################################
#
# This script queries the Gaia archive in batches of 3M results (the query limit
# without special authorization) for all Gaia targets meeting the data quality filters
# and returns the results as .csv files.  Depending on the Gaia server speed, it takes
# around 12 hours to complete as written.
#
# Requires:
#   python packages numpy, astroquery
#
# Input:
#   none
#
# Output:
#   .csv files containing the results of the queries in 3M row batches ordered
#       by parallax value descending
#
# Useage:
# download_large_query.py 

import numpy as np
from astroquery.gaia import Gaia
import time

start = time.time()

# Initial max parallax value:
para_i = 800
# this is set to 800 because it is larger than the parallax for Proxima Centauri, but smaller than parallaxes
# for the solar system objects in the catalog.


# Max interations: If we get a max of 3M results per query, and we eventually want 32M,
# count_max should be 11 to get 33M (the last one will be truncated)
count_max = 11 

count = 0
while count<count_max:
    start2 = time.time()
    print 'Performing query ',count+1

    querystring = "SELECT source_id, ref_epoch, ra, ra_error, dec, dec_error, parallax, parallax_error, parallax_over_error, parallax_error/parallax AS frac_para_error, pmra, pmra_error, pmdec, pmdec_error, astrometric_n_obs_al, astrometric_chi2_al, astrometric_excess_noise, astrometric_excess_noise_sig, visibility_periods_used, phot_g_n_obs, phot_g_mean_flux, phot_g_mean_flux_error, phot_g_mean_flux_over_error, phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag, bp_rp,bp_G,G_rp, radial_velocity,radial_velocity_error, l,b,ecl_lat,ecl_lon,priam_flags, teff_val,teff_percentile_lower,teff_percentile_upper, a_g_val, a_g_percentile_lower,a_g_percentile_upper, radius_val,radius_percentile_upper,radius_percentile_lower, lum_val,lum_percentile_upper,lum_percentile_lower \
    FROM gaiadr2.gaia_source \
    WHERE parallax < "+str(para_i)+" \
    AND parallax_over_error > 20 \
    AND phot_g_mean_flux_over_error>50 \
    AND phot_rp_mean_flux_over_error>20 \
    AND phot_bp_mean_flux_over_error>20 \
    AND phot_bp_rp_excess_factor < 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2) \
    AND phot_bp_rp_excess_factor > 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2) \
    AND visibility_periods_used>=8 AND astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))\
    ORDER BY parallax DESC"

    job = Gaia.launch_job_async(query=querystring,\
                    verbose=False, dump_to_file=True, output_format='csv')
    c=job.get_results()
    para_i = np.min(c['parallax'])
    count=count+1

    end2 = time.time()
    print 'Took ',(end2 - start2)/60./60.,' hours'
    print ''
    print job
    #jobid = raw_input('Enter Job ID')
    #j = Gaia.remove_jobs(str(jobid))

print 'Done.'
end = time.time()
print 'Took ',(end - start)/60./60.,' hours'


