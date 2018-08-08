# Breakthrough Listen: Identifying 1 Million Stars for Targeted SETI Observations with MeerKAT

Breakthrough Listen is the largest ever scientific research program aimed at finding evidence of civilizations beyond Earth  (https://breakthroughinitiatives.org/initiative/1).  BL will soon be conducting observations with the MeerKAT radio telescope in South Africa (https://www.ska.ac.za/), in which it will survey 1 million stars for evidence of technosignatures from extraterrestrial intellegence.  The observations are commensal - BL will "piggyback" on MeerKATs large survey programs (LSPs).  While LSPs are conducting their primary science - such as observations of nearby galaxies, pulsar timing, deep extra-galactic surveys - BL will be surveying other nearby stars which fall into the primary beam of telescope.  

This project is to indentify a target list of 1 million stars BL can expect to observe during commensal observations of MeerKAT.  We began with the 1.7 billion objects in the Gaia Data Release 2 catalog, and applied data quality filters to develop a subset of ~32 million high-quality Gaia objects from which to draw our targets.  We then determined the anticipated exact pointings for each LSP to the greatest extent possible, and which of those high-quality Gaia objects would fall within a the beam during those surveys.  We also included a "volume complete" sample - all high-quality Gaia objects out to 160 pc, which will accomodate any unanticipated pointings.  Lastly, we created the script "get_targets.py" which accepts a list of pointings and given parameters and returns all the targets BL could observe in those pointings.  This will be part of generating observing scripts when observations get underway.

## Products
The files within this GitHub repo consist of three general parts:
* Tutorials which demonstrate the need for making data quality filters when working with Gaia DR2 data, and justifying the quality cut decisions we made for this project
* Notebooks and LSP target lists which were used to draw the objects on our final 1.2M object target list.
* Plots and scripts for working with the final 1.2M object target list.

The complete target list .csv file is too large to hosted on GitHub.  It can be found at this public dropbox link: https://www.dropbox.com/s/yklypkckc6m2xx1/1_million_sample_complete.csv?dl=0


## Author
The project was completed by Logan Pearce, under the mentorship of Howard Isaacson, during the Berkeley SETI Research Center internship in summer 2018 (https://seti.berkeley.edu/Internship.html)


## Acknowledgments
This project made use of:
numpy, astropy, pandas, matplotlib, MySQLdb, MySQL Connector/Python

The author wishes to aknowledge to support of the Breakthrough Initiatives and BSRC for the opportunity and funding to come to Berekely and participate in this research; the mentorship of Howard Isaacson, Steve Croft, and Andrew Siemion; the assistance of Griffin Foster, Scott Ransom, and Matt Lebofsky for overcoming numerous sticking points.
