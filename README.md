# MOLUSC
MOLUSC (Multi Observational Limits on Unseen Stellar Companions) is a python tool for determining constraints on possible stellar companions to a target star by combining RV, High Resolution Imaging and Gaia data.

## Purpose


## Requirements
This code requires Python 3 and the following packages:
- Numpy
- Scipy
- Astropy
- Astroquery

The files BHAC15_2MASS.txt, BHAC15_CFHT.txt, BHAC15_GAIA.txt, gaia_contrast.txt, and RUWETableGP.txt to be present in the same folder.

## Use
MOLUSC can be used in two ways through the GUI or command line interfaces. The GUI provides a user-friendly interface with which to easily input parameters and read output, and includes slightly more functionality than the command line version, but is not as friendly to large sequences of runs or running on a cluster. The command line option makes it easy to script mutliple runs but is offers slighty fewer user options. 
Regardless of which interface is used the user will need to provide the following information about the code parameters and target star:
- RA and Dec
- Stellar Mass
- Star Age (if not provided a default age of 5 Gyr will be assumed)
- Number of companions to generate
- Output file prefix

The following information is optional, depending on which types of analysis are being used
** Radial Velocity*
- File containing RV measurements, times and errors
- Spectral resolution of RV measurements
- Added Jitter
- RV Sensitivity Floor

** High Resolution Imaging **
- File containing contrast curve
- Filter in which imaging was taken (2MASS J, H, K, Gaia G, and CFHT R, I are available)



## Results


## Plotting
