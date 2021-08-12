# MOLUSC
MOLUSC (Multi Observational Limits on Unseen Stellar Companions) is a python tool for determining constraints on possible stellar companions to a target star by combining RV, High Resolution Imaging and Gaia data.

## Purpose
The goal of MOLUSC is to generate realistic binary probabilities based on non-detections. To this end we use a Monte Carlo simulation of possible companions to a target star. We compare the generated orbital properties to observational constraints from High Resolution Imaging, Gaia Imaging, RV measurements, and RUWE to determine posteriors and detection limits on the system. 

See Wood et al. (2021) for more detail.

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

The following information is optional, depending on which types of analysis are being used. Examples of RV and High Resolution Imaging text files are provided above. Examples of running the code are shown at the end of the description.

**Radial Velocity**
- File containing RV measurements (in km/s), times (in JD) and errors
- Spectral resolution of RV measurements
- Added Jitter (in m/s)
- RV Sensitivity Floor (in m/s)

**High Resolution Imaging**
- File containing contrast curve, with Delta mag and Separation (in mas)
- Filter in which imaging was taken (2MASS J, H, K, Gaia G, Bp, Rp, and CFHT R, I are available)

**Additional Input**
There are a number of additional options allowing you to alter the generated distribution of companions. You can choose to limit any of the six orbital parameters, $P$, $q$, $e$, $i$, $\omega$, and $\varphi$, to have fixed, minimum or maximum values. The mean and standard deviation of the log-normal period distribution can be altered, as can the exponent for the power-law mass ratio distribution. Finally, you can choose to use either the more conservative 18th mag or more complete 20th mag completeness limits for Gaia imaging limits.

MOLUSC GUI
<img width="1180" alt="MOLUSC GUI" src="https://user-images.githubusercontent.com/64872115/115775317-bb648300-a380-11eb-9f91-2bf03128681a.png">

### Command Line
The command line interface does not offer the full range of user options. Specifically, using the command line interface it is **not** possible to:
- Use multiple contrast curves in one run
- Limit most orbital parameters (it is still possible to limit to only eclipsing companions)
- Adjust the Period or Mass Ratio distributions
- Use the 20th magnitude Gaia completeness limit

Optional arguments:

| Argument | Description |
| --- | ----------- |
| -h, --help | Show help message and exit |
| -v, --verbose | Turn on extra output |
| --age AGE[Gyr] | Input the age of the star in Gyr |
| --jitter JITTER[m/s] | The RV jitter to be added in quadrature \n to the measurement error in m/s |
| -- resolution | The spectral resolution of the RV data |
| --rv_floor RV FLOOR[m/s] |  The lowest RV semi-amplitude which can be rejected in m/s |
| --rv RV_PATH | The path to the file containing the RV data |
| --ao AO_PATH | The path to the file containing the HRI data |
| --filter {J,K,H,G,R,I} | The filter in which the HRI data was taken |
| --ruwe | Apply the RUWE test |
| --gaia | Apply the Gaia contrast test |
| --transit | Generate only transiting/eclipsing companions |
| -a, --all | Write out all generated companions |
                      

## Results
Results from MOLUSC are written out into one or two .csv files. The first file, which is always written out, is output_kept.csv. This includes the generated orbital parameters and calculated parameters (such as projected separation) for all simulated companions which were not ruled out by the analysis.
The second file, output_all.csv, is only written out if you choose to do so. It includes generated and calculated parameters and  rejection status for all simulated companions. For a typical star and 5 million simulated companions these files are often 1GB or more in size, so you may only want to write it out when necessary.

## Plotting
Several different types of plots to explore the posterior distributions can be made from the output files generated by MOLUSC. My plotting code, plotting_tools.py is available above. This includes code for creating corner and detection limit plots.

## Examples

### RV Analysis

GUI:
<img width="1183" alt="image" src="https://user-images.githubusercontent.com/64872115/129258407-cbd4a050-1e49-458c-9af3-36312a1d3188.png">

Command Line:

python BinaryStarGUI.py -v --rv example_rv.txt --resolution 50000 --age 0.015 rv_test 13h50m06.28s -- -40d50m08.9s 100000 1.3


### RV Analysis with Added Jitter

GUI:
<img width="1183" alt="image" src="https://user-images.githubusercontent.com/64872115/129254252-21fa2540-3ebb-47da-95dc-d594335811da.png">

Command Line:

python BinaryStarGUI.py -v --rv example_rv.txt --resolution 50000 --jitter 100 --age 0.015 rv_test 13h50m06.28s -- -40d50m08.9s 100000 1.3


### HRI Analysis

GUI:
<img width="1183" alt="image" src="https://user-images.githubusercontent.com/64872115/129256901-96d0432e-f09c-4f81-b8e4-ba32b933c4b6.png">

Command Line:

python BinaryStarGUI.py -v --ao example_contrast.txt --filter K --age 0.015 hri_test 13h50m06.28s -- -40d50m08.9s 100000 1.3

### RUWE and Gaia Contrast Analysis

GUI:
<img width="1183" alt="image" src="https://user-images.githubusercontent.com/64872115/129257602-f523dd25-b7bd-47d9-b272-9ea9aa3a6f74.png">

Command Line:

python BinaryStarGUI.py -v --age 0.015 --ruwe --gaia gaia_test 13h50m06.28s -- -40d50m08.9s 100000 1.3

### RV Analysis with Transit Limits
This shows how to limit to only transiting companions, which can be done for any of the different types of analysis.

GUI: 
<img width="1183" alt="image" src="https://user-images.githubusercontent.com/64872115/129257972-935804c0-8a0c-4b4e-b428-049b200cb5ed.png">

Command Line:

python BinaryStarGUI.py -v --rv example_rv.txt --resolution 50000 --age 0.015 --transit rv_transit_test 13h50m06.28s -- -40d50m08.9s 100000 1.3

### HRI Analysis with multiple Contrast Curves

GUI:
<img width="1183" alt="image" src="https://user-images.githubusercontent.com/64872115/129259555-d8bfa7a2-2b79-407d-b83d-9cf8b4904671.png">

Command Line:
Multiple Contrast Curves cannot be used via the command line structure. Please use the GUI for this application.
