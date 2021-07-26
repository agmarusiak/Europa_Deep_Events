# Europa_Deep_Events

Files and codes required to investigate deep events from Europa's silicate interior vs its shallow ice shell
Requires sac for matlab to read in sac files (https://home.chpc.utah.edu/~thorne/software.html)
Requires data from Panning's study on Europa's seismicity. (https://github.com/mpanning/EuropaNoise)

Records are in the form:Records/Europa_'ice_thickness'_km_'event_depth'/'distance'_MX'component'.SAC
where ice thickness is either 5 or 20, event depth is either surface or deep, distance is an integer from 1-170 representing epicentral distance and component is E,N, or Z for east, north or vertical

Detection_comp.m will read in records and noise files, plot seismograms with and without added noise using two methods. The first grabs a random section of the noise time series, the second uses a median value of background noise. Currently the preferred noise and catalog 1 of Panning et al 2018 is used. To generate additional noise files use Europa_noise_ppsd with desired noise catalog from Panning (see above for link to GitHub repository).

Records were created using internal strucutre models built from PlanetProfile (https://github.com/vancesteven/PlanetProfile) and then passed to AxiSEM (https://geodynamics.org/cig/software/axisem/) and Instaseis(https://instaseis.net/). Records are in the format of SAC files (http://ds.iris.edu/ds/nodes/dmc/software/downloads/sac/) 
