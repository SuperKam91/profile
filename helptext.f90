module helptext

  use sz_globals

! Help test rewritten for version 3

  character(len=79), dimension(17) :: ExplainString
  character(len=79), dimension(ncomm,30) :: ExplainStrings
  integer, dimension(ncomm) :: ExplainLength

  data ExplainString(1) /'General notes.'/
  data ExplainString(2) /'Profile is a program initially designed for simulating'/
  data ExplainString(3) /'Sunyaev-Zeldovich anisotropies on the CMB. It has grown considerably'/
  data ExplainString(4) /'and now can be used for more general simulation work with both '/
  data ExplainString(5) /'radio interferometers and (to a limited extent) X-ray telescopes.'/
  data ExplainString(6) /'All commands work on minimum completion, and can be heavily'/
  data ExplainString(7) /'abreviated., eg display-model-xray-profile can be abreviated to'/
  data ExplainString(8) /'d-m-x-p. If the typed command is not unique, profile will prompt for'/
  data ExplainString(9) /'this. The "run" command allows for running of scripts.'/
  data ExplainString(10) /' ' /
  data ExplainString(11) /'Common command types include: '/
  data ExplainString(12) /' "get" for defining various parameters;'/ 
  data ExplainString(13) /'"make" for simulating various different model skies; '/ 
  data ExplainString(14) /'"???-map" which allow operations on the internal map arrays.'/  
  data ExplainString(15) /'"@<command>" will allow execution of shell commands'/
 data ExplainString(16) /'"profile < massprofile.cmd > output.txt" will allow a list of profile'/
 data ExplainString(17) /'commands to be run and the output sent to output.txt (from command line)'/

! 1 'quit'
  data ExplainStrings(1,1)  /'Exit program; will prompt for "are you sure" as well.'/
  data ExplainLength(1) / 1/

! 2 'read-fits-map'
  data ExplainStrings(2,1)  /'Reads a FITS map file into the user selected array If the file doesn''t'/
  data ExplainStrings(2,2)  /'exist, then user is given another go. The pixel size for the FITS map'/
  data ExplainStrings(2,3)  /'and the arrays within Profile (define by get-geometry) do not need to '/
  data ExplainStrings(2,4)  /'be the same and Profile will attempt to interpolate and bin as necessary.' /
  data ExplainStrings(2,5)  /'For this process the user needs to specify whether the FITS map is of' /
  data ExplainStrings(2,6)  /'surface brightness (where Profile will keep the mean the same during ' /
  data ExplainStrings(2,7)  /'conversion) or of flux (where Profile will keep total number of counts ' /
  data ExplainStrings(2,8)  /'the same)'/
  data ExplainLength(2) / 8/

! 3 'display-map'
  data ExplainStrings(3,1)  /'Displays the requested Profile array to the current PGPLOT'/
  data ExplainStrings(3,2)  /'device. Various map manipulations (e.g. zoom, dump to postscript, etc) '/
  data ExplainStrings(3,3)  /'are possible if "graphics displays to be adjustable" is selected'/
  data ExplainStrings(3,4)  /'within graphics-operations.'/
  data ExplainLength(3) / 4/

! 4 'make-aperture'
  data ExplainStrings(4,1)  /'Generates the aperture plane for observations with a radio interferometer'/
  data ExplainStrings(4,2)  /'(chosen by get-telescope) of the current szsky. Write the output to'/
  data ExplainStrings(4,3)  /'arrays (R)eal and (I)maginary parts of aperture plane'/
  data ExplainLength(4) / 3/

! 5 'difference-maps'
  data ExplainStrings(5,1)  /'Differences two user defined maps and writes the output to array '/
  data ExplainStrings(5,2)  /'Di(F)ference map. If desired a goodness of fit statistic is also '/
  data ExplainStrings(5,3)  /'calculated assuming Poisson statistics. The first map should be the '/
  data ExplainStrings(5,4)  /'observed data and the second the model'/
  data ExplainLength(5) / 4/

! 6 'copy-map'
  data ExplainStrings(6,1)  /'Copies one array into another'/
  data ExplainLength(6) / 1/

! 7 'make-cluster'
  data ExplainStrings(7,1)  /'Creates a cluster model as defined in "get-cluster-model" generating the'/
  data ExplainStrings(7,2)  /'following arrays: '/
  data ExplainStrings(7,3)  /'(S)Z sky                 brightness temperature of SZ effect'/
  data ExplainStrings(7,4)  /'Compton (Y) sky          Comptonisation parameter'/
  data ExplainStrings(7,5)  /'(E)lectron density sky   Number of electrons in pixel along los'/
  data ExplainStrings(7,6)  /'(X)-ray sky              Number of X-ray photons expected in pixel'/
  data ExplainLength(7) / 7/

! 8 'help'
  data ExplainStrings(8,1)  /'Listing of available commands'/
  data ExplainLength(8) / 1/

! 9 'get-cluster-model'
  data ExplainStrings(9,1)  /'Allows selection of one of various different cluster models and the '/
  data ExplainStrings(9,2)  /'definition of the required parameters' /
  data ExplainLength(9) / 2/

! 10 'plot-map-profile'
  data ExplainStrings(10,1)  /'Displays an azimuthally smoothed profile of the user selected map' /
  data ExplainLength(10) / 1/

! 11 'histogram-map'
  data ExplainStrings(11,1)  /'Produces a histogram of all the selected map values'/
  data ExplainLength(11) / 1/

! 12 'make-model-xray-map'
  data ExplainStrings(12,1)  /'Convolves the sky map with the PSF of x-ray telescope chosen in '/
  data ExplainStrings(12,2) /'"get-xray-telescope" and adds background level. Allows Poisson noise to be ' /
  data ExplainStrings(12,3) /'added'/
  data ExplainLength(12) / 3/

! 13 'flag-region'
  data ExplainStrings(13,1)  /'Allows a region of a map to be set to zero.'/
  data ExplainLength(13) / 1/

! 14 'get-directory'
  data ExplainStrings(14,1)  /'Allows the user to specify default directories for file input /output.'/
  data ExplainStrings(14,2)  /'These can also be defined before Profile start up by using setenv'/
  data ExplainStrings(14,3)  /'on environment variables "PROFILE_DATA", "PROFILE_TELESCOPE",'/
  data ExplainStrings(14,4)  /'"PROFILE_MODEL", and "PROFILE_POSTSCRIPT"'/
  data ExplainLength(14) / 4/

! 15 'sum-maps'
  data ExplainStrings(15,1)  /'Sums two user defined maps together and allows the result to be written '/
  data ExplainStrings(15,2)  /'to the 2nd map if required'/
  data ExplainLength(15) / 2/

! 16 'Open-plot-device'
  data ExplainStrings(16,1)  /'Opens the plot device. Useful ones are "xwindow" (the default), which is the'/
  data ExplainStrings(16,2)  /'screen when running on the Suns, or "ps" which'/
  data ExplainStrings(16,3)  /'produces a file called "pgplot.ps", which can then be'/
  data ExplainStrings(16,4)  /'printed. '/
  data ExplainLength(16) / 4/

! 17 'close-plot-device'
  data ExplainStrings(17,1)  /'Closes the current PGPLOT device. Do this before opening a new device or ' / 
  data ExplainStrings(17,2)  /'before attempting to print a .ps file'/
  data ExplainLength(17) / 2/

! 18 'rescale-map'
  data ExplainStrings(18,1)  /'Allows a map to be multiplied by a user defined value '/
  data ExplainLength(18) / 1/

! 19 'smooth-map'
  data ExplainStrings(19,1)  /'Convolves a map with a Gaussian'/
  data ExplainLength(19) / 1/

! 20 'plot-model-visibility'
  data ExplainStrings(20,1)  /'Extracts visibilities from the aperture plane made by make-aperture'/
  data ExplainStrings(20,2)  /'for the current east-west interferometer, chosen in get-geometry. The'/
  data ExplainStrings(20,3)  /'routine calculates the u-v positions for a source at a specified'/
  data ExplainStrings(20,4)  /'declination and extracts the apropriate visibility. The routine also'/
  data ExplainStrings(20,5)  /'check for telescope shadowing (this routine is most appropriate for'/
  data ExplainStrings(20,6)  /'the Ryle Telescope).'/
  data ExplainLength(20) / 6/

! 21 'make-source-population'
  data ExplainStrings(21,1)  /'Generates a population of point radio sources and places them on the szsky '/
  data ExplainStrings(21,2)  /'map based on a power law source count or reads in a file with source'/
  data ExplainStrings(21,3)  /'positions, fluxes and spectra.'/
  data ExplainLength(21) / 3/

! 22 'plot-spectra'
  data ExplainStrings(22,1)  /'Plots various microwave spectra (e.g. SZ effect)'/
  data ExplainLength(22) / 1/

! 23 'get-cosmology'
  data ExplainStrings(23,1)  /'Define the cosmology to be used throughout Profile'/
  data ExplainLength(23) / 1/

! 24 'estimate-map-background'
  data ExplainStrings(24,1)  /'Allows an estimate to be made of the background level in a map '/
  data ExplainStrings(24,2)  /'by selecting a region of the map believed to contain only noise.'/
  data ExplainLength(24) / 2/

! 25 - needs work
  data ExplainStrings(25,1)  /'The meat and veg of the fitting stuff. Basically, simulates an xray'/
  data ExplainStrings(25,2)  /'observation (using make-skies [the 2D one]), calculates the likelihood'/
  data ExplainStrings(25,3)  /'of the actual data given the current set of parameters and adjusts '/
  data ExplainStrings(25,4)  /'the parameters until the likelihood is maximised. '/
  data ExplainStrings(25,5)  /'You can fix the beta parameter, but all others are floating. An'/
  data ExplainStrings(25,6)  /'ellipsoidal model needs to be inputted with get-cluster-parameters'/
  data ExplainStrings(25,7)  /'before an ellipsoidal fit is attempted. Fits can take about 700'/
  data ExplainStrings(25,8)  /'iterations, so go and get a coffee. The fitting area can be defined in'/
  data ExplainStrings(25,9)  /'3 ways, by drawing regions, or by fitting every pixel above a certain'/
  data ExplainStrings(25,10)  /'count rate, or both. It is also possible to exclude regions (eg'/
  data ExplainStrings(25,11)  /'cooling flows). If the initial model is very far from the best fit,'/
  data ExplainStrings(25,12)  /'the routine will have problems converging. I therefore suggest that'/
  data ExplainStrings(25,13)  /'the first step should be to find a spherical model in the right ball'/
  data ExplainStrings(25,14)  /'park by hand, and then run bayesian-xray-fit with a fixed'/
  data ExplainStrings(25,15)  /'beta-parameter. This fit can then be improved upon by allowing'/
  data ExplainStrings(25,16)  /'elliptical models and variable betas. Allowing beta to vary increases'/
  data ExplainStrings(25,17)  /'the time for the routine to converge significantly because beta and'/
  data ExplainStrings(25,18)  /'theta depend critically upon each other.'/
  data ExplainLength(25) / 18/

! 26 'graphics-operations'
  data ExplainStrings(26,1)  /'Allows the graphics operations to be turned on and off. When on the'/
  data ExplainStrings(26,2)  /'display can be fiddled. With the mouse over the pgplot window, the'/
  data ExplainStrings(26,3)  /'following keys do the following things. Mouse buttons are labelled'/
  data ExplainStrings(26,4)  /'from the left.  '/
  data ExplainStrings(26,5)  /''/
  data ExplainStrings(26,6)  /' r        reset plot'/
  data ExplainStrings(26,7)  /' q        quit'/
  data ExplainStrings(26,8)  /' g        grey scale'/
  data ExplainStrings(26,9)  /' c        contours'/
  data ExplainStrings(26,10)  /' <SPACE>  centre plot'/
  data ExplainStrings(26,11)  /' Mouse 1  change transfer function'/
  data ExplainStrings(26,12)  /' Mouse 2  zoom in'/
  data ExplainStrings(26,13)  /' Mouse 3  zoom out'/
  data ExplainStrings(26,14)  /''/
  data ExplainStrings(26,15)  /' Contour plotting allows 3 types of contours, linear'/
  data ExplainStrings(26,16)  /'(which autoscales), exponetial (which might not work yet...) and'/
  data ExplainStrings(26,17)  /'defined. For defined contours, the number of contours and the level'/
  data ExplainStrings(26,18)  /'of each contour is then prompted for. '/
  data ExplainLength(26) / 18/

! 27 'make-test-sky'
  data ExplainStrings(27,1)  /'Allows creation of model skies with either a constant level across'/
  data ExplainStrings(27,2)  /'the map, a point source or a gaussian source.'/
  data ExplainLength(27) / 2/

! 28 'make-pointed-observation'
  data ExplainStrings(28,1)  /'Simulates an interferometer observation by sampling from the aperture'/
  data ExplainStrings(28,2)  /'plane at uv positions determined by telescope geomerty and the HA and Dec. '/
  data ExplainStrings(28,3)  /'Gaussian random noise can be added. Visibilities can be written to '/
  data ExplainStrings(28,4)  /'either .vis or .fits files. The visibilities are available to be gridded '/
  data ExplainStrings(28,5)  /'and then mapped within Profile.'/
  data ExplainLength(28) / 5/

! 29 'plot-likelihood'
  data ExplainStrings(29,1)  /'Plots the likelihood function for a fit to an x-ray map varying either'/
  data ExplainStrings(29,2)  /'one or two parameters, holding all the rest constant. NB this is therefore'/
  data ExplainStrings(29,3)  /'*not* formal maginalisation (which would be considerably more tricky).'/
  data ExplainLength(29) / 3/

! 30 'display-sz-sky'
  data ExplainStrings(30,1)  /'Displays the current radio sky (szsky)'/
  data ExplainLength(30) / 1/

! 31 'difference-map-profiles'
  data ExplainStrings(31,1)  /'Produces and differences radial profiles from two user selected maps'/
  data ExplainLength(31) / 1/

! 32 'plot-maximum-amplitude'
    data ExplainStrings(32,1)  /'Displays a mean and maximum amplitude for the visibilities at constant '/
    data ExplainStrings(32,2)  /'radius in the model aperture plane.'/
    data ExplainLength(32) / 2/

! 33 'get-geometry'
    data ExplainStrings(33,1)  /'This routine allows the user to define the size of the maps to be used'/
    data ExplainStrings(33,2)  /'in the program (the map size must be  an integer power of 2 pixels in'/
    data ExplainStrings(33,3)  /'both width and height) and the angular size of these pixels. It allows'/
    data ExplainStrings(33,4)  /'both the pointing centre and the phase centre of radio interferometer'/
    data ExplainStrings(33,5)  /'observations to be displaced from the map centre.'/
    data ExplainLength(33) / 5/

! 34 'clear'
    data ExplainStrings(34,1)  /'Clears the currectly selected PGPLOT device.'/
    data ExplainLength(34) / 1/

! 35 'simple-clean'
    data ExplainStrings(35,1)  /'Applies the CLEAN algorithm to the dirty map generated by '/
    data ExplainStrings(35,2)  /'"map-visibility-data" (making use of the the synthesised beam) '/
    data ExplainStrings(35,3)  /'and generates maps of the restored CLEAN compenents, the residuals '/
    data ExplainStrings(35,4)  /'and their sum, the final CLEAN map. Currently the restoring beam is '/
    data ExplainStrings(35,5)  /'a circular gaussian rather than an elliptical gaussian fitted to the '/
    data ExplainStrings(35,6)  /'inner part of the synthesised beam.'/

! 36 'show-cluster-parameters'
    data ExplainStrings(36,1)  /'Shows the currently selected cluster model, its associated parameters, '/
    data ExplainStrings(36,2)  /'the cosmology, and the selected radio and X-ray telescopes.'/
    data ExplainLength(36) / 2/

! 37 'read-uv-fits'
  data ExplainStrings(37,1)  /'Reads data and uv coverage from a .fits file ready for gridding or '/
  data ExplainStrings(37,2)  /'conversion to .vis format.'/
  data ExplainLength(37) / 2/ 

! 38 'plot-cut'
    data ExplainStrings(38,1)  /'Plots a cut through a user selected map. The map can be zoomed, recentred'/
    data ExplainStrings(38,2)  /'etc before drawing the cut line'/
    data ExplainLength(38) / 2/

! 39 'mcg-xray-fit'
    data ExplainStrings(39,1)  /'Attempts to find a best fitting beta-model to the current x-ray map by'/
    data ExplainStrings(39,2)  /'a method of conjugent gradients. The fit parameter is chi-squared; as'/
    data ExplainStrings(39,3)  /'a result this routine assumes Gaussian statistics and so will only'/
    data ExplainStrings(39,4)  /'work for clusters with good x-ray detections, where poisson statistics'/
    data ExplainStrings(39,5)  /'can be well modelled by a Gaussian.'/
    data ExplainLength(39) / 5/

! 40 'combine-maps'
    data ExplainStrings(40,1)  /'Takes a set of individual FITS maps and their associated noise maps'/
    data ExplainStrings(40,2)  /'and linearly combines them in a weighted fashion to produce one'/
    data ExplainStrings(40,3)  /'or more composite FITS maps. The routine requires a .flist file'/
    data ExplainStrings(40,4)  /'containing the names of the individual maps (stripped of .fits'/
    data ExplainStrings(40,5)  /'extensions). It will then read all the headers of these FITS files'/
    data ExplainStrings(40,6)  /'and produce a .digest file summarising the important information'/
    data ExplainStrings(40,7)  /'to save time on subsequent runs of the same input files. The routine'/
    data ExplainStrings(40,8)  /'will make suggestions as to how many output maps are required in the'/
    data ExplainStrings(40,9)  /'RA and Dec directions and a map centre for the first output map (with '/
    data ExplainStrings(40,10)  /'lowest RA and Dec). Its a good idea to run get-geometry first and '/
    data ExplainStrings(40,11)  /'choose an image size appropriate for the input maps.'/
    data ExplainLength(40) / 11/

! 41 'write-model'
    data ExplainStrings(41,1)  /'Writes to file the current cluster parameters.'/
    data ExplainLength(41) / 1/

! 42 'read-model'
    data ExplainStrings(42,1)  /'Reads from file a set of cluster parameters.'/
    data ExplainLength(42) / 1/

! 43, 44 - the Explain commands themselves

! 45 'stat-map'
    data ExplainStrings(45,1)  /'Find statistics (max, min, mean, rms) of a region of a map'/
    data ExplainLength(45) / 1/

! 46 Not currently in use

! 47 'estimate-source-confusion'
    data ExplainStrings(47,1)  /'OLD TASK - NOT SUPPORTED - USE WITH CARE'/
    data ExplainLength(47) / 1/

! 48 'lensing'
    data ExplainStrings(48,1)  /'Calculates apparent positions of background objects which have been'/
    data ExplainStrings(48,2)  /'lensed by the currently defined cluster. Please refer to MEJ...'/
    data ExplainLength(48) / 2/

! 49 'adaptive-smooth'
    data ExplainStrings(49,1)  /'Attempts to get around the problem of Poisson statistics by adaptively'/
    data ExplainStrings(49,2)  /'smoothing an x-ray map and hoping that Gaussian statistics now'/
    data ExplainStrings(49,3)  /'apply. Can give odd results. I suggest you use bayesian-xray-fit instead.'/
    data ExplainLength(49) / 3/

! 50 'fill-visibilities'
    data ExplainStrings(50,1)  /'Writes out visibilities from model aperture planes to file in .vis'/
    data ExplainStrings(50,2)  /'format using an aperture plane coverage taken from another .vis '/
    data ExplainStrings(50,3)  /'file'/
    data ExplainLength(50) / 3/

! 51 'project-cluster'
    data ExplainStrings(51,1)  /'OLD TASK - NOT SUPPORTED - USE WITH CARE'/
    data ExplainLength(51) / 1/

! 52 'run'
    data ExplainStrings(52,1)  /'Reads in a .cmd file and executes the commands as if entered'/
    data ExplainStrings(52,2)  /'by hand. Lines beginning with # in .cmd file are treated as'/
    data ExplainStrings(52,3)  /'comments'/
    data ExplainLength(52) / 3/

! 53 'read-sz-sky'
    data ExplainStrings(53,1)  /'Reads in a radio sky from various different formats and converts to '/
    data ExplainStrings(53,2)  /'units of surface brightness (Kelvin of brightness temperature i.e.'/
    data ExplainStrings(53,3)  /'making the RJ assumption).'/
    data ExplainStrings(53,4)  /'One common usage of this command is to read in a map produced in'/
    data ExplainStrings(53,5)  /'AIPS from one interferometric radio telescope in order to simulate'/
    data ExplainStrings(53,6)  /'an observation with another. The units of the map are therefore'/
    data ExplainStrings(53,7)  /'likely to be Jy/beam. So to convert this to brightness temperature'/
    data ExplainStrings(53,8)  /'profile will need to know the beam size and therefore asks for'/
    data ExplainStrings(53,9)  /'the bmaj and bmin of the 2D Gaussian that approximates the'/
    data ExplainStrings(53,10) /'synthesised beam - e.g. the CLEAN beam used by IMAGR.'/
    data ExplainLength(53) / 10/

! 54 Not currently in use

! 55 Not currently in use

! 56 Not currently in use

! 57 Not currently in use

! 58 'exit' - see 1

! 59 'multiple-bayessian-xray-fit'
    data ExplainStrings(59,1)  /'Another hacky bit of code by WFG. Loops round the bayessian-xray-fit'/
    data ExplainStrings(59,2)  /'routine x times doing a fit to the Rosat map using different levels of'/
    data ExplainStrings(59,3)  /'count cutoff. Doesn''t require graphics, so can be run as a batch'/
    data ExplainStrings(59,4)  /'process overnight (as each fit will take ~10 minutes, doing 12 will'/
    data ExplainStrings(59,5)  /'take 2 hours). Prints results from each fit at the end of each fit. '/
    data ExplainLength(59) /5/

! 60 Not currently in use

! 61 'monte-carlo-xray-fit'
    data ExplainStrings(61,1)  /'OLD TASK - NOT SUPPORTED - USE WITH CARE'/
    data ExplainLength(61) / 1/

! 62 Not currently in use

! 63 'make-mosaic-observation'
 data ExplainStrings(63,1)  /'Simulates an interferometer observation by sampling from the aperture'/
  data ExplainStrings(63,2)  /'plane at uv positions determined by telescope geomerty and the HA and Dec. '/
  data ExplainStrings(63,3)  /'This process is repeated for each pointing of a mosaiced observation. '/
  data ExplainStrings(63,4)  /'Currently only hexagonal close packed mosaics are supported '/
  data ExplainStrings(63,5)  /'Gaussian random noise can be added. Visibilities can be written to '/
  data ExplainStrings(63,6)  /'either .vis or .fits files. The option for writing to a multi-source'/
  data ExplainStrings(63,7)  /'.fits file exists and almost, but not quite, works... '/
  data ExplainLength(63) / 7/

! 65 'get-telescope'
  data ExplainStrings(65,1)  /'Define a radio interferometer with which to make simulated observations.'/
  data ExplainStrings(65,2)  /'If the telescope is one where the elements are located on the ground '/
  data ExplainStrings(65,3)  /'then Profile will be able to calculate their uv coverages in '/
  data ExplainStrings(65,4)  /'"make-mosaic-observation" and "make-pointed-observation" and these '/
  data ExplainStrings(65,5)  /'telescopes can be written to a .tel file. For comounted or partially '/
  data ExplainStrings(65,6)  /'comounted designs this is not possible. Several are hardcoded into '/
  data ExplainStrings(65,7)  /'Profile and can be used with "fill-visibilities".'/
  data ExplainLength(65) / 7/

! 67 'make-spectral-index-map'
  data ExplainStrings(67,1)  /'Calculates a pixel-by-pixel spectral index for either two input map '/
  data ExplainStrings(67,2)  /'or a best fit between the multiple channels of an szsky.'/
  data ExplainLength(67) / 2/

! 71 'evaluate'
  data ExplainStrings(71,1)  /'Evaluates various terms relating to SZ surface brightness '/
  data ExplainLength(71) / 1/

! 74 'plot-scaling'
! not know about at present...

! 75 'make-cmb-sky'
  data ExplainStrings(75,1)  /'Generates a realisation of the CMB sky based on a Cl power spectrum '/
  data ExplainStrings(75,2)  /'with powers at each l either being drawn from a Gaussian distribution '/
  data ExplainStrings(75,3)  /'or exactly corresponding to the input spectrum (mostly for testing). '/
  data ExplainStrings(75,4)  /'Phases are generated randomly. An option exists for applying a high '/
  data ExplainStrings(75,5)  /'pass filter (with Hanning windowing) to the power spectrum; this is '/
  data ExplainStrings(75,6)  /'useful for reflecting the fact that an interferometer never samples '/
  data ExplainStrings(75,7)  /'the shortest possible baselines since the aperture illumination ' /
  data ExplainStrings(75,8)  /'goes to zero beyond the edge of the dish. By comparison, in the ' / 
  data ExplainStrings(75,9)  /'simulations the primary beam goes to zero beyond the edge of the array ' /
  data ExplainStrings(75,10) /'and/or the beam is approximated as a Gaussian, with the result that the ' /
  data ExplainStrings(75,11) /'corresponding aperture illumination function has infinite extent. ' /
  data ExplainStrings(75,12) /'Therefore a sampled visibility will include a contribution from the ' /
  data ExplainStrings(75,13) /'regions of the aperture plane, and for an unfiltered CMB sky this can ' /
  data ExplainStrings(75,14) /'contain a great deal of power and cause considerabl inaccuracy in the '/
  data ExplainStrings(75,15) /'simulation.' / 
  data ExplainLength(75) / 15/

! 77 'grid-visibilities'
  data ExplainStrings(77,1)  /'Grids sampled visibilities ready for mapping. '/
  data ExplainLength(77) / 1/

! 78 'read-visibilities'
  data ExplainStrings(78,1)  /'Reads data and uv coverage from a .vis file ready for gridding or '/
  data ExplainStrings(78,2)  /'conversion to .fits format.'/
  data ExplainLength(78) / 2/ 

! 79 'write-uv-fits'
  data ExplainStrings(79,1)  /'Writes out sampled visibilities to .fits format.'/
  data ExplainLength(79) / 1/

! 80 'map-visibility-data'
  data ExplainStrings(80,1)  /'Maps the gridded visibility data to produce a dirty map and a synthesised '/
  data ExplainStrings(80,2)  /'beam ready for CLEANing. '/
  data ExplainLength(80) / 2/

! 83 'write-fits-map'
  data ExplainStrings(83,1)  /'Write out any of Profiles arrays in .fits map format. '/
  data ExplainLength(83) / 1/

! 84 'correlate-maps'
  data ExplainStrings(84,1)  /'Calculates the Pearson correlation co-efficient between two input maps'/
  data ExplainLength(84) / 1/

! 89 'get-xray-telescope'
  data ExplainStrings(89,1)  /'Select an X-ray telescope for use in simulating observations in '/
  data ExplainStrings(89,1)  /'"make-model-xray-map".'/
  data ExplainLength(89) / 1/

! 90 'plot-cosmology'
  data ExplainStrings(90,1)  /'Plots the angular distance relation for the current cosmology.'/
  data ExplainLength(90) / 1/

! 91 'display-primary-beam'
  data ExplainStrings(91,1)  /'Displays the currently selected primary beam model.'/
  data ExplainLength(91) / 1/

! 92 'make-drift-observation'
 data ExplainStrings(92,1)  /'Simulates an interferometer drift scan observation by sampling from the'/
  data ExplainStrings(92,2)  /'aperture plane at uv positions determined by telescope geometry and the'/
  data ExplainStrings(92,3)  /'HA and Dec.  This process is repeated for each time sample, and samples'/
  data ExplainStrings(92,4)  /'are mapped to user-defined phase centres.  Gaussian random noise can be'/
  data ExplainStrings(92,5)  /'added. Visibilities can be written to multi-source .fits files.'/
  data ExplainLength(92) / 5/

  end module helptext


