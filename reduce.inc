************************************************************************
*
*  AMI Project - Reduce program:  reduce.inc
*
*  This file contains working arrays for the raw data analysis program.
*
*  History:
*    13/02/02 - adapted from VSA version [DJT]
*    14/05/04 - revised, frequency data array introduced [DJT]
*     2/06/05 - RA and Dec arrays introduced [djt]
*     3/08/05 - MAX_RACK increased to 16 [djt]
*     8/08/05 - rain gauge data implemented [djt]
*     1/12/05 - use MAX_SAMP1 as dimension for raw() array [djt]
*    29/03/06 - added spacing, sample selections [djt]
*    28/08/06 - added nyquist frequencies to lcal [kg]
*     8/01/07 - added array selection [djt]
*    22/01/07 - added pcal_rain [djt/nhw]
*     2/05/07 - added off_src [djt]
*    15/05/07 - added highpass arrays [MLD]
*     3/07/07 - added configuration, calibration files [djt]
*    14/12/07 - added definition of MAX_CYCLE [djt]
*    18/02/08 - added ideal_lag_sinc array [kjbg1]
*     8/04/08 - added pointing calibration files [djt]
*
************************************************************************

c   Maximum number of antennas

      integer    MAX_ANT
      parameter (MAX_ANT = 10)

c   Maximum number of racks and boards

      integer    MAX_RACK, MAX_BOARD
      parameter (MAX_RACK = 16)
      parameter (MAX_BOARD = 6)

c   Maximum number of baselines

      integer    MAX_BASE
      parameter (MAX_BASE = 45)
c     parameter (MAX_BASE = (MAX_ANT*(MAX_ANT-1))/2)

c   Maximum number of correlated spacings (i.e. 0 and 180-deg boards)

      integer    MAX_SPAC
      parameter (MAX_SPAC = 2*MAX_BASE)

c   Maximum number of lags and frequency channels

      integer    MAX_CHAN, MAX_LAGS, MAX_FREQ
      parameter (MAX_CHAN = 16)
      parameter (MAX_LAGS = MAX_CHAN)
      parameter (MAX_FREQ = MAX_LAGS/2)

c   Total number of data channels per sample

      integer    TOT_CHAN
      parameter (TOT_CHAN = MAX_SPAC*MAX_CHAN)

c   Maximum number of samples

      integer    MAX_SAMP
      integer    MAX_SAMP1
c     parameter (MAX_SAMP = 28800)     ! (8 hours of 1sec samples)
c     parameter (MAX_SAMP = 36000)     ! (10 hours of 1sec samples)
      parameter (MAX_SAMP = 43200)     ! (12 hours of 1sec samples)
c     parameter (MAX_SAMP = 230400)    ! (2 hours of 1/32sec samples)
      parameter (MAX_SAMP1 = MAX_SAMP+1)

c   Maximum number of sources

      integer    MAX_NSRC
      parameter (MAX_NSRC = 16)

c   Maximum length of offset cycle, number of offsets recorded

      integer    MAX_NSTEP
      integer    MAX_NOFFS
      parameter (MAX_NSTEP = 256)
      parameter (MAX_NOFFS = 2*MAX_ANT)

c   Maximum number of offset cycles during an observation

      integer    MAX_CYCLE
      parameter (MAX_CYCLE = 80)

c   Map size for map command
      integer    MAP_SIZE
      parameter (MAP_SIZE = 256)

c   Tenths of seconds per day
      integer    SSDAY
      parameter (SSDAY=864000)

************************************************************************
*
*   Array selection (Small Array or Large Array)
*   Array specific parameters, max antennas, max PC bits

      integer    ami_array
      integer    ami_nant, ami_pcbits
      common /array_pars/ ami_array, ami_nant, ami_pcbits

************************************************************************
*
*   Arrays containing raw data samples
*

c   Correlator sample sequence numbers

      integer*4  seqno(MAX_SAMP)

c   Time stamps for each sample (uts,sts are secs/10)

      integer*4  utod(MAX_SAMP), uts(MAX_SAMP), sts(MAX_SAMP)

c   Sample index within second (0-31)

      integer    smp_insec(MAX_SAMP)

c   Correlator, drive, monitor status words

      integer    corr_stat(MAX_SAMP)
      integer    drv_stat(MAX_SAMP)
c     integer    mon_stat(MAX_SAMP)

c   Index of current pointing centre

      integer    off_src(MAX_SAMP)

c   Offset phase for each sample, offsets applied during each phase

      integer*2  off_phase(MAX_SAMP)
      integer*2  off_values(MAX_NOFFS, MAX_NSTEP)

c   RA and Declination for each raw sample

      real*8     smp_ra(MAX_SAMP)
      real*8     smp_dec(MAX_SAMP)

c   Arrays containing raw correlator data and weights

c     integer*4  corr_data(TOT_CHAN, MAX_SAMP)
c     byte       corr_wt(TOT_CHAN, MAX_SAMP)

c   Rain gauge data

c     byte       rain_data(MAX_ANT, MAX_SAMP)
c     byte       rain_wt(MAX_ANT, MAX_SAMP)
      real*4     rain_data(MAX_ANT, MAX_SAMP)


c   Applied PC values and offsets

      integer*2  pcs(MAX_ANT, MAX_SAMP)
      integer*2  pc_error(MAX_ANT, MAX_SAMP)
      integer*2  pc_offset(MAX_ANT, MAX_SAMP)

c   Recorded AGC values (cos and sin channels for each antenna)

c     integer*2  agcs(2, MAX_ANT, MAX_SAMP)

c   HA and Dec requested pointing, errors and offsets (encoder units)

      integer*4  ha_reqd(MAX_ANT, MAX_SAMP)
      integer*4  dec_reqd(MAX_ANT, MAX_SAMP)
      integer*4  ha_error(MAX_ANT, MAX_SAMP)
      integer*4  dec_error(MAX_ANT, MAX_SAMP)
      integer*4  ha_offset(MAX_ANT, MAX_SAMP)
      integer*4  dec_offset(MAX_ANT, MAX_SAMP)

c   Cryo temperatures and vacuum levels (voltages)

c     integer*2  tcryo(MAX_ANT, MAX_SAMP), vacuum(MAX_ANT, MAX_SAMP)

c   Number of raw samples read

      integer    no_samp

      common /raw_data/ seqno, utod, uts, sts, smp_insec,
c    :                  corr_stat, drv_stat, mon_stat,
     :                  corr_stat, drv_stat,
     :                  off_src, off_phase, smp_ra, smp_dec,
     :                  rain_data,
     :                  pcs, pc_error, pc_offset, off_values,
     :                  ha_reqd, ha_error, ha_offset,
     :                  dec_reqd, dec_error, dec_offset, no_samp


************************************************************************
*
*   Arrays containing reduced data samples
*

c   Raw sample index for each reduced sample

      integer*4  raw(MAX_SAMP1)

c   Time stamps for reduced data samples (secs/10)

      integer*4  utss(MAX_SAMP), stss(MAX_SAMP)

c   Arrays containing lag data, frequency data and weights

      real*4     data(MAX_LAGS, MAX_SPAC, MAX_SAMP)
      complex    zdata(MAX_FREQ, MAX_SPAC, MAX_SAMP)
      real*4     weight(MAX_CHAN, MAX_SPAC, MAX_SAMP)
      
c   Number of reduced data samples

      integer    nsamp


c   Smoothing (integration) applied to reduced data samples
c     nsmooth1 is first stage of integration (before fringe rotation),
c     nsmooth2 is subsequent integration (after fringe rotation),
c     nsmooth is total integration (nsmooth1*nsmooth2).

      integer    nsmooth1, nsmooth2, nsmooth

      common /red_data/ raw, utss, stss,
     :                  data, weight,
     :                  nsamp, nsmooth, nsmooth1, nsmooth2

      equivalence  (data, zdata)

************************************************************************
*
*   Arrays containing simulated data
*

      real*4     sdata(MAX_LAGS, MAX_SPAC, MAX_SAMP)
      complex    zsdata(MAX_FREQ, MAX_SPAC, MAX_SAMP)
      logical    sim_done

      common /simdata/ sdata, sim_done

      equivalence  (sdata, zsdata)

************************************************************************
*
*   Arrays containing PC drift response data
*

      integer    REF_LAG
      parameter (REF_LAG=8)
      integer    MAX_PCDRIFT
      parameter (MAX_PCDRIFT=10001)
      real*4     lag_derr(MAX_PCDRIFT)
      real*4     lag_data(MAX_PCDRIFT, MAX_LAGS, MAX_SPAC)
      complex    lag_sinc(MAX_PCDRIFT, MAX_LAGS, MAX_SPAC)
      complex    ideal_lag_sinc(MAX_PCDRIFT, MAX_LAGS, MAX_SPAC)
      logical    lag_done
      integer    nderr

      common /lagdata/ lag_derr, lag_data, lag_sinc, ideal_lag_sinc,
     :                 lag_done, nderr

************************************************************************
*
*   Arrays containing lag gains (primary calibration)
*

      character*64 lcal_file
      real*4       lag_gain(MAX_LAGS, MAX_SPAC)
      real*4       lag_spac(MAX_LAGS, MAX_SPAC)
      real*4       f_nyq(MAX_SPAC)
      real*4       mean_f_nyq

      common /lcaldata/ lcal_file, lag_gain, lag_spac, f_nyq, mean_f_nyq

************************************************************************
*
*   Arrays containing baseline gains (primary calibration)
*

      character*64 pcal_file
      real*4       pcal_rain(MAX_ANT)
      real*4       pcal_gain(MAX_FREQ,MAX_SPAC)
      real*4       pcal_phase(MAX_FREQ,MAX_SPAC)

      common /pcaldata/ pcal_file, pcal_rain, pcal_gain, pcal_phase

************************************************************************
*
*   Default directories and filenames
*

c   Directory containing calibration gains files

      character*64  gains_dir

c   Directories containing configuration files, raw and reduced data

      character*64  ami_dir, data_dir, etc_dir
      character*64  cmd_dir, log_dir, flag_dir, red_dir
      character*64  fits_dir, vis_dir

c   Directory containing output plot (PostScript) files

      character*64  plot_dir

      common /def_dir/ ami_dir, data_dir, etc_dir, gains_dir, cmd_dir,
     :        log_dir, red_dir, flag_dir, fits_dir, vis_dir, plot_dir

************************************************************************
*
*   Default selections for sample range, spacings, channels, etc
*

c    Strings containing spacing, channel, antenna selections

      character def_spac*128, def_chan*32, def_ant*32, def_skey*2

c    Integers containing default sample range selection

      integer   def_samp1, def_samp2

      common /def_select/ def_spac, def_chan, def_ant, def_skey
      common /def_samp/ def_samp1, def_samp2

************************************************************************
*
*   Calibration and configuration files, defaults in /ami/include/files.inc
*

c   Telescope configuration file

      character*64  tel_file

c   Frequency calibration file

      character*64  freq_file

c   Geometry calibration file

      character*64  geom_file

c   Correlator board zeros file

      character*64  corr_zfile

c   Correlator board lag gains file

      character*64  corr_lfile

c   Correlator baseline gains file

      character*64  corr_pfile

c   PC drift reference data file

      character*64  pcdrift_file

c   PC delay values reference file

      character*64  delays_file

c   5-point analysis data file

      character*64  fivept_file

c   Pointing correction array file

      character*64  point_file

c   Rain gauge nominal values

      character*64  rain_file

c   Calibrator flux archive

      character*64  flux_file

      common /def_files/ tel_file, freq_file, geom_file,
     :                   corr_lfile, corr_pfile, corr_zfile,
     :                   delays_file, rain_file, flux_file,
     :                   pcdrift_file, fivept_file, point_file

************************************************************************
*
*   Sample file selection control
*

c   Filename pattern match for sample files

      character*64  file_match

      common /file_pattern/ file_match

************************************************************************
*
*   Automatic flagging and calibration control
*

c   Control switches for actions when loading a sample file

      logical auto_flag, auto_freq, auto_geom, auto_pcal

      common /auto_pars/ auto_flag, auto_freq, auto_geom, auto_pcal


************************************************************************
*
*   Logging control
*

c   Data reduction log file name, log switch

      character*128 log_file
      logical       logging, applog, filelog, histlog

      common /log_pars/ log_file, logging, applog, filelog, histlog

************************************************************************
*
*   Flagging control
*

c   Arrays containing details of interactive flagging

      integer FL_MAX
      parameter (FL_MAX=256)

      integer      fl_no
      character*4  fl_fkey(FL_MAX)
      character*4  fl_ckey(FL_MAX)
      character*4  fl_skey(FL_MAX)
      character*16 fl_spac(FL_MAX)
      character*16 fl_chan(FL_MAX)
      real*8       fl_st1(FL_MAX), fl_st2(FL_MAX)

      common /flag_table/ fl_fkey, fl_ckey, fl_skey,
     :                    fl_spac, fl_chan, fl_st1, fl_st2, fl_no

c   Percentage of flagged visibilities, for each baseline

      real*4     flagvis(MAX_SPAC)

      common /flag_pars/ flagvis

************************************************************************
*
*   Rain gauge quality flags
*

c   Applied rain gauge correction (amplitudes and weights)
c   Logical array indicating good/bad rain gauge data
c   Logical flag indicating 'update rain' has been run

      real*8     rain_asum, rain_wsum
      logical    rain_ok(MAX_ANT)
      logical    rain_done

      common /rain_pars/ rain_asum, rain_wsum, rain_ok, rain_done

************************************************************************
*
*   Arrays for highpass filter.
*

      real       hpstartf(MAX_SPAC,MAX_CHAN)
      real       hpstopf(MAX_SPAC,MAX_CHAN)
      logical    hpfilt(MAX_SPAC,MAX_CHAN)
      logical    hpass_done

      common /hpass_data/ hpstartf, hpstopf, hpfilt, hpass_done

************************************************************************
*
*   Plot control
*

      real       contrast, brightness
      integer    palette, pflag
      integer    plot_col, line_col, line2_col
      integer    cos_col, sin_col, amp_col, phi_col

c   Graphical output device

      character  plot_device*64

      common /plot_pars/ plot_col, line_col, line2_col,
     :                   cos_col, sin_col, amp_col, phi_col,
     :                   plot_device
      common /image_pars/ palette, pflag, contrast, brightness

************************************************************************
*
*   Working arrays for display, etc
*

      real*4     xx(MAX_SAMP)
      real*4     yy(MAX_SAMP, TOT_CHAN)
      real*4     yy2(MAX_SAMP, TOT_CHAN)
      real*4     yyw(MAX_SAMP, TOT_CHAN)
      real*4     yylag(MAX_SAMP, TOT_CHAN)

      common /work_arrays/ xx, yy, yy2

      equivalence (yy2, yyw, yylag)

